# tests/testthat/test-user-tdnn-np.R
#
# Behavioral tests for np_tdnn(): shapes, identities, jackknife modes,
# return extras, and input validation.
#
# Interface-equivalence (matrix vs data-frame) lives in the companion file
# test-user-tdnn-np-data-interface.R.

skip_if_not(
  all(vapply(
    c("np_tdnn", "weights_tdnn"),
    exists, logical(1L),
    mode = "function"
  )),
  "np_tdnn() or required helpers not found; skipping file"
)

set.seed(123)

.make_tdnn_data <- function(n = 50L, p = 2L, m = 4L) {
  X  <- matrix(rnorm(n * p), n, p)
  X0 <- matrix(rnorm(m * p), m, p)
  y  <- X[, 1L] - 0.5 * X[, min(2L, p)] + rnorm(n, sd = 0.25)
  list(X = X, X0 = X0, y = y)
}

# ============================================================================
# 1. Basic shape and weighted-mean identity
# ============================================================================

test_that("np_tdnn returns valid shapes and weighted-mean identity", {
  d <- .make_tdnn_data()
  out <- np_tdnn(
    X = d$X, X0 = d$X0, Y = d$y,
    distance = "l2", method = "exact",
    s1 = 2L, s2 = 4L, se = TRUE,
    return = c("W", "ess")
  )

  expect_s3_class(out, "tdnn_np_result")
  expect_s3_class(out, "np_result", exact = FALSE)
  expect_true(is.numeric(out$mu))
  expect_length(out$mu, nrow(d$X0))
  expect_true(is.numeric(out$se))
  expect_length(out$se, nrow(d$X0))
  expect_true(all(is.finite(out$se)))
  expect_equal(dim(out$W), c(nrow(d$X), nrow(d$X0)))
  # W can have negative entries; the identity mu = y'W still holds exactly.
  expect_equal(out$mu, as.numeric(crossprod(d$y, out$W)), tolerance = 1e-12)
  expect_true(all(is.finite(out$ess)))
  expect_identical(out$meta$n, nrow(d$X))
  expect_identical(out$meta$m, nrow(d$X0))
  expect_identical(out$meta$s1, 2L)
  expect_identical(out$meta$s2, 4L)
})

# ============================================================================
# 2. x0 path matches corresponding single-row X0 result
# ============================================================================

test_that("x0 path matches corresponding single-row X0 result", {
  d <- .make_tdnn_data(n = 40L, p = 3L, m = 5L)
  j <- 2L

  out_X0 <- np_tdnn(
    X = d$X, X0 = d$X0, Y = d$y,
    s1 = 2L, s2 = 4L, method = "exact", se = FALSE
  )
  out_x0 <- np_tdnn(
    X = d$X, x0 = d$X0[j, ], Y = d$y,
    s1 = 2L, s2 = 4L, method = "exact", se = FALSE
  )

  expect_equal(out_x0$mu, out_X0$mu[j], tolerance = 1e-12)
})

# ============================================================================
# 3. Jackknife: "none" omits se; constant y gives se = 0
# ============================================================================

test_that("jackknife='none' omits se; constant y gives se = 0", {
  n  <- 30L; p <- 2L; m <- 3L
  X  <- matrix(rnorm(n * p), n, p)
  X0 <- matrix(rnorm(m * p), m, p)
  y_const <- rep(3.5, n)

  out_none <- np_tdnn(
    X = X, X0 = X0, Y = y_const,
    s1 = 2L, s2 = 4L, se = TRUE, jackknife = "none"
  )
  expect_null(out_none$se)

  out_loo <- np_tdnn(
    X = X, X0 = X0, Y = y_const,
    s1 = 2L, s2 = 4L, se = TRUE, jackknife = "loo"
  )
  expect_true(all(abs(out_loo$mu - 3.5) < 1e-12))
  expect_equal(out_loo$se, rep(0, m), tolerance = 1e-10)
})

# ============================================================================
# 4. fast and loo jackknife paths agree (exact method)
# ============================================================================

test_that("fast and loo jackknife paths agree (exact)", {
  d <- .make_tdnn_data(n = 35L, m = 3L)

  out_fast <- np_tdnn(
    X = d$X, X0 = d$X0, Y = d$y,
    s1 = 2L, s2 = 4L, method = "exact", se = TRUE, jackknife = "fast"
  )
  out_loo <- np_tdnn(
    X = d$X, X0 = d$X0, Y = d$y,
    s1 = 2L, s2 = 4L, method = "exact", se = TRUE, jackknife = "loo"
  )

  expect_equal(out_fast$mu, out_loo$mu, tolerance = 1e-12)
  expect_equal(out_fast$se, out_loo$se, tolerance = 1e-12)
})

# ============================================================================
# 5. fast and loo jackknife paths agree (asymp method)
# ============================================================================

test_that("fast and loo jackknife paths agree (asymp)", {
  d <- .make_tdnn_data(n = 38L, m = 3L)

  out_fast <- np_tdnn(
    X = d$X, X0 = d$X0, Y = d$y,
    s1 = 2L, s2 = 4L, method = "asymp", se = TRUE, jackknife = "fast"
  )
  out_loo <- np_tdnn(
    X = d$X, X0 = d$X0, Y = d$y,
    s1 = 2L, s2 = 4L, method = "asymp", se = TRUE, jackknife = "loo"
  )

  expect_equal(out_fast$mu, out_loo$mu, tolerance = 1e-12)
  expect_equal(out_fast$se, out_loo$se, tolerance = 1e-12)
})

# ============================================================================
# 6. Default jackknife mode is "fast"
# ============================================================================

test_that("default jackknife mode matches jackknife='fast' explicitly", {
  d <- .make_tdnn_data(n = 32L, m = 3L)

  out_default <- np_tdnn(
    X = d$X, X0 = d$X0, Y = d$y,
    s1 = 2L, s2 = 4L, method = "exact", se = TRUE
  )
  out_fast <- np_tdnn(
    X = d$X, X0 = d$X0, Y = d$y,
    s1 = 2L, s2 = 4L, method = "exact", se = TRUE, jackknife = "fast"
  )

  expect_equal(out_default$mu, out_fast$mu, tolerance = 1e-12)
  expect_equal(out_default$se, out_fast$se, tolerance = 1e-12)
})

# ============================================================================
# 7. Return extras populated consistently
# ============================================================================

test_that("return extras are populated with correct shapes", {
  d <- .make_tdnn_data(n = 40L, p = 3L, m = 4L)
  out <- np_tdnn(
    X = d$X, X0 = d$X0, Y = d$y,
    distance = "l2", method = "exact",
    s1 = 2L, s2 = 4L, se = FALSE,
    return = c("W", "ranks", "ord", "dx", "y", "ess", "neighbors")
  )

  expect_equal(dim(out$W),     c(nrow(d$X), nrow(d$X0)))
  if (!is.null(out$ranks))
    expect_equal(dim(out$ranks), c(nrow(d$X), nrow(d$X0)))
  if (!is.null(out$ord))
    expect_equal(dim(out$ord),   c(nrow(d$X), nrow(d$X0)))
  if (!is.null(out$dx))
    expect_equal(dim(out$dx),    c(nrow(d$X), nrow(d$X0)))
  expect_equal(out$y,   d$y, tolerance = 1e-12)
  expect_equal(out$ess, 1 / colSums(out$W^2), tolerance = 1e-12)
  # TDNN k = 1: neighbors is a 1 x m matrix of nearest-neighbor indices.
  expect_equal(dim(out$neighbors), c(1L, nrow(d$X0)))
  expect_true(all(out$neighbors >= 1L & out$neighbors <= nrow(d$X)))
})

# ============================================================================
# 8. Signed weights: W can contain negatives, columns sum to 1
# ============================================================================

test_that("TDNN weight columns sum to 1 and may contain negative entries", {
  d   <- .make_tdnn_data(n = 60L, m = 5L)
  out <- np_tdnn(
    X = d$X, X0 = d$X0, Y = d$y,
    s1 = 3L, s2 = 6L, method = "exact", se = FALSE,
    return = "W"
  )

  expect_equal(colSums(out$W), rep(1, ncol(out$W)), tolerance = 1e-12)
  # With these parameters, at least some negative entries should exist.
  expect_true(
    any(out$W < 0),
    info = "TDNN weights should include negative entries (bias correction)"
  )
})

# ============================================================================
# 9. dim defaults to ncol(X) and can be overridden
# ============================================================================

test_that("dim = NULL defaults to ncol(X) and gives same result as explicit dim", {
  d <- .make_tdnn_data(n = 40L, p = 3L, m = 3L)

  out_null <- np_tdnn(
    X = d$X, X0 = d$X0, Y = d$y,
    s1 = 2L, s2 = 4L, method = "exact", se = FALSE, dim = NULL
  )
  out_explicit <- np_tdnn(
    X = d$X, X0 = d$X0, Y = d$y,
    s1 = 2L, s2 = 4L, method = "exact", se = FALSE, dim = 3L
  )

  expect_equal(out_null$mu, out_explicit$mu, tolerance = 1e-12)
  expect_identical(out_null$meta$dim, 3L)
})

# ============================================================================
# 10. Input validation
# ============================================================================

test_that("input validation catches common misuse", {
  d  <- .make_tdnn_data(n = 25L, p = 2L, m = 2L)
  df <- data.frame(X1 = d$X[, 1], X2 = d$X[, 2], Y = d$y)

  # Mutually exclusive X0/x0
  expect_error(
    np_tdnn(X = d$X, X0 = d$X0, x0 = d$X0[1L, ], Y = d$y, s1 = 1L, s2 = 2L),
    "either `X0` or `x0`",
    ignore.case = TRUE
  )
  # Matrix API requires Y
  expect_error(
    np_tdnn(X = d$X, X0 = d$X0, s1 = 1L, s2 = 2L),
    "requires `Y`|`Y` is required",
    ignore.case = TRUE
  )
  # y only supported in data= path
  expect_error(
    np_tdnn(X = d$X, X0 = d$X0, y = d$y, s1 = 1L, s2 = 2L),
    "`y` is only supported with `data=`",
    ignore.case = TRUE
  )
  # data= requires exactly one of y_col/y
  expect_error(
    np_tdnn(data = df, x_cols = c("X1", "X2"), X0 = d$X0, s1 = 1L, s2 = 2L),
    "exactly one of `y_col` or `y`",
    ignore.case = TRUE
  )
  # s1 >= s2
  expect_error(
    np_tdnn(X = d$X, X0 = d$X0, Y = d$y, s1 = 3L, s2 = 2L, se = FALSE),
    "s1 < s2",
    ignore.case = TRUE
  )
  expect_error(
    np_tdnn(X = d$X, X0 = d$X0, Y = d$y, s1 = 3L, s2 = 3L, se = FALSE),
    "s1 < s2",
    ignore.case = TRUE
  )
  # s1 must be a positive integer
  expect_error(
    np_tdnn(X = d$X, X0 = d$X0, Y = d$y, s1 = 0L, s2 = 2L, se = FALSE),
    "positive integer",
    ignore.case = TRUE
  )
  expect_error(
    np_tdnn(X = d$X, X0 = d$X0, Y = d$y, s1 = 1.5, s2 = 3L, se = FALSE),
    "integer scalar",
    ignore.case = TRUE
  )
  # Jackknife feasibility: s2 >= n
  expect_error(
    np_tdnn(
      X = d$X, X0 = d$X0, Y = d$y,
      s1 = 1L, s2 = nrow(d$X), se = TRUE, jackknife = "loo"
    ),
    "s2 < n",
    ignore.case = TRUE
  )
  # data= requires x_cols
  expect_error(
    np_tdnn(data = df, y_col = "Y", X0 = d$X0, s1 = 1L, s2 = 2L),
    "x_cols",
    ignore.case = TRUE
  )
  # data= must be a data.frame
  expect_error(
    np_tdnn(
      data = as.matrix(df), x_cols = c("X1", "X2"), y_col = "Y",
      X0 = d$X0, s1 = 1L, s2 = 2L
    ),
    "data.frame",
    ignore.case = TRUE
  )
})

# ============================================================================
# 11. Warns on unknown return values
# ============================================================================

test_that("np_tdnn warns on unknown return values", {
  d <- .make_tdnn_data(n = 25L, m = 2L)
  expect_warning(
    np_tdnn(
      X = d$X, X0 = d$X0, Y = d$y,
      s1 = 1L, s2 = 2L, se = FALSE, return = c("W", "typo")
    ),
    "Unknown `return` value"
  )
  expect_no_warning(
    np_tdnn(
      X = d$X, X0 = d$X0, Y = d$y,
      s1 = 1L, s2 = 2L, se = FALSE,
      return = c("W", "ranks", "ord", "dx", "y", "ess", "neighbors")
    )
  )
})

test_that("np_tdnn jackknife feasibility diagnostics include expected and actual values", {
  d <- .make_tdnn_data(n = 25L, p = 2L, m = 2L)
  msg <- tryCatch(
    {
      np_tdnn(
        X = d$X, X0 = d$X0, Y = d$y,
        s1 = 1L, s2 = nrow(d$X), se = TRUE, jackknife = "loo"
      )
      ""
    },
    error = function(e) conditionMessage(e)
  )
  expect_match(msg, "Expected:", ignore.case = FALSE)
  expect_match(msg, "Got:", ignore.case = FALSE)
})

# ============================================================================
# 12. print method runs without error and shows expected fields
# ============================================================================

test_that("print.tdnn_np_result runs and shows key fields", {
  d   <- .make_tdnn_data(n = 20L, m = 2L)
  out <- np_tdnn(X = d$X, X0 = d$X0, Y = d$y, s1 = 1L, s2 = 2L, se = FALSE)
  txt <- capture.output(print(out))
  expect_true(any(grepl("<tdnn_np_result>", txt, fixed = TRUE)))
  expect_true(any(grepl("TDNN", txt)))
  expect_true(any(grepl("s1=", txt)))
  expect_true(any(grepl("s2=", txt)))
})
