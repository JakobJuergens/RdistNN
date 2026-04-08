skip_if_not(
  all(vapply(
    c("np_dnn", "weights_dnn"),
    exists, logical(1),
    mode = "function"
  )),
  "np_dnn() or required helpers not found; skipping file"
)

set.seed(123)

.make_np_data <- function(n = 40L, p = 2L, m = 4L) {
  X <- matrix(rnorm(n * p), n, p)
  X0 <- matrix(rnorm(m * p), m, p)
  y <- X[, 1] - 0.5 * X[, min(2L, p)] + rnorm(n, sd = 0.25)
  list(X = X, X0 = X0, y = y)
}

test_that("np_dnn basic matrix API returns valid shapes and weighted mean identity", {
  d <- .make_np_data()
  out <- np_dnn(
    X = d$X, X0 = d$X0, Y = d$y,
    distance = "l2", method = "exact",
    k = 2L, s = 3L, se = TRUE,
    return = c("W", "ess")
  )

  expect_s3_class(out, "np_result")
  expect_true(is.numeric(out$mu))
  expect_length(out$mu, nrow(d$X0))
  expect_true(is.numeric(out$se))
  expect_length(out$se, nrow(d$X0))
  expect_true(all(is.finite(out$se)))
  expect_equal(dim(out$W), c(nrow(d$X), nrow(d$X0)))
  expect_equal(out$mu, as.numeric(crossprod(d$y, out$W)), tolerance = 1e-12)
  expect_true(all(is.finite(out$ess)))
  expect_identical(out$meta$n, nrow(d$X))
  expect_identical(out$meta$m, nrow(d$X0))
})

test_that("x0 path matches corresponding single-row X0 result", {
  d <- .make_np_data(n = 35L, p = 3L, m = 5L)
  j <- 3L

  out_X0 <- np_dnn(
    X = d$X, X0 = d$X0, Y = d$y,
    k = 2L, s = 4L, method = "exact", se = FALSE
  )
  out_x0 <- np_dnn(
    X = d$X, x0 = d$X0[j, ], Y = d$y,
    k = 2L, s = 4L, method = "exact", se = FALSE
  )

  expect_equal(out_x0$mu, out_X0$mu[j], tolerance = 1e-12)
})

test_that("data-frame API supports y_col and external y paths", {
  d <- .make_np_data(n = 30L, p = 2L, m = 3L)
  df <- data.frame(X1 = d$X[, 1], X2 = d$X[, 2], Y = d$y)

  out_col <- np_dnn(
    data = df, x_cols = c("X1", "X2"), y_col = "Y",
    X0 = d$X0, k = 2L, s = 3L, se = FALSE
  )
  out_vec <- np_dnn(
    data = df, x_cols = 1:2, y = d$y,
    X0 = d$X0, k = 2L, s = 3L, se = FALSE
  )
  out_col_idx <- np_dnn(
    data = df, x_cols = c("X1", "X2"), y_col = 3L,
    X0 = d$X0, k = 2L, s = 3L, se = FALSE
  )

  expect_equal(out_col$mu, out_vec$mu, tolerance = 1e-12)
  expect_equal(out_col$mu, out_col_idx$mu, tolerance = 1e-12)
})

test_that("jackknife behavior: none omits se; constant outcomes give zero se", {
  n <- 24L
  p <- 2L
  m <- 3L
  X <- matrix(rnorm(n * p), n, p)
  X0 <- matrix(rnorm(m * p), m, p)
  y_const <- rep(4.25, n)

  out_none <- np_dnn(
    X = X, X0 = X0, Y = y_const,
    k = 2L, s = 4L, se = TRUE, jackknife = "none"
  )
  expect_null(out_none$se)

  out_loo <- np_dnn(
    X = X, X0 = X0, Y = y_const,
    k = 2L, s = 4L, se = TRUE, jackknife = "loo"
  )
  expect_true(all(abs(out_loo$mu - 4.25) < 1e-12))
  expect_equal(out_loo$se, rep(0, m), tolerance = 1e-10)
})

test_that("default jackknife mode uses fast path and matches loo reference", {
  d <- .make_np_data(n = 28L, p = 2L, m = 3L)

  out_default <- np_dnn(
    X = d$X, X0 = d$X0, Y = d$y,
    k = 2L, s = 4L, method = "exact", se = TRUE
  )
  out_fast <- np_dnn(
    X = d$X, X0 = d$X0, Y = d$y,
    k = 2L, s = 4L, method = "exact", se = TRUE, jackknife = "fast"
  )
  out_loo <- np_dnn(
    X = d$X, X0 = d$X0, Y = d$y,
    k = 2L, s = 4L, method = "exact", se = TRUE, jackknife = "loo"
  )

  expect_equal(out_default$mu, out_fast$mu, tolerance = 1e-12)
  expect_equal(out_default$se, out_fast$se, tolerance = 1e-12)
  expect_equal(out_fast$se, out_loo$se, tolerance = 1e-12)
})

test_that("asymp method fast and loo jackknife paths agree", {
  d <- .make_np_data(n = 32L, p = 2L, m = 3L)

  out_fast <- np_dnn(
    X = d$X, X0 = d$X0, Y = d$y,
    method = "asymp",
    k = 3L, s = 4L, se = TRUE, jackknife = "fast"
  )
  out_loo <- np_dnn(
    X = d$X, X0 = d$X0, Y = d$y,
    method = "asymp",
    k = 3L, s = 4L, se = TRUE, jackknife = "loo"
  )

  expect_equal(out_fast$mu, out_loo$mu, tolerance = 1e-12)
  expect_equal(out_fast$se, out_loo$se, tolerance = 1e-12)
})

test_that("return extras are populated consistently", {
  d <- .make_np_data(n = 36L, p = 3L, m = 4L)
  out <- np_dnn(
    X = d$X, X0 = d$X0, Y = d$y,
    distance = "l2", method = "exact",
    k = 3L, s = 4L, se = FALSE,
    return = c("W", "ranks", "ord", "dx", "y", "ess", "neighbors")
  )

  expect_equal(dim(out$W), c(nrow(d$X), nrow(d$X0)))
  if (!is.null(out$ranks)) expect_equal(dim(out$ranks), c(nrow(d$X), nrow(d$X0)))
  if (!is.null(out$ord)) expect_equal(dim(out$ord), c(nrow(d$X), nrow(d$X0)))
  if (!is.null(out$dx)) expect_equal(dim(out$dx), c(nrow(d$X), nrow(d$X0)))
  expect_equal(out$y, d$y, tolerance = 1e-12)
  expect_equal(out$ess, 1 / colSums(out$W^2), tolerance = 1e-12)
  expect_equal(dim(out$neighbors), c(3L, nrow(d$X0)))
  expect_true(all(out$neighbors >= 1 & out$neighbors <= nrow(d$X)))
})

test_that("input validation catches common misuse", {
  d <- .make_np_data(n = 20L, p = 2L, m = 2L)
  df <- data.frame(X1 = d$X[, 1], X2 = d$X[, 2], Y = d$y)

  expect_error(
    np_dnn(X = d$X, X0 = d$X0, x0 = d$X0[1, ], Y = d$y, k = 1L, s = 1L),
    "either `X0` or `x0`",
    ignore.case = TRUE
  )
  expect_error(
    np_dnn(X = d$X, X0 = d$X0, k = 1L, s = 1L),
    "requires `Y`",
    ignore.case = TRUE
  )
  expect_error(
    np_dnn(X = d$X, X0 = d$X0, y = d$y, k = 1L, s = 1L),
    "`y` is only supported with `data=`",
    ignore.case = TRUE
  )
  expect_error(
    np_dnn(data = df, x_cols = c("X1", "X2"), X0 = d$X0, k = 1L, s = 1L),
    "exactly one of `y_col` or `y`",
    ignore.case = TRUE
  )
  expect_error(
    np_dnn(data = df, x_cols = c("X1", "X2"), y_col = "Y", y = d$y, X0 = d$X0, k = 1L, s = 1L),
    "exactly one of `y_col` or `y`",
    ignore.case = TRUE
  )
  expect_error(
    np_dnn(X = d$X, X0 = d$X0, Y = d$y, k = 1L, s = nrow(d$X), se = TRUE, jackknife = "loo"),
    "s < n",
    ignore.case = TRUE
  )
  expect_error(
    np_dnn(data = as.matrix(df), x_cols = c("X1", "X2"), y_col = "Y", X0 = d$X0),
    "data.frame",
    ignore.case = TRUE
  )
  expect_error(
    np_dnn(data = df, y_col = "Y", X0 = d$X0),
    "x_cols",
    ignore.case = TRUE
  )
  expect_error(
    np_dnn(data = df, x_cols = c("X1", "X_missing"), y_col = "Y", X0 = d$X0, k = 1L, s = 2L, se = FALSE),
    "specified by `x_cols`",
    ignore.case = TRUE
  )
  df_bad_x <- df
  df_bad_x$X1 <- as.character(df_bad_x$X1)
  expect_error(
    np_dnn(data = df_bad_x, x_cols = c("X1", "X2"), y_col = "Y", X0 = d$X0, k = 1L, s = 2L, se = FALSE),
    "X.*numeric",
    ignore.case = TRUE
  )
  # k and s validation at entry point
  expect_error(
    np_dnn(X = d$X, X0 = d$X0, Y = d$y, k = 0L, s = 3L, se = FALSE),
    "positive integer"
  )
  expect_error(
    np_dnn(X = d$X, X0 = d$X0, Y = d$y, k = 3L, s = 2L, se = FALSE),
    "k <= s"
  )
  expect_error(
    np_dnn(X = d$X, X0 = d$X0, Y = d$y, k = 1.5, s = 3L, se = FALSE),
    "integer scalar"
  )
  expect_error(
    np_dnn(X = d$X, X0 = d$X0, Y = d$y, k = "foo", s = 3L, se = FALSE),
    "integer scalar"
  )
})

test_that("np_dnn warns on unknown return values", {
  d <- .make_np_data(n = 20L, p = 2L, m = 2L)
  expect_warning(
    np_dnn(X = d$X, X0 = d$X0, Y = d$y, k = 1L, s = 2L, se = FALSE, return = c("W", "typo")),
    "Unknown `return` value"
  )
  expect_no_warning(
    np_dnn(X = d$X, X0 = d$X0, Y = d$y, k = 1L, s = 2L, se = FALSE,
           return = c("W", "ranks", "ord", "dx", "y", "ess", "neighbors"))
  )
})

test_that("print.np_result runs without error", {
  d <- .make_np_data(n = 16L, p = 2L, m = 2L)
  out <- np_dnn(
    X = d$X, X0 = d$X0, Y = d$y,
    k = 1L, s = 2L, se = FALSE
  )
  txt <- capture.output(print(out))
  expect_true(any(grepl("<np_result>", txt, fixed = TRUE)))
})
