# tests/testthat/test-user-np-dnn-data-interface.R
#
# Interface-equivalence tests for np_dnn().
#
# These tests verify three things:
#
#   1. The data-frame interface and the matrix interface produce numerically
#      identical estimates (mu) and standard errors (se) for the same problem.
#
#   2. All sub-variants of the data-frame interface (column names, integer
#      indices, external y vector) agree with each other and with the matrix
#      interface.
#
#   3. Input pathologies that are caught at the data-frame boundary produce
#      precise, distinguishable diagnostics.
#
# Design note ref: note 13 §10.2 ("Add interface-equivalence tests").

skip_if_not(
  all(vapply(c("np_dnn"), exists, logical(1L), mode = "function")),
  "np_dnn() not found; skipping file"
)

# ---- shared helper -----------------------------------------------------------

# Produces a deterministic dataset with a named data.frame suitable for both
# the matrix interface (X, X0, Y) and the data-frame interface (df, x_cols).
.make_np_iface_data <- function(n = 50L, p = 2L, m = 4L, seed = 42L) {
  set.seed(seed)
  X  <- matrix(rnorm(n * p), n, p)
  X0 <- matrix(rnorm(m * p), m, p)
  y  <- X[, 1L] - 0.5 * X[, min(2L, p)] + rnorm(n, sd = 0.25)
  x_cols <- paste0("X", seq_len(p))
  df <- as.data.frame(X)
  names(df) <- x_cols
  df$Y <- y
  list(X = X, X0 = X0, y = y, df = df, x_cols = x_cols)
}

# ============================================================================
# 1. Core mu equivalence: exact method, no SE
# ============================================================================

test_that("mu is identical across all input variants (exact, no SE)", {
  d   <- .make_np_iface_data()
  nms <- names(d$df)

  out_mat <- np_dnn(
    X = d$X, X0 = d$X0, Y = d$y,
    k = 2L, s = 3L, method = "exact", se = FALSE
  )
  # data-frame: column names for x_cols and y_col
  out_df_name <- np_dnn(
    data = d$df, x_cols = d$x_cols, y_col = "Y",
    X0 = d$X0, k = 2L, s = 3L, method = "exact", se = FALSE
  )
  # data-frame: integer indices for x_cols and y_col
  out_df_idx <- np_dnn(
    data    = d$df,
    x_cols  = match(d$x_cols, nms),
    y_col   = match("Y", nms),
    X0 = d$X0, k = 2L, s = 3L, method = "exact", se = FALSE
  )
  # data-frame: external y vector instead of y_col
  out_df_yvec <- np_dnn(
    data = d$df, x_cols = d$x_cols, y = d$y,
    X0 = d$X0, k = 2L, s = 3L, method = "exact", se = FALSE
  )

  expect_equal(out_df_name$mu,  out_mat$mu, tolerance = 1e-12)
  expect_equal(out_df_idx$mu,   out_mat$mu, tolerance = 1e-12)
  expect_equal(out_df_yvec$mu,  out_mat$mu, tolerance = 1e-12)
})

# ============================================================================
# 2. mu and SE equivalence: exact method, fast jackknife
# ============================================================================

test_that("mu and se are identical across interfaces (exact, fast jackknife)", {
  d <- .make_np_iface_data(n = 45L, m = 3L, seed = 7L)

  out_mat <- np_dnn(
    X = d$X, X0 = d$X0, Y = d$y,
    k = 2L, s = 3L, method = "exact", se = TRUE, jackknife = "fast"
  )
  out_df_name <- np_dnn(
    data = d$df, x_cols = d$x_cols, y_col = "Y",
    X0 = d$X0, k = 2L, s = 3L, method = "exact", se = TRUE, jackknife = "fast"
  )
  out_df_idx <- np_dnn(
    data   = d$df,
    x_cols = match(d$x_cols, names(d$df)),
    y_col  = match("Y", names(d$df)),
    X0 = d$X0, k = 2L, s = 3L, method = "exact", se = TRUE, jackknife = "fast"
  )
  out_df_yvec <- np_dnn(
    data = d$df, x_cols = d$x_cols, y = d$y,
    X0 = d$X0, k = 2L, s = 3L, method = "exact", se = TRUE, jackknife = "fast"
  )

  expect_equal(out_df_name$mu,  out_mat$mu, tolerance = 1e-12)
  expect_equal(out_df_name$se,  out_mat$se, tolerance = 1e-12)
  expect_equal(out_df_idx$mu,   out_mat$mu, tolerance = 1e-12)
  expect_equal(out_df_idx$se,   out_mat$se, tolerance = 1e-12)
  expect_equal(out_df_yvec$mu,  out_mat$mu, tolerance = 1e-12)
  expect_equal(out_df_yvec$se,  out_mat$se, tolerance = 1e-12)
})

# ============================================================================
# 3. mu and SE equivalence: exact method, loo jackknife
# ============================================================================

test_that("mu and se are identical across interfaces (exact, loo jackknife)", {
  d <- .make_np_iface_data(n = 30L, m = 2L, seed = 13L)

  out_mat <- np_dnn(
    X = d$X, X0 = d$X0, Y = d$y,
    k = 1L, s = 2L, method = "exact", se = TRUE, jackknife = "loo"
  )
  out_df <- np_dnn(
    data = d$df, x_cols = d$x_cols, y_col = "Y",
    X0 = d$X0, k = 1L, s = 2L, method = "exact", se = TRUE, jackknife = "loo"
  )

  expect_equal(out_df$mu, out_mat$mu, tolerance = 1e-12)
  expect_equal(out_df$se, out_mat$se, tolerance = 1e-12)
})

# ============================================================================
# 4. mu and SE equivalence: asymp method
# ============================================================================

test_that("mu and se are identical across interfaces (asymp, fast jackknife)", {
  d <- .make_np_iface_data(n = 48L, m = 3L, seed = 17L)

  out_mat <- np_dnn(
    X = d$X, X0 = d$X0, Y = d$y,
    k = 2L, s = 4L, method = "asymp", se = TRUE, jackknife = "fast"
  )
  out_df_name <- np_dnn(
    data = d$df, x_cols = d$x_cols, y_col = "Y",
    X0 = d$X0, k = 2L, s = 4L, method = "asymp", se = TRUE, jackknife = "fast"
  )
  out_df_idx <- np_dnn(
    data   = d$df,
    x_cols = match(d$x_cols, names(d$df)),
    y_col  = match("Y", names(d$df)),
    X0 = d$X0, k = 2L, s = 4L, method = "asymp", se = TRUE, jackknife = "fast"
  )

  expect_equal(out_df_name$mu, out_mat$mu, tolerance = 1e-12)
  expect_equal(out_df_name$se, out_mat$se, tolerance = 1e-12)
  expect_equal(out_df_idx$mu,  out_mat$mu, tolerance = 1e-12)
  expect_equal(out_df_idx$se,  out_mat$se, tolerance = 1e-12)
})

# ============================================================================
# 5. x0 path (single target): data-frame interface matches matrix interface
# ============================================================================

test_that("x0 path via data-frame interface matches matrix interface", {
  d <- .make_np_iface_data(n = 40L, p = 3L, m = 5L, seed = 55L)
  j <- 3L  # pick an interior target

  out_mat <- np_dnn(
    X = d$X, x0 = d$X0[j, ], Y = d$y,
    k = 2L, s = 3L, method = "exact", se = FALSE
  )
  out_df_name <- np_dnn(
    data = d$df, x_cols = d$x_cols, y_col = "Y",
    x0 = d$X0[j, ],
    k = 2L, s = 3L, method = "exact", se = FALSE
  )
  out_df_yvec <- np_dnn(
    data = d$df, x_cols = d$x_cols, y = d$y,
    x0 = d$X0[j, ],
    k = 2L, s = 3L, method = "exact", se = FALSE
  )

  expect_equal(out_df_name$mu, out_mat$mu, tolerance = 1e-12)
  expect_equal(out_df_yvec$mu, out_mat$mu, tolerance = 1e-12)
})

# ============================================================================
# 6. Return extras: W and ranks are identical across interfaces
# ============================================================================

test_that("return extras W and ranks are identical across interfaces", {
  d <- .make_np_iface_data(n = 42L, m = 3L, seed = 77L)

  out_mat <- np_dnn(
    X = d$X, X0 = d$X0, Y = d$y,
    k = 2L, s = 3L, method = "exact", se = FALSE,
    return = c("W", "ranks")
  )
  out_df <- np_dnn(
    data = d$df, x_cols = d$x_cols, y_col = "Y",
    X0 = d$X0, k = 2L, s = 3L, method = "exact", se = FALSE,
    return = c("W", "ranks")
  )

  expect_equal(out_df$W,     out_mat$W,     tolerance = 1e-12)
  expect_equal(out_df$ranks, out_mat$ranks)
})

# ============================================================================
# 7. Custom distance: equivalence is preserved under non-default distance
# ============================================================================

test_that("mu is identical across interfaces for non-default distance (l1)", {
  d <- .make_np_iface_data(n = 40L, m = 3L, seed = 33L)

  out_mat <- np_dnn(
    X = d$X, X0 = d$X0, Y = d$y,
    distance = "l1", k = 2L, s = 3L, method = "exact", se = FALSE
  )
  out_df <- np_dnn(
    data = d$df, x_cols = d$x_cols, y_col = "Y",
    X0 = d$X0, distance = "l1", k = 2L, s = 3L, method = "exact", se = FALSE
  )

  expect_equal(out_df$mu, out_mat$mu, tolerance = 1e-12)
})

# ============================================================================
# 8. Tibble acceptance: same result as data.frame
# ============================================================================

test_that("tibble input produces identical mu to data.frame input", {
  skip_if_not_installed("tibble")
  d <- .make_np_iface_data(seed = 88L)

  out_df  <- np_dnn(
    data = d$df, x_cols = d$x_cols, y_col = "Y",
    X0 = d$X0, k = 2L, s = 3L, method = "exact", se = FALSE
  )
  out_tbl <- np_dnn(
    data = tibble::as_tibble(d$df), x_cols = d$x_cols, y_col = "Y",
    X0 = d$X0, k = 2L, s = 3L, method = "exact", se = FALSE
  )

  expect_equal(out_tbl$mu, out_df$mu, tolerance = 1e-12)
})

# ============================================================================
# 9. Meta fields are consistent across interfaces
# ============================================================================

test_that("meta fields are identical across interfaces", {
  d <- .make_np_iface_data()

  out_mat <- np_dnn(
    X = d$X, X0 = d$X0, Y = d$y,
    k = 2L, s = 3L, method = "exact", se = FALSE
  )
  out_df <- np_dnn(
    data = d$df, x_cols = d$x_cols, y_col = "Y",
    X0 = d$X0, k = 2L, s = 3L, method = "exact", se = FALSE
  )

  expect_identical(out_df$meta$n,      out_mat$meta$n)
  expect_identical(out_df$meta$m,      out_mat$meta$m)
  expect_identical(out_df$meta$p,      out_mat$meta$p)
  expect_identical(out_df$meta$k,      out_mat$meta$k)
  expect_identical(out_df$meta$s,      out_mat$meta$s)
  expect_identical(out_df$meta$method, out_mat$meta$method)
})

# ============================================================================
# 10. x_cols index diagnostics: NA and non-finite are distinguished
# ============================================================================

test_that("NA integer index in x_cols produces NA-specific error", {
  d <- .make_np_iface_data()

  expect_error(
    np_dnn(
      data   = d$df,
      x_cols = c(1L, NA_integer_),
      y_col  = "Y",
      X0 = d$X0, k = 2L, s = 3L, se = FALSE
    ),
    regexp = "NA",
    ignore.case = TRUE
  )
})

test_that("non-finite index in x_cols produces non-finite-specific error", {
  d <- .make_np_iface_data()

  expect_error(
    np_dnn(
      data   = d$df,
      x_cols = c(1L, Inf),
      y_col  = "Y",
      X0 = d$X0, k = 2L, s = 3L, se = FALSE
    ),
    regexp = "non-finite",
    ignore.case = TRUE
  )
})

test_that("NA-index error is distinct from non-finite-index error", {
  d    <- .make_np_iface_data()
  base <- function(x_cols) {
    tryCatch(
      np_dnn(data = d$df, x_cols = x_cols, y_col = "Y",
             X0 = d$X0, k = 2L, s = 3L, se = FALSE),
      error = function(e) conditionMessage(e)
    )
  }

  msg_na  <- base(c(1L, NA_integer_))
  msg_inf <- base(c(1L, Inf))

  expect_false(
    identical(msg_na, msg_inf),
    info = "NA-index and non-finite-index errors should have distinct messages"
  )
  expect_match(msg_na,  "NA",         ignore.case = TRUE)
  expect_match(msg_inf, "non-finite", ignore.case = TRUE)
})

test_that("data-frame covariate values distinguish NA and non-finite errors", {
  d <- .make_np_iface_data()
  base <- function(df) {
    tryCatch(
      np_dnn(
        data = df, x_cols = d$x_cols, y_col = "Y",
        X0 = d$X0, k = 2L, s = 3L, se = FALSE
      ),
      error = conditionMessage
    )
  }

  df_na <- d$df
  df_na$X1[2L] <- NA_real_
  df_inf <- d$df
  df_inf$X1[2L] <- Inf

  msg_na <- base(df_na)
  msg_inf <- base(df_inf)

  expect_match(msg_na, "X.*NA", ignore.case = TRUE)
  expect_match(msg_inf, "X.*non-finite", ignore.case = TRUE)
  expect_false(identical(msg_na, msg_inf))
})

# ============================================================================
# 11. Missing column names produce informative errors
# ============================================================================

test_that("unknown column in x_cols names the missing column in the error", {
  d <- .make_np_iface_data()

  expect_error(
    np_dnn(
      data   = d$df,
      x_cols = c(d$x_cols[1L], "X_missing"),
      y_col  = "Y",
      X0 = d$X0, k = 2L, s = 3L, se = FALSE
    ),
    regexp = "X_missing",
    ignore.case = FALSE
  )
})

test_that("unknown y_col column name produces informative error", {
  d <- .make_np_iface_data()

  expect_error(
    np_dnn(
      data   = d$df,
      x_cols = d$x_cols,
      y_col  = "Y_typo",
      X0 = d$X0, k = 2L, s = 3L, se = FALSE
    ),
    regexp = "Y_typo",
    ignore.case = FALSE
  )
})

# ============================================================================
# 12. Out-of-range column indices produce informative errors
# ============================================================================

test_that("out-of-range integer index in x_cols produces informative error", {
  d <- .make_np_iface_data()
  p_too_large <- ncol(d$df) + 10L

  expect_error(
    np_dnn(
      data   = d$df,
      x_cols = c(1L, p_too_large),
      y_col  = "Y",
      X0 = d$X0, k = 2L, s = 3L, se = FALSE
    ),
    regexp = "out-of-range",
    ignore.case = TRUE
  )
})

# ============================================================================
# 13. Extended extras parity across interfaces
# ============================================================================

test_that("extended return extras are identical across interfaces", {
  d <- .make_np_iface_data(n = 44L, m = 3L, seed = 112L)

  out_mat <- np_dnn(
    X = d$X, X0 = d$X0, Y = d$y,
    k = 2L, s = 4L, method = "exact", se = FALSE,
    return = c("W", "ranks", "ord", "dx", "y", "ess", "neighbors")
  )
  out_df <- np_dnn(
    data = d$df, x_cols = d$x_cols, y_col = "Y",
    X0 = d$X0, k = 2L, s = 4L, method = "exact", se = FALSE,
    return = c("W", "ranks", "ord", "dx", "y", "ess", "neighbors")
  )

  expect_equal(out_df$mu, out_mat$mu, tolerance = 1e-12)
  expect_equal(out_df$W, out_mat$W, tolerance = 1e-12)
  expect_equal(out_df$ranks, out_mat$ranks)
  expect_equal(out_df$ord, out_mat$ord)
  expect_equal(out_df$dx, out_mat$dx, tolerance = 1e-12)
  expect_equal(out_df$y, out_mat$y, tolerance = 1e-12)
  expect_equal(out_df$ess, out_mat$ess, tolerance = 1e-12)
  expect_equal(out_df$neighbors, out_mat$neighbors)
})

# ============================================================================
# 14. normalize=FALSE parity
# ============================================================================

test_that("normalize=FALSE preserves matrix/data-frame equivalence", {
  d <- .make_np_iface_data(n = 46L, m = 3L, seed = 113L)

  out_mat <- np_dnn(
    X = d$X, X0 = d$X0, Y = d$y,
    k = 2L, s = 4L, method = "asymp", normalize = FALSE, se = FALSE,
    return = "W"
  )
  out_df <- np_dnn(
    data = d$df, x_cols = d$x_cols, y_col = "Y",
    X0 = d$X0, k = 2L, s = 4L, method = "asymp", normalize = FALSE, se = FALSE,
    return = "W"
  )

  expect_equal(out_df$mu, out_mat$mu, tolerance = 1e-12)
  expect_equal(out_df$W, out_mat$W, tolerance = 1e-12)
})

# ============================================================================
# 15. Distance spec object parity
# ============================================================================

test_that("distance spec object preserves matrix/data-frame equivalence", {
  d <- .make_np_iface_data(seed = 114L)
  dist_spec <- distance_l2()

  out_mat <- np_dnn(
    X = d$X, X0 = d$X0, Y = d$y,
    distance = dist_spec, k = 2L, s = 3L, method = "exact", se = FALSE
  )
  out_df <- np_dnn(
    data = d$df, x_cols = d$x_cols, y_col = "Y",
    X0 = d$X0, distance = dist_spec, k = 2L, s = 3L, method = "exact", se = FALSE
  )

  expect_equal(out_df$mu, out_mat$mu, tolerance = 1e-12)
})

# ============================================================================
# 16. Degenerate shape parity (p=1, m=1)
# ============================================================================

test_that("p=1,m=1 edge case remains equivalent across interfaces", {
  d <- .make_np_iface_data(n = 36L, p = 1L, m = 1L, seed = 115L)

  out_mat <- np_dnn(
    X = d$X, X0 = d$X0, Y = d$y,
    k = 1L, s = 2L, method = "exact", se = TRUE, jackknife = "fast"
  )
  out_df <- np_dnn(
    data = d$df, x_cols = d$x_cols, y_col = "Y",
    X0 = d$X0, k = 1L, s = 2L, method = "exact", se = TRUE, jackknife = "fast"
  )

  expect_equal(out_df$mu, out_mat$mu, tolerance = 1e-12)
  expect_equal(out_df$se, out_mat$se, tolerance = 1e-12)
})
