# tests/testthat/test-user-tdnn-np-data-interface.R
#
# Interface-equivalence tests for np_tdnn().
#
# Mirrors test-user-np-dnn-data-interface.R with TDNN parameters.
# These tests verify that the data-frame interface and the matrix
# interface produce numerically identical estimates and standard errors,
# and that the index-level diagnostics (NA, non-finite, out-of-range)
# are precise and distinguishable.

# ---- shared helper -----------------------------------------------------------

.make_tdnn_iface_data <- function(n = 55L, p = 2L, m = 4L, seed = 42L) {
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
# 1. Mean and metadata equivalence: exact method, no SE
# ============================================================================

test_that("means and metadata agree across input variants (exact, no SE)", {
  d   <- .make_tdnn_iface_data()
  nms <- names(d$df)

  out_mat <- np_tdnn(
    X = d$X, X0 = d$X0, Y = d$y,
    s1 = 2L, s2 = 4L, method = "exact", se = FALSE
  )
  # data-frame: column names for x_cols and y_col
  out_df_name <- np_tdnn(
    data = d$df, x_cols = d$x_cols, y_col = "Y",
    X0 = d$X0, s1 = 2L, s2 = 4L, method = "exact", se = FALSE
  )
  # data-frame: integer indices for x_cols and y_col
  out_df_idx <- np_tdnn(
    data   = d$df,
    x_cols = match(d$x_cols, nms),
    y_col  = match("Y", nms),
    X0 = d$X0, s1 = 2L, s2 = 4L, method = "exact", se = FALSE
  )
  # data-frame: external y vector instead of y_col
  out_df_yvec <- np_tdnn(
    data = d$df, x_cols = d$x_cols, y = d$y,
    X0 = d$X0, s1 = 2L, s2 = 4L, method = "exact", se = FALSE
  )

  expect_equal(out_df_name$mu,  out_mat$mu, tolerance = 1e-12)
  expect_equal(out_df_idx$mu,   out_mat$mu, tolerance = 1e-12)
  expect_equal(out_df_yvec$mu,  out_mat$mu, tolerance = 1e-12)

  expect_identical(out_df_name$meta$n,      out_mat$meta$n)
  expect_identical(out_df_name$meta$m,      out_mat$meta$m)
  expect_identical(out_df_name$meta$p,      out_mat$meta$p)
  expect_identical(out_df_name$meta$s1,     out_mat$meta$s1)
  expect_identical(out_df_name$meta$s2,     out_mat$meta$s2)
  expect_identical(out_df_name$meta$dim,    out_mat$meta$dim)
  expect_identical(out_df_name$meta$method, out_mat$meta$method)
})

# ============================================================================
# 2. mu and SE equivalence: exact method, fast jackknife
# ============================================================================

test_that("mu and se are identical across interfaces (exact, fast jackknife)", {
  d <- .make_tdnn_iface_data(n = 48L, m = 3L, seed = 7L)

  out_mat <- np_tdnn(
    X = d$X, X0 = d$X0, Y = d$y,
    s1 = 2L, s2 = 4L, method = "exact", se = TRUE, jackknife = "fast"
  )
  out_df_name <- np_tdnn(
    data = d$df, x_cols = d$x_cols, y_col = "Y",
    X0 = d$X0, s1 = 2L, s2 = 4L, method = "exact", se = TRUE, jackknife = "fast"
  )
  out_df_idx <- np_tdnn(
    data   = d$df,
    x_cols = match(d$x_cols, names(d$df)),
    y_col  = match("Y", names(d$df)),
    X0 = d$X0, s1 = 2L, s2 = 4L, method = "exact", se = TRUE, jackknife = "fast"
  )
  out_df_yvec <- np_tdnn(
    data = d$df, x_cols = d$x_cols, y = d$y,
    X0 = d$X0, s1 = 2L, s2 = 4L, method = "exact", se = TRUE, jackknife = "fast"
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
  d <- .make_tdnn_iface_data(n = 32L, m = 2L, seed = 13L)

  out_mat <- np_tdnn(
    X = d$X, X0 = d$X0, Y = d$y,
    s1 = 1L, s2 = 3L, method = "exact", se = TRUE, jackknife = "loo"
  )
  out_df <- np_tdnn(
    data = d$df, x_cols = d$x_cols, y_col = "Y",
    X0 = d$X0, s1 = 1L, s2 = 3L, method = "exact", se = TRUE, jackknife = "loo"
  )

  expect_equal(out_df$mu, out_mat$mu, tolerance = 1e-12)
  expect_equal(out_df$se, out_mat$se, tolerance = 1e-12)
})

# ============================================================================
# 4. mu and SE equivalence: asymp method
# ============================================================================

test_that("mu and se are identical across interfaces (asymp, fast jackknife)", {
  d <- .make_tdnn_iface_data(n = 50L, m = 3L, seed = 17L)

  out_mat <- np_tdnn(
    X = d$X, X0 = d$X0, Y = d$y,
    s1 = 2L, s2 = 5L, method = "asymp", se = TRUE, jackknife = "fast"
  )
  out_df_name <- np_tdnn(
    data = d$df, x_cols = d$x_cols, y_col = "Y",
    X0 = d$X0, s1 = 2L, s2 = 5L, method = "asymp", se = TRUE, jackknife = "fast"
  )
  out_df_idx <- np_tdnn(
    data   = d$df,
    x_cols = match(d$x_cols, names(d$df)),
    y_col  = match("Y", names(d$df)),
    X0 = d$X0, s1 = 2L, s2 = 5L, method = "asymp", se = TRUE, jackknife = "fast"
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
  d <- .make_tdnn_iface_data(n = 45L, p = 3L, m = 5L, seed = 55L)
  j <- 3L

  out_mat <- np_tdnn(
    X = d$X, x0 = d$X0[j, ], Y = d$y,
    s1 = 2L, s2 = 4L, method = "exact", se = FALSE
  )
  out_df_name <- np_tdnn(
    data = d$df, x_cols = d$x_cols, y_col = "Y",
    x0 = d$X0[j, ],
    s1 = 2L, s2 = 4L, method = "exact", se = FALSE
  )
  out_df_yvec <- np_tdnn(
    data = d$df, x_cols = d$x_cols, y = d$y,
    x0 = d$X0[j, ],
    s1 = 2L, s2 = 4L, method = "exact", se = FALSE
  )

  expect_equal(out_df_name$mu, out_mat$mu, tolerance = 1e-12)
  expect_equal(out_df_yvec$mu, out_mat$mu, tolerance = 1e-12)
})

# ============================================================================
# 6. Return extras W and ranks are identical across interfaces
# ============================================================================

test_that("return extras W and ranks are identical across interfaces", {
  d <- .make_tdnn_iface_data(n = 45L, m = 3L, seed = 77L)

  out_mat <- np_tdnn(
    X = d$X, X0 = d$X0, Y = d$y,
    s1 = 2L, s2 = 4L, method = "exact", se = FALSE,
    return = c("W", "ranks")
  )
  out_df <- np_tdnn(
    data = d$df, x_cols = d$x_cols, y_col = "Y",
    X0 = d$X0, s1 = 2L, s2 = 4L, method = "exact", se = FALSE,
    return = c("W", "ranks")
  )

  expect_equal(out_df$W,     out_mat$W,     tolerance = 1e-12)
  expect_equal(out_df$ranks, out_mat$ranks)
})

# ============================================================================
# 7. Custom distance: equivalence preserved under non-default distance
# ============================================================================

test_that("mu is identical across interfaces for non-default distance (l1)", {
  d <- .make_tdnn_iface_data(n = 45L, m = 3L, seed = 33L)

  out_mat <- np_tdnn(
    X = d$X, X0 = d$X0, Y = d$y,
    distance = "l1", s1 = 2L, s2 = 4L, method = "exact", se = FALSE
  )
  out_df <- np_tdnn(
    data = d$df, x_cols = d$x_cols, y_col = "Y",
    X0 = d$X0, distance = "l1", s1 = 2L, s2 = 4L, method = "exact", se = FALSE
  )

  expect_equal(out_df$mu, out_mat$mu, tolerance = 1e-12)
})

# ============================================================================
# 8. Tibble acceptance: same result as data.frame
# ============================================================================

test_that("tibble input produces identical mu to data.frame input", {
  skip_if_not_installed("tibble")
  d <- .make_tdnn_iface_data(seed = 88L)

  out_df  <- np_tdnn(
    data = d$df, x_cols = d$x_cols, y_col = "Y",
    X0 = d$X0, s1 = 2L, s2 = 4L, method = "exact", se = FALSE
  )
  out_tbl <- np_tdnn(
    data = tibble::as_tibble(d$df), x_cols = d$x_cols, y_col = "Y",
    X0 = d$X0, s1 = 2L, s2 = 4L, method = "exact", se = FALSE
  )

  expect_equal(out_tbl$mu, out_df$mu, tolerance = 1e-12)
})

# ============================================================================
# 9. x_cols index diagnostics: NA and non-finite are distinguished
# ============================================================================

test_that("NA-index error is distinct from non-finite-index error", {
  d <- .make_tdnn_iface_data()
  run <- function(x_cols) {
    np_tdnn(data = d$df, x_cols = x_cols, y_col = "Y",
            X0 = d$X0, s1 = 2L, s2 = 4L, se = FALSE)
  }

  err_na <- expect_error(run(c(1L, NA_integer_)), "NA", ignore.case = TRUE)
  err_inf <- expect_error(run(c(1L, Inf)), "non-finite", ignore.case = TRUE)
  expect_false(identical(conditionMessage(err_na), conditionMessage(err_inf)))
})

# ============================================================================
# 10. Missing column names produce informative errors
# ============================================================================

test_that("unknown column in x_cols names the missing column", {
  d <- .make_tdnn_iface_data()

  expect_error(
    np_tdnn(
      data   = d$df,
      x_cols = c(d$x_cols[1L], "X_missing"),
      y_col  = "Y",
      X0 = d$X0, s1 = 2L, s2 = 4L, se = FALSE
    ),
    regexp = "X_missing",
    ignore.case = FALSE
  )
})

test_that("unknown y_col column name produces informative error", {
  d <- .make_tdnn_iface_data()

  expect_error(
    np_tdnn(
      data   = d$df,
      x_cols = d$x_cols,
      y_col  = "Y_typo",
      X0 = d$X0, s1 = 2L, s2 = 4L, se = FALSE
    ),
    regexp = "Y_typo",
    ignore.case = FALSE
  )
})

# ============================================================================
# 11. Out-of-range column indices produce informative errors
# ============================================================================

test_that("out-of-range integer index in x_cols produces informative error", {
  d           <- .make_tdnn_iface_data()
  p_too_large <- ncol(d$df) + 10L

  expect_error(
    np_tdnn(
      data   = d$df,
      x_cols = c(1L, p_too_large),
      y_col  = "Y",
      X0 = d$X0, s1 = 2L, s2 = 4L, se = FALSE
    ),
    regexp = "out-of-range",
    ignore.case = TRUE
  )
})

# ============================================================================
# 12. Extended extras parity across interfaces
# ============================================================================

test_that("extended return extras are identical across interfaces", {
  d <- .make_tdnn_iface_data(n = 46L, m = 3L, seed = 116L)

  out_mat <- np_tdnn(
    X = d$X, X0 = d$X0, Y = d$y,
    s1 = 2L, s2 = 4L, method = "exact", se = FALSE,
    return = c("W", "ranks", "ord", "dx", "y", "ess", "neighbors")
  )
  out_df <- np_tdnn(
    data = d$df, x_cols = d$x_cols, y_col = "Y",
    X0 = d$X0, s1 = 2L, s2 = 4L, method = "exact", se = FALSE,
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
# 13. normalize=FALSE parity
# ============================================================================

test_that("normalize=FALSE preserves matrix/data-frame equivalence", {
  d <- .make_tdnn_iface_data(n = 50L, m = 3L, seed = 117L)

  out_mat <- np_tdnn(
    X = d$X, X0 = d$X0, Y = d$y,
    s1 = 2L, s2 = 5L, method = "asymp", normalize = FALSE, se = FALSE,
    return = "W"
  )
  out_df <- np_tdnn(
    data = d$df, x_cols = d$x_cols, y_col = "Y",
    X0 = d$X0, s1 = 2L, s2 = 5L, method = "asymp", normalize = FALSE, se = FALSE,
    return = "W"
  )

  expect_equal(out_df$mu, out_mat$mu, tolerance = 1e-12)
  expect_equal(out_df$W, out_mat$W, tolerance = 1e-12)
})

# ============================================================================
# 14. Distance spec object parity
# ============================================================================

test_that("distance spec object preserves matrix/data-frame equivalence", {
  d <- .make_tdnn_iface_data(seed = 118L)
  dist_spec <- distance_l2()

  out_mat <- np_tdnn(
    X = d$X, X0 = d$X0, Y = d$y,
    distance = dist_spec, s1 = 2L, s2 = 4L, method = "exact", se = FALSE
  )
  out_df <- np_tdnn(
    data = d$df, x_cols = d$x_cols, y_col = "Y",
    X0 = d$X0, distance = dist_spec, s1 = 2L, s2 = 4L, method = "exact", se = FALSE
  )

  expect_equal(out_df$mu, out_mat$mu, tolerance = 1e-12)
})

# ============================================================================
# 15. Degenerate shape parity (p=1, m=1)
# ============================================================================

test_that("p=1,m=1 edge case remains equivalent across interfaces", {
  d <- .make_tdnn_iface_data(n = 38L, p = 1L, m = 1L, seed = 119L)

  out_mat <- np_tdnn(
    X = d$X, X0 = d$X0, Y = d$y,
    s1 = 1L, s2 = 3L, method = "exact", se = TRUE, jackknife = "fast"
  )
  out_df <- np_tdnn(
    data = d$df, x_cols = d$x_cols, y_col = "Y",
    X0 = d$X0, s1 = 1L, s2 = 3L, method = "exact", se = TRUE, jackknife = "fast"
  )

  expect_equal(out_df$mu, out_mat$mu, tolerance = 1e-12)
  expect_equal(out_df$se, out_mat$se, tolerance = 1e-12)
})
