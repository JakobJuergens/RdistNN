# Estimator identities, retained extras and jackknife contracts.

test_that("np_tdnn returns valid shapes and weighted-mean identity", {
  .test_rng_scope({
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
  }, seed = 123L)
})

test_that("x0 path matches corresponding single-row X0 result", {
  .test_rng_scope({
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
  }, seed = 123L)
})

test_that("jackknife='none' omits se; constant y gives se = 0", {
  .test_rng_scope({
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
  }, seed = 123L)
})

test_that("fast and loo jackknife paths agree (exact)", {
  .test_rng_scope({
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
  }, seed = 123L)
})

test_that("fast and loo jackknife paths agree (asymp)", {
  .test_rng_scope({
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
  }, seed = 123L)
})

test_that("default jackknife mode matches jackknife='fast' explicitly", {
  .test_rng_scope({
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
  }, seed = 123L)
})

test_that("return extras are populated with correct shapes", {
  .test_rng_scope({
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
  }, seed = 123L)
})

test_that("TDNN weight columns sum to 1 and may contain negative entries", {
  .test_rng_scope({
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
  }, seed = 123L)
})

test_that("dim = NULL defaults to ncol(X) and gives same result as explicit dim", {
  .test_rng_scope({
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
  }, seed = 123L)
})
