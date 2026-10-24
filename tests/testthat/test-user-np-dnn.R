# Estimator identities, retained extras and jackknife contracts.

test_that("np_dnn basic matrix API returns valid shapes and weighted mean identity", {
  .test_rng_scope({
    d <- .make_np_data()
    out <- np_dnn(
      X = d$X, X0 = d$X0, Y = d$y,
      distance = "l2", method = "exact",
      k = 2L, s = 3L, se = TRUE,
      return = c("W", "ess")
    )

    expect_s3_class(out, "rdistnn_np_result")
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
  }, seed = 123L)
})

test_that("x0 path matches corresponding single-row X0 result", {
  .test_rng_scope({
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
  }, seed = 123L)
})

test_that("jackknife behavior: none omits se; constant outcomes give zero se", {
  .test_rng_scope({
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
  }, seed = 123L)
})

test_that("default jackknife mode uses fast path and matches loo reference", {
  .test_rng_scope({
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
  }, seed = 123L)
})

test_that("asymp method fast and loo jackknife paths agree", {
  .test_rng_scope({
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
  }, seed = 123L)
})

test_that("return extras are populated consistently", {
  .test_rng_scope({
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
  }, seed = 123L)
})

test_that("neighbors extra is nearest-k, not the full DNN weight support", {
  .test_rng_scope({
    d <- .make_np_data(n = 36L, p = 2L, m = 2L)
    k <- 1L
    s <- 8L
    out <- np_dnn(
      X = d$X, X0 = d$X0, Y = d$y,
      distance = "l2", method = "exact",
      k = k, s = s, se = FALSE,
      return = c("W", "ranks", "neighbors")
    )
    support_edge <- nrow(d$X) - s + k

    expect_equal(nrow(out$neighbors), k)
    for (j in seq_len(nrow(d$X0))) {
      expect_equal(out$neighbors[, j], order(out$ranks[, j])[seq_len(k)])
      expect_equal(sum(out$W[, j] > 0), support_edge)
      expect_gt(sum(out$W[, j] > 0), nrow(out$neighbors))
    }
  }, seed = 123L)
})

test_that("exact DNN collapses to equal nearest-k weights only at s = n", {
  .test_rng_scope({
    set.seed(125)
    d <- .make_np_data(n = 30L, p = 2L, m = 3L)
    k <- 3L

    boundary <- np_dnn(
      X = d$X, X0 = d$X0, Y = d$y,
      distance = "l2", method = "exact",
      k = k, s = nrow(d$X), se = FALSE,
      return = c("W", "ranks", "neighbors")
    )
    interior <- np_dnn(
      X = d$X, X0 = d$X0, Y = d$y,
      distance = "l2", method = "exact",
      k = k, s = 8L, se = FALSE,
      return = "W"
    )

    for (j in seq_len(nrow(d$X0))) {
      nearest <- order(boundary$ranks[, j])[seq_len(k)]
      positive <- which(boundary$W[, j] > 0)
      expect_equal(sort(positive), sort(nearest))
      expect_equal(boundary$W[nearest, j], rep(1 / k, k), tolerance = 1e-12)
      expect_equal(boundary$mu[[j]], mean(d$y[nearest]), tolerance = 1e-12)
      expect_equal(boundary$neighbors[, j], nearest)
      expect_gt(sum(interior$W[, j] > 0), k)
    }
  }, seed = 123L)
})
