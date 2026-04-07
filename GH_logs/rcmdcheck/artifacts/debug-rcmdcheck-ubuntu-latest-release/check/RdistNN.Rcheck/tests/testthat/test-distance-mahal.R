# tests/testthat/test-distance-mahal.R
# Tests for distance_mahalanobis() in dist_mahal.R.
#
# Coverage targets:
#   - correctness: identity and non-identity covariance
#   - setup() / eval() / eval_batch() contract
#   - eval matches column of eval_batch
#   - singular/near-singular covariance falls back to pseudoinverse
#   - kx < 2 rejected
#   - NA inputs rejected by eval and eval_batch
#   - "mahalanobis" string alias is registered at package load
#   - integration with cate_dnn() / weights_dnn()

skip_if_not(
  exists("distance_mahalanobis", mode = "function"),
  "distance_mahalanobis() not found; skipping file"
)

# ── correctness ───────────────────────────────────────────────────────────────

test_that("distance_mahalanobis: identity covariance matches L2", {
  # When cov(X) ≈ I, Mahalanobis = Euclidean
  set.seed(1)
  n  <- 120L
  X  <- matrix(c(rnorm(n), rnorm(n)), n, 2)   # independent columns  → cov ≈ I
  x0 <- c(0, 0)
  spec <- distance_mahalanobis()
  ctx  <- spec$setup(X)
  d_mahal <- spec$eval(X, x0, ctx)
  d_l2    <- sqrt(rowSums((X - matrix(x0, n, 2L, byrow = TRUE))^2))
  # Not exact (sample cov ≠ pop cov), but should be well within 15%
  expect_true(max(abs(d_mahal - d_l2) / pmax(d_l2, 1e-8)) < 0.15,
              info = "Mahalanobis should be close to L2 under identity covariance")
})

test_that("distance_mahalanobis: known solution with diagonal covariance", {
  # Construct X so that stats::cov(X) = diag(c(4, 1), 2) exactly.
  # stats::cov uses the n-1 (Bessel-corrected) denominator, so we need n=5:
  #   col1 = [2, 2, -2, -2, 0]: SS = 16, var = 16/4 = 4  ✓
  #   col2 = [1,-1,  1, -1, 0]: SS =  4, var =  4/4 = 1  ✓
  X <- matrix(c(2, 1, 2, -1, -2, 1, -2, -1, 0, 0), nrow = 5L, byrow = TRUE)
  # confirm:
  S <- stats::cov(X)
  expect_equal(S, diag(c(4, 1), 2), tolerance = 1e-12)

  spec <- distance_mahalanobis()
  ctx  <- spec$setup(X)
  x0   <- c(0, 0)
  d    <- spec$eval(X, x0, ctx)

  # expected: sqrt((xi1^2/4) + (xi2^2/1))
  expected <- sqrt(X[, 1]^2 / 4 + X[, 2]^2 / 1)
  expect_equal(d, expected, tolerance = 1e-10)
})

test_that("distance_mahalanobis: correlated design gives different distances than L2", {
  set.seed(42)
  n   <- 200L
  rho <- 0.9
  Z   <- matrix(rnorm(n * 2L), n, 2L)
  X   <- cbind(Z[, 1L], rho * Z[, 1L] + sqrt(1 - rho^2) * Z[, 2L])
  x0  <- c(0, 0)

  spec_mahal <- distance_mahalanobis()
  spec_l2    <- distance_l2()
  ctx_m      <- spec_mahal$setup(X)
  ctx_l2     <- spec_l2$setup(X)

  d_m <- spec_mahal$eval(X, x0, ctx_m)
  d_e <- spec_l2$eval(X, x0, ctx_l2)

  # With rho=0.9 the off-diagonal is large; Mahal and L2 should differ meaningfully
  expect_gt(mean(abs(d_m - d_e)), 0.01)
})

# ── setup / eval / eval_batch contract ───────────────────────────────────────

test_that("distance_mahalanobis: eval matches column 1 of eval_batch", {
  set.seed(2)
  X  <- matrix(rnorm(60L), 30L, 2L)
  X0 <- matrix(rnorm(10L), 5L, 2L)
  spec <- distance_mahalanobis()
  ctx  <- spec$setup(X)

  d_eval  <- spec$eval(X, X0[1L, ], ctx)
  D_batch <- spec$eval_batch(X, X0, ctx)

  expect_equal(as.numeric(d_eval), as.numeric(D_batch[, 1L]), tolerance = 1e-12)
})

test_that("distance_mahalanobis: eval_batch columns equal loop over eval", {
  set.seed(3)
  X  <- matrix(rnorm(80L), 40L, 2L)
  X0 <- matrix(rnorm(8L),  4L,  2L)
  spec <- distance_mahalanobis()
  ctx  <- spec$setup(X)

  D_batch <- spec$eval_batch(X, X0, ctx)
  D_loop  <- sapply(seq_len(nrow(X0)),
                    function(j) spec$eval(X, X0[j, ], ctx))

  expect_equal(D_batch, as.matrix(D_loop), tolerance = 1e-12)
})

test_that("distance_mahalanobis: setup called once is equivalent to NULL ctx", {
  set.seed(7)
  X  <- matrix(rnorm(60L), 30L, 2L)
  x0 <- c(0.5, -0.5)
  spec <- distance_mahalanobis()
  ctx  <- spec$setup(X)

  d_with_ctx  <- spec$eval(X, x0, ctx)
  d_null_ctx  <- spec$eval(X, x0, NULL)   # recomputes Si internally

  expect_equal(d_with_ctx, d_null_ctx, tolerance = 1e-10)
})

test_that("distance_mahalanobis: output is non-negative finite", {
  set.seed(9)
  X    <- matrix(rnorm(100L), 50L, 2L)
  x0   <- c(1, -1)
  spec <- distance_mahalanobis()
  ctx  <- spec$setup(X)
  d    <- spec$eval(X, x0, ctx)

  expect_true(all(is.finite(d)))
  expect_true(all(d >= 0))
})

test_that("distance_mahalanobis: self-distance is zero", {
  set.seed(11)
  X    <- matrix(rnorm(40L), 20L, 2L)
  spec <- distance_mahalanobis()
  ctx  <- spec$setup(X)

  for (i in c(1L, 5L, 20L)) {
    d <- spec$eval(X, X[i, ], ctx)
    expect_equal(d[i], 0, tolerance = 1e-10)
  }
})

# ── singular/near-singular covariance ─────────────────────────────────────────

test_that("distance_mahalanobis: near-singular covariance (constant column) falls back to pseudoinverse", {
  # X with one constant column → cov is singular
  set.seed(5)
  X <- matrix(rnorm(60L), 30L, 2L)
  X[, 2L] <- 3.0   # constant column → Sigma is rank-1, det = 0
  x0 <- c(0, 3)

  spec <- distance_mahalanobis()
  # Should not error; pseudoinverse handles the singularity
  expect_no_error({
    ctx <- spec$setup(X)
    d   <- spec$eval(X, x0, ctx)
  })
  expect_true(all(is.finite(d)))
  expect_true(all(d >= 0))
})

# ── input validation ──────────────────────────────────────────────────────────

test_that("distance_mahalanobis: rejects kx = 1", {
  spec <- distance_mahalanobis()
  X    <- matrix(rnorm(20L), 20L, 1L)
  expect_error(spec$setup(X), "at least 2 covariates")
})

test_that("distance_mahalanobis: rejects NA in X at eval", {
  set.seed(6)
  X    <- matrix(rnorm(40L), 20L, 2L)
  spec <- distance_mahalanobis()
  ctx  <- spec$setup(X)

  X_na <- X
  X_na[3L, 1L] <- NA_real_
  x0 <- c(0, 0)
  expect_error(spec$eval(X_na, x0, ctx))
})

test_that("distance_mahalanobis: rejects NA in X at eval_batch", {
  set.seed(6)
  X    <- matrix(rnorm(40L), 20L, 2L)
  spec <- distance_mahalanobis()
  ctx  <- spec$setup(X)

  X_na <- X
  X_na[1L, 2L] <- NA_real_
  X0 <- matrix(c(0, 0, 1, 1), 2L, 2L, byrow = TRUE)
  expect_error(spec$eval_batch(X_na, X0, ctx))
})

test_that("distance_mahalanobis: rejects NA in X0 at eval_batch", {
  set.seed(6)
  X    <- matrix(rnorm(40L), 20L, 2L)
  spec <- distance_mahalanobis()
  ctx  <- spec$setup(X)

  X0 <- matrix(c(NA_real_, 0), 1L, 2L)
  expect_error(spec$eval_batch(X, X0, ctx))
})

# ── string alias ──────────────────────────────────────────────────────────────

test_that("'mahalanobis' string alias is registered and evaluates correctly", {
  expect_true("mahalanobis" %in% distance_names())
  spec <- get_distance("mahalanobis")
  expect_true(is_distance(spec))
  expect_equal(spec$name, "mahalanobis")

  # Verify it produces finite, non-negative distances
  set.seed(8)
  X  <- matrix(rnorm(60L), 30L, 2L)
  x0 <- c(0, 0)
  ctx <- spec$setup(X)
  d   <- distance_eval(spec, X, x0, ctx)
  expect_true(all(is.finite(d)))
  expect_true(all(d >= 0))
})

# ── integration ───────────────────────────────────────────────────────────────

test_that("cate_dnn() accepts distance = 'mahalanobis'", {
  skip_if_not(exists("cate_dnn", mode = "function"), "cate_dnn() not found")
  set.seed(100)
  n   <- 120L
  rho <- 0.8
  Z   <- matrix(rnorm(n * 2L), n, 2L)
  X   <- cbind(Z[, 1L], rho * Z[, 1L] + sqrt(1 - rho^2) * Z[, 2L])
  psi <- X[, 1L] + rnorm(n, sd = 0.5)
  x0  <- c(0.5, 0.5)

  fit_mahal <- cate_dnn(
    X = X, x0 = x0, pseudo_outcome = psi,
    distance = "mahalanobis",
    k = 3L, s = 10L, method = "exact",
    se = TRUE, jackknife = "fast"
  )
  expect_true(is.numeric(fit_mahal$tau))
  expect_true(is.finite(fit_mahal$tau))
  expect_true(is.numeric(fit_mahal$se))
  expect_true(is.finite(fit_mahal$se))
  expect_gt(fit_mahal$se, 0)
})

test_that("cate_dnn() Mahalanobis vs L2: finite and structurally consistent", {
  skip_if_not(exists("cate_dnn", mode = "function"), "cate_dnn() not found")
  set.seed(101)
  n   <- 150L
  rho <- 0.85
  Z   <- matrix(rnorm(n * 2L), n, 2L)
  X   <- cbind(Z[, 1L], rho * Z[, 1L] + sqrt(1 - rho^2) * Z[, 2L])
  psi <- sin(2 * X[, 1L]) + rnorm(n, sd = 0.3)
  x0  <- c(0.5, 0)

  fit_l2    <- cate_dnn(X = X, x0 = x0, pseudo_outcome = psi,
                        distance = "l2",          k = 3L, s = 12L,
                        method = "exact", se = TRUE)
  fit_mahal <- cate_dnn(X = X, x0 = x0, pseudo_outcome = psi,
                        distance = "mahalanobis", k = 3L, s = 12L,
                        method = "exact", se = TRUE)

  # Both should return finite estimates; they will differ due to different neighborhoods
  expect_true(is.finite(fit_l2$tau))
  expect_true(is.finite(fit_mahal$tau))
  expect_true(is.finite(fit_l2$se))
  expect_true(is.finite(fit_mahal$se))
})

test_that("dnn_oracle_mahal benchmark method runs through sim_cate_one_rep", {
  skip_if_not(exists("sim_cate_one_rep", mode = "function"), "sim_cate_one_rep() not found")
  r <- sim_cate_one_rep(
    n = 120L, kx = 2L, k = 3L,
    s = floor(120L^(3/4)),
    design = "dgp5_correlated",
    methods = "dnn_oracle_mahal",
    seed = 42L
  )
  expect_s3_class(r, "data.frame")
  expect_equal(nrow(r), 1L)
  expect_equal(r$method, "dnn_oracle_mahal")
  expect_true(is.finite(r$tau_hat))
  expect_true(r$se > 0)
  expect_true(r$cover %in% c(0, 1))
})
