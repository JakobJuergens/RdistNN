# tests/testthat/test-tdnn-jackknife.R
# Tests for .tdnn_jackknife_se_loo_matrix() and .tdnn_jackknife_se_loo_full()
# in R/user-tdnn-core.R.
#
# Coverage targets:
#   [se-shape]    SE matrix is q x m with finite, non-negative entries
#   [constant]    constant signal -> SE = 0 (all engines)
#   [agreement]   fast path agrees with loo fallback
#   [full-path]   _loo_full agrees with _loo_matrix (fast) for both exact/asymp
#   [q>1]         vector signal (q=2) works in all engines
#   [precomp]     precomputed_ranks produces same result as internal ranking
#   [s2-guard]    s2 >= n rejected by both engines
#   [method]      exact and asymp paths both work end-to-end

# ── helpers ───────────────────────────────────────────────────────────────────

.make_X_jk <- function(n, p = 2L, seed = 1L) {
  set.seed(seed)
  matrix(rnorm(n * p), n, p)
}

.make_x0_jk <- function(p = 2L) rep(0, p)

# thin wrapper: call the internal function with scalar x0 wrapped as x0_mat
.se_matrix <- function(X, x0_mat, signal, s1, s2, d,
                        method = "exact", jackknife = "fast",
                        precomputed_ranks = NULL) {
  RdistNN:::.tdnn_jackknife_se_loo_matrix(
    x_mat             = X,
    x0_mat            = x0_mat,
    signal            = signal,
    distance          = "l2",
    method            = method,
    s1                = s1,
    s2                = s2,
    dim               = d,
    normalize         = TRUE,
    ctx               = NULL,
    signal_arg        = "signal",
    jackknife         = jackknife,
    precomputed_ranks = precomputed_ranks
  )
}

.se_full <- function(X, x0_mat, signal, s1, s2, d,
                      method = "exact",
                      precomputed_ranks = NULL) {
  RdistNN:::.tdnn_jackknife_se_loo_full(
    x_mat             = X,
    x0_mat            = x0_mat,
    signal            = signal,
    distance          = "l2",
    method            = method,
    s1                = s1,
    s2                = s2,
    dim               = d,
    normalize         = TRUE,
    ctx               = NULL,
    signal_arg        = "signal",
    precomputed_ranks = precomputed_ranks
  )
}

# ── [se-shape] ────────────────────────────────────────────────────────────────

test_that("[se-shape] _loo_matrix fast: SE matrix is 1 x m, finite, non-negative", {
  n <- 25L; p <- 2L; m <- 3L; s1 <- 3L; s2 <- 8L; d <- 2L
  X    <- .make_X_jk(n, p, seed = 1L)
  set.seed(2L); X0 <- matrix(rnorm(m * p), m, p)
  set.seed(3L); psi <- rnorm(n)
  se <- .se_matrix(X, X0, psi, s1, s2, d, method = "exact", jackknife = "fast")
  expect_equal(dim(se), c(1L, m))
  expect_true(all(is.finite(se)))
  expect_true(all(se >= 0))
})

test_that("[se-shape] _loo_full: SE matrix is 1 x m, finite, non-negative", {
  n <- 25L; p <- 2L; m <- 3L; s1 <- 3L; s2 <- 8L; d <- 2L
  X    <- .make_X_jk(n, p, seed = 4L)
  set.seed(5L); X0 <- matrix(rnorm(m * p), m, p)
  set.seed(6L); psi <- rnorm(n)
  se <- .se_full(X, X0, psi, s1, s2, d, method = "exact")
  expect_equal(dim(se), c(1L, m))
  expect_true(all(is.finite(se)))
  expect_true(all(se >= 0))
})

test_that("[se-shape] _loo_matrix loo fallback: SE is 1 x m, finite, non-negative", {
  n <- 20L; p <- 2L; m <- 2L; s1 <- 2L; s2 <- 6L; d <- 2L
  X  <- .make_X_jk(n, p, seed = 7L)
  set.seed(8L); X0 <- matrix(rnorm(m * p), m, p)
  set.seed(9L); psi <- rnorm(n)
  se <- .se_matrix(X, X0, psi, s1, s2, d, method = "exact", jackknife = "loo")
  expect_equal(dim(se), c(1L, m))
  expect_true(all(is.finite(se)))
  expect_true(all(se >= 0))
})

# ── [constant] ────────────────────────────────────────────────────────────────

test_that("[constant] constant signal: fast path SE = 0", {
  n <- 20L; p <- 2L; m <- 2L; s1 <- 3L; s2 <- 8L; d <- 2L; c0 <- 4.5
  X  <- .make_X_jk(n, p, seed = 10L)
  set.seed(11L); X0 <- matrix(rnorm(m * p), m, p)
  psi <- rep(c0, n)
  se  <- .se_matrix(X, X0, psi, s1, s2, d, method = "exact", jackknife = "fast")
  expect_equal(as.numeric(se), rep(0, m), tolerance = 1e-10)
})

test_that("[constant] constant signal: loo fallback SE = 0", {
  n <- 20L; p <- 2L; m <- 2L; s1 <- 3L; s2 <- 8L; d <- 2L; c0 <- -1.5
  X  <- .make_X_jk(n, p, seed = 12L)
  set.seed(13L); X0 <- matrix(rnorm(m * p), m, p)
  psi <- rep(c0, n)
  se  <- .se_matrix(X, X0, psi, s1, s2, d, method = "exact", jackknife = "loo")
  expect_equal(as.numeric(se), rep(0, m), tolerance = 1e-10)
})

test_that("[constant] constant signal: _loo_full SE = 0", {
  n <- 20L; p <- 2L; m <- 2L; s1 <- 3L; s2 <- 8L; d <- 2L; c0 <- 2.0
  X  <- .make_X_jk(n, p, seed = 14L)
  set.seed(15L); X0 <- matrix(rnorm(m * p), m, p)
  psi <- rep(c0, n)
  se  <- .se_full(X, X0, psi, s1, s2, d, method = "exact")
  expect_equal(as.numeric(se), rep(0, m), tolerance = 1e-10)
})

# ── [agreement] ───────────────────────────────────────────────────────────────

test_that("[agreement] fast path and loo fallback agree (exact, q=1)", {
  n <- 20L; p <- 2L; m <- 2L; s1 <- 3L; s2 <- 7L; d <- 2L
  X  <- .make_X_jk(n, p, seed = 20L)
  set.seed(21L); X0 <- matrix(rnorm(m * p), m, p)
  set.seed(22L); psi <- rnorm(n)
  se_fast <- .se_matrix(X, X0, psi, s1, s2, d, method = "exact", jackknife = "fast")
  se_loo  <- .se_matrix(X, X0, psi, s1, s2, d, method = "exact", jackknife = "loo")
  expect_equal(se_fast, se_loo, tolerance = 1e-10)
})

test_that("[agreement] fast path and loo fallback agree (asymp, q=1)", {
  n <- 20L; p <- 2L; m <- 2L; s1 <- 3L; s2 <- 7L; d <- 2L
  X  <- .make_X_jk(n, p, seed = 23L)
  set.seed(24L); X0 <- matrix(rnorm(m * p), m, p)
  set.seed(25L); psi <- rnorm(n)
  se_fast <- .se_matrix(X, X0, psi, s1, s2, d, method = "asymp", jackknife = "fast")
  se_loo  <- .se_matrix(X, X0, psi, s1, s2, d, method = "asymp", jackknife = "loo")
  expect_equal(se_fast, se_loo, tolerance = 1e-10)
})

# ── [full-path] ───────────────────────────────────────────────────────────────

test_that("[full-path] _loo_full agrees with fast matrix path (exact)", {
  n <- 25L; p <- 2L; m <- 3L; s1 <- 3L; s2 <- 9L; d <- 2L
  X  <- .make_X_jk(n, p, seed = 30L)
  set.seed(31L); X0 <- matrix(rnorm(m * p), m, p)
  set.seed(32L); psi <- rnorm(n)
  se_matrix <- .se_matrix(X, X0, psi, s1, s2, d, method = "exact", jackknife = "fast")
  se_full   <- .se_full(X, X0, psi, s1, s2, d, method = "exact")
  expect_equal(se_full, se_matrix, tolerance = 1e-10)
})

test_that("[full-path] _loo_full agrees with fast matrix path (asymp)", {
  n <- 25L; p <- 2L; m <- 3L; s1 <- 3L; s2 <- 9L; d <- 2L
  X  <- .make_X_jk(n, p, seed = 33L)
  set.seed(34L); X0 <- matrix(rnorm(m * p), m, p)
  set.seed(35L); psi <- rnorm(n)
  se_matrix <- .se_matrix(X, X0, psi, s1, s2, d, method = "asymp", jackknife = "fast")
  se_full   <- .se_full(X, X0, psi, s1, s2, d, method = "asymp")
  expect_equal(se_full, se_matrix, tolerance = 1e-10)
})

# ── [q>1] ─────────────────────────────────────────────────────────────────────

test_that("[q>1] vector signal (q=2) fast path: SE is 2 x m, finite, non-negative", {
  n <- 25L; p <- 2L; m <- 3L; q <- 2L; s1 <- 3L; s2 <- 8L; d <- 2L
  X  <- .make_X_jk(n, p, seed = 40L)
  set.seed(41L); X0 <- matrix(rnorm(m * p), m, p)
  set.seed(42L); sig <- matrix(rnorm(n * q), n, q)
  se <- .se_matrix(X, X0, sig, s1, s2, d, method = "exact", jackknife = "fast")
  expect_equal(dim(se), c(q, m))
  expect_true(all(is.finite(se)))
  expect_true(all(se >= 0))
})

test_that("[q>1] vector signal (q=2): constant columns have SE=0", {
  n <- 20L; p <- 2L; m <- 2L; s1 <- 3L; s2 <- 7L; d <- 2L
  X  <- .make_X_jk(n, p, seed = 43L)
  set.seed(44L); X0 <- matrix(rnorm(m * p), m, p)
  sig <- cbind(rep(1.5, n), rnorm(n))   # first column constant, second not
  se  <- .se_matrix(X, X0, sig, s1, s2, d, method = "exact", jackknife = "fast")
  expect_equal(dim(se), c(2L, m))
  expect_equal(as.numeric(se[1L, ]), rep(0, m), tolerance = 1e-10)
  # second row SE should be positive (non-constant signal)
  expect_true(all(se[2L, ] > 0))
})

test_that("[q>1] _loo_full with q=2 agrees with fast path", {
  n <- 20L; p <- 2L; m <- 2L; q <- 2L; s1 <- 3L; s2 <- 7L; d <- 2L
  X  <- .make_X_jk(n, p, seed = 45L)
  set.seed(46L); X0 <- matrix(rnorm(m * p), m, p)
  set.seed(47L); sig <- matrix(rnorm(n * q), n, q)
  se_fast <- .se_matrix(X, X0, sig, s1, s2, d, method = "exact", jackknife = "fast")
  se_full <- .se_full(X, X0, sig, s1, s2, d, method = "exact")
  expect_equal(se_full, se_fast, tolerance = 1e-10)
})

# ── [precomp] ─────────────────────────────────────────────────────────────────

test_that("[precomp] precomputed_ranks gives same SE as internal ranking (matrix fast)", {
  n <- 20L; p <- 2L; m <- 3L; s1 <- 3L; s2 <- 7L; d <- 2L
  X  <- .make_X_jk(n, p, seed = 50L)
  set.seed(51L); X0 <- matrix(rnorm(m * p), m, p)
  set.seed(52L); psi <- rnorm(n)
  # Precompute ranks via weights_tdnn
  preranks <- weights_tdnn(X, X0 = X0, s1 = s1, s2 = s2, dim = d,
                           method = "exact", return = "ranks")$ranks
  se_internal <- .se_matrix(X, X0, psi, s1, s2, d, method = "exact",
                             jackknife = "fast", precomputed_ranks = NULL)
  se_precomp  <- .se_matrix(X, X0, psi, s1, s2, d, method = "exact",
                             jackknife = "fast", precomputed_ranks = preranks)
  expect_equal(se_precomp, se_internal, tolerance = 1e-12)
})

test_that("[precomp] precomputed_ranks gives same SE as internal ranking (_loo_full)", {
  n <- 20L; p <- 2L; m <- 2L; s1 <- 3L; s2 <- 7L; d <- 2L
  X  <- .make_X_jk(n, p, seed = 53L)
  set.seed(54L); X0 <- matrix(rnorm(m * p), m, p)
  set.seed(55L); psi <- rnorm(n)
  preranks <- weights_tdnn(X, X0 = X0, s1 = s1, s2 = s2, dim = d,
                           method = "exact", return = "ranks")$ranks
  se_internal <- .se_full(X, X0, psi, s1, s2, d, method = "exact",
                          precomputed_ranks = NULL)
  se_precomp  <- .se_full(X, X0, psi, s1, s2, d, method = "exact",
                          precomputed_ranks = preranks)
  expect_equal(se_precomp, se_internal, tolerance = 1e-12)
})

# ── [s2-guard] ────────────────────────────────────────────────────────────────

test_that("[s2-guard] s2 >= n is rejected by _loo_matrix", {
  n <- 10L; p <- 2L; m <- 1L
  X  <- .make_X_jk(n, p, seed = 60L)
  X0 <- .make_X0(m, p, seed = 61L)
  psi <- rnorm(n)
  expect_error(
    .se_matrix(X, X0, psi, s1 = 2L, s2 = n, d = 2L),
    "s2 < n"
  )
  expect_error(
    .se_matrix(X, X0, psi, s1 = 2L, s2 = n + 1L, d = 2L),
    "s2 < n"
  )
})

test_that("[s2-guard] s2 >= n is rejected by _loo_full", {
  n <- 10L; p <- 2L; m <- 1L
  X  <- .make_X_jk(n, p, seed = 62L)
  X0 <- .make_X0(m, p, seed = 63L)
  psi <- rnorm(n)
  expect_error(
    .se_full(X, X0, psi, s1 = 2L, s2 = n, d = 2L),
    "s2 < n"
  )
})

# ── [method] ──────────────────────────────────────────────────────────────────

test_that("[method] exact and asymp engines agree for large n", {
  # For large n the exact and asymp rank tables are close; SE should be close too
  n <- 200L; p <- 2L; m <- 2L; s1 <- 5L; s2 <- 15L; d <- 2L
  X  <- .make_X_jk(n, p, seed = 70L)
  set.seed(71L); X0 <- matrix(rnorm(m * p), m, p)
  set.seed(72L); psi <- rnorm(n)
  se_exact <- .se_matrix(X, X0, psi, s1, s2, d, method = "exact",  jackknife = "fast")
  se_asymp <- .se_matrix(X, X0, psi, s1, s2, d, method = "asymp",  jackknife = "fast")
  # Not identical, but should be close
  expect_lt(max(abs(se_exact - se_asymp)), 0.05 * max(se_exact))
})

test_that("[method] asymp _loo_full returns finite non-negative SE", {
  n <- 25L; p <- 2L; m <- 2L; s1 <- 3L; s2 <- 8L; d <- 2L
  X  <- .make_X_jk(n, p, seed = 73L)
  set.seed(74L); X0 <- matrix(rnorm(m * p), m, p)
  set.seed(75L); psi <- rnorm(n)
  se <- .se_full(X, X0, psi, s1, s2, d, method = "asymp")
  expect_true(all(is.finite(se)))
  expect_true(all(se >= 0))
})
