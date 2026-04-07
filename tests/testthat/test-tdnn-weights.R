# tests/testthat/test-tdnn-weights.R
# Tests for weights_tdnn_exact(), weights_tdnn_asymp(), weights_tdnn() in R/api-weights.R
# and jackknife_rank_full() in R/jackknife-rank-full.R.
#
# Coverage targets:
#   [exact]     weights_tdnn_exact: shape, colSums, signed, extras, x0 vs X0, dim default
#   [asymp]     weights_tdnn_asymp: same battery; normalize flag
#   [dispatch]  weights_tdnn: exact/asymp dispatch agrees with direct calls
#   [constant]  constant signal → estimate equals constant (all methods)
#   [jkfull]    jackknife_rank_full: shape, constant-signal, brute-force parity, validation

# ── helpers ───────────────────────────────────────────────────────────────────

.make_X <- function(n, p, seed = 1L) {
  set.seed(seed)
  matrix(rnorm(n * p), n, p)
}

.make_X0 <- function(m, p, seed = 99L) {
  set.seed(seed)
  matrix(rnorm(m * p), m, p)
}

# ── [exact] ───────────────────────────────────────────────────────────────────

test_that("[exact] weights_tdnn_exact returns W with correct shape", {
  n <- 20L; p <- 2L; m <- 3L
  X  <- .make_X(n, p, seed = 1L)
  X0 <- .make_X0(m, p, seed = 2L)
  out <- weights_tdnn_exact(X, X0 = X0, s1 = 3L, s2 = 8L)
  expect_true(is.list(out))
  expect_true(!is.null(out$W))
  expect_equal(dim(out$W), c(n, m))
})

test_that("[exact] weights_tdnn_exact colSums = 1", {
  n <- 25L; p <- 2L
  X  <- .make_X(n, p, seed = 3L)
  x0 <- c(0, 0)
  out <- weights_tdnn_exact(X, x0 = x0, s1 = 4L, s2 = 10L)
  expect_equal(colSums(out$W), 1, tolerance = 1e-10)
})

test_that("[exact] weights_tdnn_exact has negative weights (bias-correction)", {
  n <- 30L; p <- 2L
  X  <- .make_X(n, p, seed = 4L)
  x0 <- c(0.5, -0.5)
  out <- weights_tdnn_exact(X, x0 = x0, s1 = 3L, s2 = 9L)
  expect_true(any(out$W < 0))
})

test_that("[exact] weights_tdnn_exact x0 and X0 interfaces agree", {
  n <- 20L; p <- 2L
  X  <- .make_X(n, p, seed = 5L)
  x0 <- c(0.1, -0.2)
  X0 <- matrix(x0, nrow = 1L)
  out_x0 <- weights_tdnn_exact(X, x0 = x0, s1 = 3L, s2 = 7L)
  out_X0 <- weights_tdnn_exact(X, X0 = X0, s1 = 3L, s2 = 7L)
  expect_equal(out_x0$W, out_X0$W, tolerance = 0)
})

test_that("[exact] weights_tdnn_exact dim defaults to ncol(X)", {
  n <- 18L; p <- 3L
  X  <- .make_X(n, p, seed = 6L)
  x0 <- rep(0, p)
  out_default  <- weights_tdnn_exact(X, x0 = x0, s1 = 2L, s2 = 6L)
  out_explicit <- weights_tdnn_exact(X, x0 = x0, s1 = 2L, s2 = 6L, dim = p)
  expect_equal(out_default$W, out_explicit$W, tolerance = 0)
})

test_that("[exact] weights_tdnn_exact with NULL X and NULL dim gives a clear error", {
  # Before the fix, NULL X with dim = NULL produced integer(0) for dim (via
  # ncol(NULL) = NULL -> as.integer(NULL) = integer(0)), which cascaded into a
  # confusing sprintf warning or deferred failure deep in eval_batch.
  # After the fix, an explicit NULL guard fires immediately.
  expect_error(
    weights_tdnn_exact(NULL, x0 = c(0, 0), s1 = 2L, s2 = 6L),
    "`dim`"
  )
})

test_that("[exact] weights_tdnn_exact return='ranks' includes rank matrix", {
  n <- 15L; p <- 2L; m <- 2L
  X  <- .make_X(n, p, seed = 7L)
  X0 <- .make_X0(m, p, seed = 8L)
  out <- weights_tdnn_exact(X, X0 = X0, s1 = 2L, s2 = 5L, return = c("W", "ranks"))
  expect_equal(dim(out$ranks), c(n, m))
  expect_true(is.integer(out$ranks))
})

test_that("[exact] weights_tdnn_exact return='ord' and 'dx' included when requested", {
  n <- 15L; p <- 2L
  X  <- .make_X(n, p, seed = 9L)
  x0 <- c(0, 0)
  out <- weights_tdnn_exact(X, x0 = x0, s1 = 2L, s2 = 6L, return = c("W", "ord", "dx"))
  expect_false(is.null(out$ord))
  expect_false(is.null(out$dx))
  expect_equal(dim(out$ord), c(n, 1L))
  expect_equal(dim(out$dx),  c(n, 1L))
})

# ── [asymp] ───────────────────────────────────────────────────────────────────

test_that("[asymp] weights_tdnn_asymp returns W with correct shape", {
  n <- 20L; p <- 2L; m <- 3L
  X  <- .make_X(n, p, seed = 10L)
  X0 <- .make_X0(m, p, seed = 11L)
  out <- weights_tdnn_asymp(X, X0 = X0, s1 = 3L, s2 = 8L)
  expect_equal(dim(out$W), c(n, m))
})

test_that("[asymp] weights_tdnn_asymp colSums = 1 (normalize=TRUE)", {
  n <- 25L; p <- 2L
  X  <- .make_X(n, p, seed = 12L)
  x0 <- c(0, 0)
  out <- weights_tdnn_asymp(X, x0 = x0, s1 = 4L, s2 = 10L)
  expect_equal(colSums(out$W), 1, tolerance = 1e-10)
})

test_that("[asymp] weights_tdnn_asymp has negative weights", {
  n <- 30L; p <- 2L
  X  <- .make_X(n, p, seed = 13L)
  x0 <- c(0.5, -0.5)
  out <- weights_tdnn_asymp(X, x0 = x0, s1 = 3L, s2 = 9L)
  expect_true(any(out$W < 0))
})

test_that("[asymp] weights_tdnn_asymp normalize=FALSE sign structure preserved", {
  n <- 20L; p <- 2L
  X  <- .make_X(n, p, seed = 14L)
  x0 <- c(0, 0)
  out_norm <- weights_tdnn_asymp(X, x0 = x0, s1 = 3L, s2 = 7L, normalize = TRUE)
  out_raw  <- weights_tdnn_asymp(X, x0 = x0, s1 = 3L, s2 = 7L, normalize = FALSE)
  expect_equal(sign(out_norm$W), sign(out_raw$W))
  expect_equal(colSums(out_norm$W), 1, tolerance = 1e-10)
})

test_that("[asymp] weights_tdnn_asymp x0 and X0 interfaces agree", {
  n <- 20L; p <- 2L
  X  <- .make_X(n, p, seed = 15L)
  x0 <- c(0.1, -0.2)
  X0 <- matrix(x0, nrow = 1L)
  out_x0 <- weights_tdnn_asymp(X, x0 = x0, s1 = 3L, s2 = 7L)
  out_X0 <- weights_tdnn_asymp(X, X0 = X0, s1 = 3L, s2 = 7L)
  expect_equal(out_x0$W, out_X0$W, tolerance = 0)
})

test_that("[asymp] weights_tdnn_asymp dim defaults to ncol(X)", {
  n <- 18L; p <- 3L
  X  <- .make_X(n, p, seed = 16L)
  x0 <- rep(0, p)
  out_default  <- weights_tdnn_asymp(X, x0 = x0, s1 = 2L, s2 = 6L)
  out_explicit <- weights_tdnn_asymp(X, x0 = x0, s1 = 2L, s2 = 6L, dim = p)
  expect_equal(out_default$W, out_explicit$W, tolerance = 0)
})

test_that("[asymp] weights_tdnn_asymp with NULL X and NULL dim gives a clear error", {
  expect_error(
    weights_tdnn_asymp(NULL, x0 = c(0, 0), s1 = 2L, s2 = 6L),
    "`dim`"
  )
})

# ── [dispatch] ────────────────────────────────────────────────────────────────

test_that("[dispatch] weights_tdnn(method='exact') equals weights_tdnn_exact", {
  n <- 20L; p <- 2L; m <- 2L
  X  <- .make_X(n, p, seed = 20L)
  X0 <- .make_X0(m, p, seed = 21L)
  out_direct <- weights_tdnn_exact(X, X0 = X0, s1 = 3L, s2 = 8L)
  out_unified <- weights_tdnn(X, X0 = X0, s1 = 3L, s2 = 8L, method = "exact")
  expect_equal(out_unified$W, out_direct$W, tolerance = 0)
})

test_that("[dispatch] weights_tdnn(method='asymp') equals weights_tdnn_asymp", {
  n <- 20L; p <- 2L; m <- 2L
  X  <- .make_X(n, p, seed = 22L)
  X0 <- .make_X0(m, p, seed = 23L)
  out_direct  <- weights_tdnn_asymp(X, X0 = X0, s1 = 3L, s2 = 8L)
  out_unified <- weights_tdnn(X, X0 = X0, s1 = 3L, s2 = 8L, method = "asymp")
  expect_equal(out_unified$W, out_direct$W, tolerance = 0)
})

test_that("[dispatch] weights_tdnn default method is 'exact'", {
  n <- 20L; p <- 2L
  X  <- .make_X(n, p, seed = 24L)
  x0 <- c(0, 0)
  out_default <- weights_tdnn(X, x0 = x0, s1 = 3L, s2 = 8L)
  out_exact   <- weights_tdnn(X, x0 = x0, s1 = 3L, s2 = 8L, method = "exact")
  expect_equal(out_default$W, out_exact$W, tolerance = 0)
})

# ── [constant] ────────────────────────────────────────────────────────────────

test_that("[constant] constant signal: exact weighted mean = constant", {
  set.seed(30L)
  n <- 20L; p <- 2L; m <- 3L; c0 <- 7.5
  X     <- .make_X(n, p, seed = 30L)
  X0    <- .make_X0(m, p, seed = 31L)
  out   <- weights_tdnn_exact(X, X0 = X0, s1 = 3L, s2 = 8L)
  psi   <- rep(c0, n)
  ests  <- as.numeric(crossprod(psi, out$W))   # W^T psi (length m)
  expect_equal(ests, rep(c0, m), tolerance = 1e-10)
})

test_that("[constant] constant signal: asymp weighted mean = constant", {
  n <- 20L; p <- 2L; m <- 3L; c0 <- -2.1
  X     <- .make_X(n, p, seed = 32L)
  X0    <- .make_X0(m, p, seed = 33L)
  out   <- weights_tdnn_asymp(X, X0 = X0, s1 = 3L, s2 = 8L)
  psi   <- rep(c0, n)
  ests  <- as.numeric(crossprod(psi, out$W))
  expect_equal(ests, rep(c0, m), tolerance = 1e-10)
})

# ── [jkfull] ──────────────────────────────────────────────────────────────────

test_that("[jkfull] jackknife_rank_full returns correct structure", {
  set.seed(40L)
  n <- 15L; m <- 3L; q <- 2L
  ranks <- apply(matrix(runif(n * m), n, m), 2L, rank)
  storage.mode(ranks) <- "integer"
  PHI  <- matrix(rnorm(n * q), n, q)
  wn   <- rep(1 / n, n)
  wnm1 <- rep(1 / (n - 1L), n - 1L)
  out  <- jackknife_rank_full(ranks, PHI, wn, wnm1)

  expect_type(out, "list")
  expect_named(out, c("theta0", "theta_loo", "info"))
  expect_equal(dim(out$theta0),    c(q, m))
  expect_equal(dim(out$theta_loo), c(n, q, m))
  expect_equal(out$info$n, n)
  expect_equal(out$info$m, m)
  expect_equal(out$info$q, q)
})

test_that("[jkfull] jackknife_rank_full with q=1 still returns array dim c(n,1,m)", {
  set.seed(41L)
  n <- 12L; m <- 2L
  ranks <- apply(matrix(runif(n * m), n, m), 2L, rank)
  storage.mode(ranks) <- "integer"
  PHI  <- matrix(rnorm(n), n, 1L)
  wn   <- rep(1 / n, n)
  wnm1 <- rep(1 / (n - 1L), n - 1L)
  out  <- jackknife_rank_full(ranks, PHI, wn, wnm1)
  expect_equal(dim(out$theta_loo), c(n, 1L, m))
})

test_that("[jkfull] constant signal: theta0 = constant, all loo = constant", {
  set.seed(42L)
  n <- 15L; m <- 3L; c0 <- 5.0
  ranks <- apply(matrix(runif(n * m), n, m), 2L, rank)
  storage.mode(ranks) <- "integer"
  PHI  <- matrix(c0, n, 1L)
  wn   <- rep(1 / n, n)
  wnm1 <- rep(1 / (n - 1L), n - 1L)
  out  <- jackknife_rank_full(ranks, PHI, wn, wnm1)
  expect_equal(as.numeric(out$theta0), rep(c0, m), tolerance = 1e-10)
  expect_equal(as.numeric(out$theta_loo), rep(c0, n * m), tolerance = 1e-10)
})

test_that("[jkfull] theta0 equals naive sum (uniform weights)", {
  # With uniform wn, theta0[t] = mean(PHI), independent of ranks
  set.seed(43L)
  n <- 15L; m <- 3L; q <- 2L
  ranks <- apply(matrix(runif(n * m), n, m), 2L, rank)
  storage.mode(ranks) <- "integer"
  PHI  <- matrix(rnorm(n * q), n, q)
  wn   <- rep(1 / n, n)
  wnm1 <- rep(1 / (n - 1L), n - 1L)
  out  <- jackknife_rank_full(ranks, PHI, wn, wnm1)
  # theta0[qq, t] = sum_i wn[r_it] * PHI[i, qq] = (1/n) sum_i PHI[i,qq] = mean(PHI[,qq])
  expected <- matrix(colMeans(PHI), q, m)
  expect_equal(out$theta0, expected, tolerance = 1e-12)
})

test_that("[jkfull] brute-force LOO parity (small n, q=1)", {
  # Verify against a direct loop for a non-trivial weight vector
  set.seed(44L)
  n <- 8L; m <- 2L
  ranks <- apply(matrix(runif(n * m), n, m), 2L, rank)
  storage.mode(ranks) <- "integer"
  PHI  <- matrix(rnorm(n), n, 1L)
  # Use TDNN-asymp rank table as wn/wnm1
  wn   <- tdnn_rank_table_asymp(n,       2L, 5L, 2L)
  wnm1 <- tdnn_rank_table_asymp(n - 1L,  2L, 5L, 2L)
  out  <- jackknife_rank_full(ranks, PHI, wn, wnm1)

  # Brute-force: delete obs i, rebuild weighted estimate
  loo_ref <- matrix(NA_real_, n, m)
  for (t in seq_len(m)) {
    for (i in seq_len(n)) {
      keep <- seq_len(n)[-i]
      r_keep <- ranks[keep, t]
      # Ranks in reduced dataset: decrement ranks above i's rank
      r_i    <- ranks[i, t]
      r_red  <- r_keep - (r_keep > r_i)
      loo_ref[i, t] <- sum(wnm1[r_red] * PHI[keep, 1L])
    }
  }
  expect_equal(out$theta_loo[, 1L, ], loo_ref, tolerance = 1e-12)
})

test_that("[jkfull] input validation: wn wrong length", {
  n <- 10L; m <- 2L
  ranks <- apply(matrix(runif(n * m), n, m), 2L, rank)
  storage.mode(ranks) <- "integer"
  PHI  <- matrix(rnorm(n), n, 1L)
  wn   <- rep(1 / n, n + 1L)   # wrong length
  wnm1 <- rep(1 / (n - 1L), n - 1L)
  expect_error(jackknife_rank_full(ranks, PHI, wn, wnm1), "wn")
})

test_that("[jkfull] input validation: wnm1 wrong length", {
  n <- 10L; m <- 2L
  ranks <- apply(matrix(runif(n * m), n, m), 2L, rank)
  storage.mode(ranks) <- "integer"
  PHI  <- matrix(rnorm(n), n, 1L)
  wn   <- rep(1 / n, n)
  wnm1 <- rep(1 / (n - 1L), n)  # wrong length (should be n-1)
  expect_error(jackknife_rank_full(ranks, PHI, wn, wnm1), "wnm1")
})

test_that("[jkfull] input validation: PHI wrong nrow", {
  n <- 10L; m <- 2L
  ranks <- apply(matrix(runif(n * m), n, m), 2L, rank)
  storage.mode(ranks) <- "integer"
  PHI  <- matrix(rnorm(n + 1L), n + 1L, 1L)  # wrong nrow
  wn   <- rep(1 / n, n)
  wnm1 <- rep(1 / (n - 1L), n - 1L)
  expect_error(jackknife_rank_full(ranks, PHI, wn, wnm1), "PHI")
})

test_that("[jkfull] R fallback produces same result as dispatch", {
  set.seed(45L)
  n <- 12L; m <- 2L; q <- 2L
  ranks <- apply(matrix(runif(n * m), n, m), 2L, rank)
  storage.mode(ranks) <- "integer"
  PHI  <- matrix(rnorm(n * q), n, q)
  wn   <- tdnn_rank_table_asymp(n,       2L, 5L, 2L)
  wnm1 <- tdnn_rank_table_asymp(n - 1L,  2L, 5L, 2L)
  out_dispatch <- jackknife_rank_full(ranks, PHI, wn, wnm1)
  out_r        <- RdistNN:::.jackknife_rank_full_R(ranks, PHI, wn, wnm1)
  expect_equal(out_dispatch$theta0,    out_r$theta0,    tolerance = 1e-14)
  expect_equal(out_dispatch$theta_loo, out_r$theta_loo, tolerance = 1e-14)
})
