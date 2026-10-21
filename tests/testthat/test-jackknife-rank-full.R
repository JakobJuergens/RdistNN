# Generic full rank jackknife: structure, explicit deletions and validation.

test_that("[jkfull] jackknife_rank_full returns correct structure", {
  .test_rng_scope({
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
  }, seed = 44L)
})

test_that("[jkfull] jackknife_rank_full with q=1 still returns array dim c(n,1,m)", {
  .test_rng_scope({
    set.seed(41L)
    n <- 12L; m <- 2L
    ranks <- apply(matrix(runif(n * m), n, m), 2L, rank)
    storage.mode(ranks) <- "integer"
    PHI  <- matrix(rnorm(n), n, 1L)
    wn   <- rep(1 / n, n)
    wnm1 <- rep(1 / (n - 1L), n - 1L)
    out  <- jackknife_rank_full(ranks, PHI, wn, wnm1)
    expect_equal(dim(out$theta_loo), c(n, 1L, m))
  }, seed = 44L)
})

test_that("[jkfull] constant signal: theta0 = constant, all loo = constant", {
  .test_rng_scope({
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
  }, seed = 44L)
})

test_that("[jkfull] theta0 equals naive sum (uniform weights)", {
  .test_rng_scope({
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
  }, seed = 44L)
})

test_that("[jkfull] brute-force LOO parity (small n, q=1)", {
  .test_rng_scope({
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
  }, seed = 44L)
})

test_that("[jkfull] input validation: wn wrong length", {
  .test_rng_scope({
    n <- 10L; m <- 2L
    ranks <- apply(matrix(runif(n * m), n, m), 2L, rank)
    storage.mode(ranks) <- "integer"
    PHI  <- matrix(rnorm(n), n, 1L)
    wn   <- rep(1 / n, n + 1L)   # wrong length
    wnm1 <- rep(1 / (n - 1L), n - 1L)
    expect_error(jackknife_rank_full(ranks, PHI, wn, wnm1), "wn")
  }, seed = 44L)
})

test_that("[jkfull] input validation: wnm1 wrong length", {
  .test_rng_scope({
    n <- 10L; m <- 2L
    ranks <- apply(matrix(runif(n * m), n, m), 2L, rank)
    storage.mode(ranks) <- "integer"
    PHI  <- matrix(rnorm(n), n, 1L)
    wn   <- rep(1 / n, n)
    wnm1 <- rep(1 / (n - 1L), n)  # wrong length (should be n-1)
    expect_error(jackknife_rank_full(ranks, PHI, wn, wnm1), "wnm1")
  }, seed = 44L)
})

test_that("[jkfull] input validation: PHI wrong nrow", {
  .test_rng_scope({
    n <- 10L; m <- 2L
    ranks <- apply(matrix(runif(n * m), n, m), 2L, rank)
    storage.mode(ranks) <- "integer"
    PHI  <- matrix(rnorm(n + 1L), n + 1L, 1L)  # wrong nrow
    wn   <- rep(1 / n, n)
    wnm1 <- rep(1 / (n - 1L), n - 1L)
    expect_error(jackknife_rank_full(ranks, PHI, wn, wnm1), "PHI")
  }, seed = 44L)
})

test_that("[jkfull] input validation rejects malformed ranks and non-finite inputs", {
  .test_rng_scope({
    n <- 10L; m <- 2L
    ranks <- apply(matrix(runif(n * m), n, m), 2L, rank)
    storage.mode(ranks) <- "integer"
    PHI  <- matrix(rnorm(n), n, 1L)
    wn   <- rep(1 / n, n)
    wnm1 <- rep(1 / (n - 1L), n - 1L)

    ranks_dup <- ranks
    ranks_dup[1L, 1L] <- ranks_dup[2L, 1L]
    expect_error(jackknife_rank_full(ranks_dup, PHI, wn, wnm1), "permutation")

    ranks_frac <- ranks
    storage.mode(ranks_frac) <- "double"
    ranks_frac[1L, 1L] <- 1.5
    expect_error(jackknife_rank_full(ranks_frac, PHI, wn, wnm1), "integer ranks")

    PHI_bad <- PHI
    PHI_bad[1L, 1L] <- Inf
    expect_error(jackknife_rank_full(ranks, PHI_bad, wn, wnm1), "PHI.*non-finite")

    wn_bad <- wn
    wn_bad[1L] <- NA_real_
    expect_error(jackknife_rank_full(ranks, PHI, wn_bad, wnm1), "wn.*NA")

    wnm1_bad <- wnm1
    wnm1_bad[1L] <- Inf
    expect_error(jackknife_rank_full(ranks, PHI, wn, wnm1_bad), "wnm1.*non-finite")
  }, seed = 44L)
})

test_that("[jkfull] R fallback produces same result as dispatch", {
  .test_rng_scope({
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
  }, seed = 44L)
})
