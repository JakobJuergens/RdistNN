# Extracted from test-tdnn-jackknife.R:258

# setup ------------------------------------------------------------------------
library(testthat)
test_env <- simulate_test_env(package = "RdistNN", path = "..")
attach(test_env, warn.conflicts = FALSE)

# prequel ----------------------------------------------------------------------
.make_X_jk <- function(n, p = 2L, seed = 1L) {
  set.seed(seed)
  matrix(rnorm(n * p), n, p)
}
.make_x0_jk <- function(p = 2L) rep(0, p)
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

# test -------------------------------------------------------------------------
n <- 10L
p <- 2L
m <- 1L
X  <- .make_X_jk(n, p, seed = 62L)
X0 <- .make_X0(m, p, seed = 63L)
