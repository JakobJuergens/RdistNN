# Extracted from test-rank-weight-tdnn-asymp.R:130

# setup ------------------------------------------------------------------------
library(testthat)
test_env <- simulate_test_env(package = "RdistNN", path = "..")
attach(test_env, warn.conflicts = FALSE)

# prequel ----------------------------------------------------------------------
.make_rank_matrix <- function(n, m, seed = 1L) {
  set.seed(seed)
  R <- replicate(m, sample.int(n))
  storage.mode(R) <- "integer"
  R
}

# test -------------------------------------------------------------------------
n <- 15L
s1 <- 2L
s2 <- 6L
d <- 3L
ranks <- sample.int(n)
out <- dnn_weights_tdnn_asymp_full_from_ranks(ranks, s1, s2, d)
expect_identical(out$meta$n, n)
expect_identical(out$meta$s1, as.integer(s1))
expect_identical(out$meta$s2, as.integer(s2))
expect_identical(out$meta$dim, as.integer(d))
expect_true(out$meta$normalize)
expect_named(out$meta$scalars, c("w1", "w2"))
expect_equal(out$meta$scalars["w1"] + out$meta$scalars["w2"], 1, tolerance = 1e-14)
