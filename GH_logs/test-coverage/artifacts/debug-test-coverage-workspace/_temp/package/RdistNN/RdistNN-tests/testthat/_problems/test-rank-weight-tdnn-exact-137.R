# Extracted from test-rank-weight-tdnn-exact.R:137

# setup ------------------------------------------------------------------------
library(testthat)
test_env <- simulate_test_env(package = "RdistNN", path = "..")
attach(test_env, warn.conflicts = FALSE)

# prequel ----------------------------------------------------------------------
.tdnn_exact_brute <- function(n, s1, s2, dim) {
  # Compute the TDNN exact rank table directly from the definition:
  #   w1* * w_DNN(s1) + w2* * w_DNN(s2)
  # using .dnn_rank_table_exact_k1 via dnn_weights_exact_full_from_ranks.
  sc <- RdistNN:::.tdnn_scalars(s1, s2, dim)
  t1 <- dnn_weights_exact_full_from_ranks(seq_len(n), k = 1L, s = s1)$w_rank
  t2 <- dnn_weights_exact_full_from_ranks(seq_len(n), k = 1L, s = s2)$w_rank
  sc["w1"] * t1 + sc["w2"] * t2
}

# test -------------------------------------------------------------------------
n <- 14L
s1 <- 2L
s2 <- 5L
d <- 3L
ranks <- sample.int(n)
out <- dnn_weights_tdnn_exact_full_from_ranks(ranks, s1, s2, d)
expect_identical(out$meta$n, n)
expect_identical(out$meta$s1, as.integer(s1))
expect_identical(out$meta$s2, as.integer(s2))
expect_identical(out$meta$dim, as.integer(d))
expect_named(out$meta$scalars, c("w1", "w2"))
expect_equal(out$meta$scalars["w1"] + out$meta$scalars["w2"], 1, tolerance = 1e-14)
