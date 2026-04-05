# Extracted from test-rank-weight-tdnn-exact.R:59

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
for (params in list(c(2, 4, 2), c(3, 9, 3), c(5, 10, 1))) {
    s1 <- params[1]; s2 <- params[2]; d <- params[3]
    sc <- RdistNN:::.tdnn_scalars(s1, s2, d)
    expect_equal(sc["w1"] + sc["w2"], 1, tolerance = 1e-14)
    expect_lt(sc["w1"], 0)
    expect_gt(sc["w2"], 1)
  }
