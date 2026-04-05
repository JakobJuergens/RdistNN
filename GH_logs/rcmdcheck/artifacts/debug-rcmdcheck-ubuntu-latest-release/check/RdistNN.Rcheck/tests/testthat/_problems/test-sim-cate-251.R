# Extracted from test-sim-cate.R:251

# setup ------------------------------------------------------------------------
library(testthat)
test_env <- simulate_test_env(package = "RdistNN", path = "..")
attach(test_env, warn.conflicts = FALSE)

# test -------------------------------------------------------------------------
dgp <- sim_cate_data("dgp5_correlated", n = 50L, kx = 2L, seed = 1L)
expect_equal(dgp$x0[1L], 0.75, tolerance = 1e-12)
