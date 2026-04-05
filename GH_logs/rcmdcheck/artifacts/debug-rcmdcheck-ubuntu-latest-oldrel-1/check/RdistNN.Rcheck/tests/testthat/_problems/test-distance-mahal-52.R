# Extracted from test-distance-mahal.R:52

# setup ------------------------------------------------------------------------
library(testthat)
test_env <- simulate_test_env(package = "RdistNN", path = "..")
attach(test_env, warn.conflicts = FALSE)

# prequel ----------------------------------------------------------------------
skip_if_not(
  exists("distance_mahalanobis", mode = "function"),
  "distance_mahalanobis() not found; skipping file"
)

# test -------------------------------------------------------------------------
X <- matrix(c(2, 1, 2, -1, -2, 1, -2, -1), nrow = 4, byrow = TRUE)
S <- stats::cov(X)
expect_equal(S, diag(c(4, 1), 2), tolerance = 1e-12)
spec <- distance_mahalanobis()
ctx  <- spec$setup(X)
x0   <- c(0, 0)
d    <- spec$eval(X, x0, ctx)
expected <- sqrt(X[, 1]^2 / 4 + X[, 2]^2 / 1)
expect_equal(d, expected, tolerance = 1e-10)
