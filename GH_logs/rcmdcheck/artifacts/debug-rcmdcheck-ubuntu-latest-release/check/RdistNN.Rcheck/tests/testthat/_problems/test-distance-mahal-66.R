# Extracted from test-distance-mahal.R:66

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
