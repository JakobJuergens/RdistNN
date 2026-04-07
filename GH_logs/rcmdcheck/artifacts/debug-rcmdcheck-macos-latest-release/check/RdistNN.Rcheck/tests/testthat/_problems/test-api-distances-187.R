# Extracted from test-api-distances.R:187

# setup ------------------------------------------------------------------------
library(testthat)
test_env <- simulate_test_env(package = "RdistNN", path = "..")
attach(test_env, warn.conflicts = FALSE)

# prequel ----------------------------------------------------------------------
skip_if_not(
  all(vapply(
    c(
      "is_distance", "as_distance", "register_distance", "get_distance",
      "distance_names", "distance_eval", "distance_eval_batch"
    ),
    exists, logical(1),
    mode = "function"
  )),
  "distance API not loaded"
)

# test -------------------------------------------------------------------------
bad <- as_distance(function(X, x0, ctx) 1, name = "bad_len")
X <- matrix(0, 3, 2)
x0 <- c(0, 0)
expect_error(distance_eval(bad, X, x0), "must return numeric length nrow\\(X\\)")
