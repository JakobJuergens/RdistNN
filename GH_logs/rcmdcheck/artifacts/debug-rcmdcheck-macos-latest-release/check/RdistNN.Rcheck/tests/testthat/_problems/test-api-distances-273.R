# Extracted from test-api-distances.R:273

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
X <- matrix(1:12, nrow = 4)
X0 <- matrix(c(0, 0, 0, 1, 1, 1), nrow = 2, byrow = TRUE)
bad_shape <- list(
    name = "bad_shape",
    eval_batch = function(X, X0, ctx) rep(0, nrow(X) * nrow(X0))
  )
class(bad_shape) <- "dnn_distance"
expect_error(
    distance_eval_batch(bad_shape, X, X0),
    "must return a finite numeric matrix with shape nrow\\(X\\) x nrow\\(X0\\)"
  )
