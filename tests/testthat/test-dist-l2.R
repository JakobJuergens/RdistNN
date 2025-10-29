# -------------------------------------------------------------------
# L2 distance implementation tests
# Tags:
#   [engine]  Shapes/finite & squared toggle
#   [compiled]Parity vs R fallback when compiled present
#   [batch]   distance_l2 batch parity
#   [order]   Deterministic ordering & squared order equivalence
# -------------------------------------------------------------------

test_that("[engine] dist_l2 returns correct shape; squared toggle works", {
  set.seed(1)
  X <- matrix(rnorm(30), 10, 3)
  x0 <- rnorm(3)
  d_sq <- dist_l2(X, x0, squared = TRUE)
  d <- dist_l2(X, x0, squared = FALSE)
  expect_length(d_sq, nrow(X))
  expect_true(all(is.finite(d_sq)))
  expect_true(all(is.finite(d)))
  expect_equal(d_sq, d^2, tolerance = 1e-12)
})

test_that("[compiled] compiled dist_l2 parity with R fallback", {
  set.seed(42)
  X <- matrix(rnorm(40), 20, 2)
  x0 <- rnorm(2)
  skip_if_not(exists("rustat_dist_l2", mode = "function"))
  # Force R fallback for reference
  ref_sq <- {
    dx2 <- rowSums((X - matrix(x0, nrow(X), ncol(X), TRUE))^2)
    dx2
  }
  ref <- sqrt(ref_sq)
  expect_equal(dist_l2(X, x0, squared = TRUE), ref_sq, tolerance = 1e-12)
  expect_equal(dist_l2(X, x0, squared = FALSE), ref, tolerance = 1e-12)
})

test_that("[batch] distance_l2 batch equals columnwise calls", {
  set.seed(7)
  X <- matrix(rnorm(24), 8, 3)
  X0 <- matrix(rnorm(6), 2, 3)
  dspec <- distance_l2(squared = TRUE)
  M <- distance_eval_batch(dspec, X, X0)
  M2 <- sapply(seq_len(nrow(X0)), function(j) distance_eval(dspec, X, X0[j, ]))
  expect_equal(M, M2, tolerance = 1e-12)
})

test_that("[order] order_by_l2 is deterministic and squared/unsquared yield same ord", {
  set.seed(99)
  X <- rbind(c(0, 0), c(1, 0), c(1, 0)) # tie on distances between rows 2 and 3
  x0 <- c(0, 0)
  o1 <- order_by_l2(X, x0, squared = FALSE)$ord
  o2 <- order_by_l2(X, x0, squared = TRUE)$ord
  expect_equal(o1, o2) # identical ordering
  expect_equal(o1, order(c(0, 1, 1), 1:3)) # tie broken by index
})

test_that("[registry] get_distance('l2') works and evaluates", {
  skip_if_not(exists("get_distance", mode = "function"))
  dspec <- get_distance("l2")
  expect_s3_class(dspec, "dnn_distance")
  X <- matrix(rnorm(12), 6, 2)
  x0 <- rnorm(2)
  out <- distance_eval(dspec, X, x0)
  expect_length(out, nrow(X))
  expect_true(all(is.finite(out)))
  expect_false(is_compiled_distance(dspec)) # ptr is NULL; compiled used inside eval()
})
