# -------------------------------------------------------------------
# L1 distance implementation tests
# Tags:
#   [engine]  Shapes/finite
#   [compiled]Parity vs R fallback when compiled present
#   [batch]   distance_l1 batch parity
#   [order]   Deterministic ordering
#   [registry]Registry wiring
# -------------------------------------------------------------------

test_that("[engine] dist_l1 returns correct shape and finite values", {
  set.seed(1)
  X  <- matrix(rnorm(30), 10, 3)
  x0 <- rnorm(3)
  d  <- dist_l1(X, x0)
  expect_length(d, nrow(X))
  expect_true(all(is.finite(d)))
})

test_that("[compiled] compiled dist_l1 parity with R fallback", {
  set.seed(42)
  X  <- matrix(rnorm(40), 20, 2)
  x0 <- rnorm(2)
  skip_if_not(exists("rustat_dist_l1", mode = "function"))
  ref <- rowSums(abs(X - matrix(x0, nrow(X), ncol(X), TRUE)))
  expect_equal(dist_l1(X, x0), ref, tolerance = 0)
})

test_that("[batch] distance_l1 batch equals columnwise calls", {
  set.seed(7)
  X  <- matrix(rnorm(24), 8, 3)
  X0 <- matrix(rnorm(6), 2, 3)
  dspec <- distance_l1()
  M  <- distance_eval_batch(dspec, X, X0)
  M2 <- sapply(seq_len(nrow(X0)), function(j) distance_eval(dspec, X, X0[j, ]))
  expect_equal(M, M2, tolerance = 0)
})

test_that("[order] order_by_l1 is deterministic (index tie-break)", {
  X  <- rbind(c(0,0), c(1,0), c(1,0))
  x0 <- c(0,0)
  ob <- order_by_l1(X, x0)$ord
  expect_equal(ob, order(c(0,1,1), 1:3))
})

test_that("[registry] get_distance('l1') works and evaluates", {
  skip_if_not(exists("get_distance", mode = "function"))
  dspec <- get_distance("l1")
  expect_s3_class(dspec, "dnn_distance")
  X  <- matrix(rnorm(12), 6, 2)
  x0 <- rnorm(2)
  out <- distance_eval(dspec, X, x0)
  expect_length(out, nrow(X))
  expect_true(all(is.finite(out)))
})
