# tests/testthat/test-wrank-DNN-wrapper.R
# -------------------------------------------------------------------
# DNN rank-weight wrappers (full U, WOR)
# Tags:
#   [construct] S3 class & naming
#   [eval]      single-column eval shapes, sums, finiteness
#   [batch]     batch equals columnwise eval
#   [edge]      s=1,k=1 uniform; s=n top-k uniform (exact)
#   [errors]    invalid ranks; k > s -> error on eval
# -------------------------------------------------------------------

test_that("[construct] wrappers build rank_weight objects with informative names", {
  rw1 <- rank_weight_dnn_exact(k = 3, s = 8)
  rw2 <- rank_weight_dnn_asymp(k = 3, s = 8)
  expect_s3_class(rw1, "rank_weight")
  expect_s3_class(rw2, "rank_weight")
  expect_match(rw1$name, "^dnn_exact\\[k=3,s=8\\]$")
  expect_match(rw2$name, "^dnn_asymp\\[k=3,s=8\\]$")
})

test_that("[eval] eval() returns numeric(n), finite, sums to 1", {
  set.seed(1)
  n <- 15; s <- 7; k <- 3
  ranks <- sample.int(n)
  rw_ex <- rank_weight_dnn_exact(k, s)
  rw_ap <- rank_weight_dnn_asymp(k, s)

  w_ex <- rank_weight_eval(rw_ex, ranks)
  w_ap <- rank_weight_eval(rw_ap, ranks)

  expect_type(w_ex, "double"); expect_length(w_ex, n)
  expect_type(w_ap, "double"); expect_length(w_ap, n)
  expect_true(all(is.finite(w_ex)))
  expect_true(all(is.finite(w_ap)))
  expect_equal(sum(w_ex), 1, tolerance = 1e-12)
  expect_equal(sum(w_ap), 1, tolerance = 1e-12)
})

test_that("[batch] eval_batch equals columnwise eval", {
  set.seed(2)
  n <- 20; m <- 4; s <- 9; k <- 2
  R <- replicate(m, sample.int(n))
  rw_ex <- rank_weight_dnn_exact(k, s)
  rw_ap <- rank_weight_dnn_asymp(k, s)

  W_ex <- rank_weight_eval_batch(rw_ex, R)
  W_ap <- rank_weight_eval_batch(rw_ap, R)

  # columnwise recompute
  W_ex2 <- sapply(seq_len(m), function(j) rank_weight_eval(rw_ex, R[, j]))
  W_ap2 <- sapply(seq_len(m), function(j) rank_weight_eval(rw_ap, R[, j]))

  expect_equal(W_ex, W_ex2, tolerance = 0)
  expect_equal(W_ap, W_ap2, tolerance = 0)

  # column sums are 1
  expect_equal(colSums(W_ex), rep(1, m), tolerance = 1e-12)
  expect_equal(colSums(W_ap), rep(1, m), tolerance = 1e-12)
})

test_that("[edge] s = 1, k = 1 yields uniform weights (exact & asymp)", {
  set.seed(3)
  n <- 25; s <- 1; k <- 1
  ranks <- sample.int(n)

  rw_ex <- rank_weight_dnn_exact(k, s)
  rw_ap <- rank_weight_dnn_asymp(k, s)

  # data-order weights
  w_ex <- rank_weight_eval(rw_ex, ranks)
  w_ap <- rank_weight_eval(rw_ap, ranks)

  # both should be exactly/unambiguously uniform
  expect_equal(w_ex, rep(1/n, n), tolerance = 1e-12)
  expect_equal(w_ap, rep(1/n, n), tolerance = 1e-12)
})

test_that("[edge] s = n: exact puts 1/k on top-k ranks (zero else)", {
  set.seed(4)
  n <- 30; s <- n; k <- 5
  ranks <- sample.int(n)
  rw_ex <- rank_weight_dnn_exact(k, s)
  w <- rank_weight_eval(rw_ex, ranks)

  # Convert to rank order: w_rank[j] = weight for j-th nearest
  idx_rank <- order(ranks)   # data indices ordered by rank
  w_rank   <- w[idx_rank]

  expect_equal(w_rank[1:k], rep(1/k, k), tolerance = 0)
  expect_equal(w_rank[(k+1):n], rep(0, n-k), tolerance = 0)
})

test_that("[errors] invalid ranks or k > s raise errors on eval()", {
  n <- 8
  bad_ranks1 <- c(1L,2L,2L,4L,5L,6L,7L,8L)  # duplicate
  bad_ranks2 <- c(0L,2L,3L,4L,5L,6L,7L,8L)  # out of range
  rw_ex <- rank_weight_dnn_exact(k = 3, s = 2)  # k > s triggers inside exact fn

  expect_error(rank_weight_eval(rw_ex, bad_ranks1), "permutation")
  expect_error(rank_weight_eval(rw_ex, bad_ranks2), "permutation")

  ranks <- sample.int(n)
  # The underlying exact function checks k <= s and should error here:
  expect_error(rank_weight_eval(rw_ex, ranks), "k .*<=.* s|k .*must satisfy")
})

test_that("[batch][errors] eval_batch propagates input rank errors per column", {
  n <- 6; m <- 3; s <- 4; k <- 2
  R <- replicate(m, sample.int(n))
  R[, 2] <- c(1,1,3,4,5,6)  # bad column (duplicate)
  rw_ap <- rank_weight_dnn_asymp(k, s)

  expect_error(rank_weight_eval_batch(rw_ap, R), "permutation")
})
