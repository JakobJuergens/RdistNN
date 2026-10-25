# tests/testthat/test-rank-weight-tdnn-asymp.R
# Tests for R/wrank-tdnn_approx.R
#
# Coverage targets:
#   [table]       tdnn_rank_table_asymp: length, sum, signed, full support
#   [vs-exact]    asymp converges to exact as n grows
#   [spec]        rank_weight_tdnn_asymp: class, caching, eval_batch, colSums
#   [parity]      asymp spec agrees with from_ranks
#   [errors]      input validation for table and from_ranks

# ── helpers ───────────────────────────────────────────────────────────────────

.make_rank_matrix <- function(n, m, seed = 1L) {
  .test_rng_scope({
    R <- replicate(m, sample.int(n))
    storage.mode(R) <- "integer"
    R
  }, seed = seed)
}

# ── [table] ───────────────────────────────────────────────────────────────────

test_that("[table] tdnn_rank_table_asymp returns length-n numeric summing to 1", {
  n <- 25L; s1 <- 3L; s2 <- 8L; d <- 2L
  w <- tdnn_rank_table_asymp(n, s1, s2, d)
  expect_type(w, "double")
  expect_length(w, n)
  expect_equal(sum(w), 1, tolerance = 1e-10)
  expect_true(is.finite(sum(w)))
})

test_that("[table] tdnn_rank_table_asymp has signed weights (bias-correction)", {
  n <- 30L; s1 <- 4L; s2 <- 10L; d <- 2L
  w <- tdnn_rank_table_asymp(n, s1, s2, d)
  expect_true(any(w < 0))
  expect_equal(sum(w), 1, tolerance = 1e-10)
})

test_that("[table] tdnn_rank_table_asymp has full support (all entries non-zero)", {
  # Unlike the exact variant, asymp assigns non-zero weight to all n ranks
  n <- 20L; s1 <- 3L; s2 <- 7L; d <- 2L
  w <- tdnn_rank_table_asymp(n, s1, s2, d)
  expect_true(all(w != 0))
})

test_that("[table] normalize=FALSE preserves sign structure but relaxes sum", {
  n <- 20L; s1 <- 3L; s2 <- 7L; d <- 2L
  w_norm  <- tdnn_rank_table_asymp(n, s1, s2, d, normalize = TRUE)
  w_raw   <- tdnn_rank_table_asymp(n, s1, s2, d, normalize = FALSE)
  # Signs agree
  expect_equal(sign(w_norm), sign(w_raw))
  # Normalized version sums to 1; raw should be close but not guaranteed
  expect_equal(sum(w_norm), 1, tolerance = 1e-10)
})

test_that("[table] tdnn_rank_table_asymp equals w1*·DNN_asymp(s1) + w2*·DNN_asymp(s2)", {
  n <- 20L; s1 <- 3L; s2 <- 8L; d <- 2L
  sc <- RdistNN:::.tdnn_scalars(s1, s2, d)
  t1 <- RdistNN:::.dnn_rank_table_asymp_k1(n, s1, normalize = FALSE)
  t2 <- RdistNN:::.dnn_rank_table_asymp_k1(n, s2, normalize = FALSE)
  w_ref  <- sc["w1"] * t1 + sc["w2"] * t2
  w_ref  <- w_ref / sum(w_ref)
  w_tdnn <- tdnn_rank_table_asymp(n, s1, s2, d)
  expect_equal(w_tdnn, w_ref, tolerance = 1e-12)
})

test_that("[table] tdnn_rank_table_asymp validates inputs", {
  expect_error(tdnn_rank_table_asymp(0L,  2L, 4L, 2L), "positive integer")
  expect_error(tdnn_rank_table_asymp(10L, 4L, 4L, 2L), "s1 < s2")
  expect_error(tdnn_rank_table_asymp(10L, 5L, 4L, 2L), "s1 < s2")
  expect_error(tdnn_rank_table_asymp(10L, 2L, 11L, 2L), "s2 <= n")
  expect_error(tdnn_rank_table_asymp(10L, 2L, 4L, 0L), "positive integer")
  expect_error(tdnn_rank_table_asymp(10L, 2L, 4L, 2L, normalize = NA), "normalize.*TRUE or FALSE")
  expect_error(tdnn_rank_table_asymp(10L, 2L, 4L, 2L, normalize = "yes"), "normalize.*TRUE or FALSE")
  expect_error(
    tdnn_rank_table_asymp(10L, 2L, 4L, 2L, normalize = c(TRUE, FALSE)),
    "normalize.*TRUE or FALSE"
  )
})

# ── [vs-exact] ────────────────────────────────────────────────────────────────

test_that("[vs-exact] asymp converges to exact as n grows", {
  # For n >> s2 the Beta approximation is accurate; max entry difference shrinks
  s1 <- 3L; s2 <- 8L; d <- 2L
  diffs <- vapply(c(50L, 200L, 1000L), function(n) {
    we <- tdnn_rank_table_exact(n, s1, s2, d)
    wa <- tdnn_rank_table_asymp(n, s1, s2, d)
    max(abs(we - wa))
  }, numeric(1))
  # Differences should shrink monotonically (or at least the large-n diff < small-n diff)
  expect_lt(diffs[3L], diffs[1L])
})

# ── [from_ranks] ──────────────────────────────────────────────────────────────








# ── [spec-errors] rank_weight_tdnn_asymp constructor ─────────────────────────

test_that("[spec-errors] rank_weight_tdnn_asymp validates s1, s2, dim at construction time", {
  # non-positive arguments
  expect_error(rank_weight_tdnn_asymp(s1 = 0L, s2 = 4L, dim = 2L), "positive integer")
  expect_error(rank_weight_tdnn_asymp(s1 = 2L, s2 = 0L, dim = 2L), "positive integer")
  expect_error(rank_weight_tdnn_asymp(s1 = 2L, s2 = 4L, dim = 0L), "positive integer")
  # s1 >= s2
  expect_error(rank_weight_tdnn_asymp(s1 = 4L, s2 = 4L, dim = 2L), "s1 < s2")
  expect_error(rank_weight_tdnn_asymp(s1 = 5L, s2 = 4L, dim = 2L), "s1 < s2")
  # non-integer scalars (coercion trap)
  expect_error(rank_weight_tdnn_asymp(s1 = 1.5, s2 = 4L, dim = 2L), "integer scalar")
  expect_error(rank_weight_tdnn_asymp(s1 = "foo", s2 = 4L, dim = 2L), "integer scalar")
  expect_error(rank_weight_tdnn_asymp(s1 = 2L, s2 = 3.9, dim = 2L), "integer scalar")
  expect_error(rank_weight_tdnn_asymp(s1 = 2L, s2 = 4L, dim = 1.5), "integer scalar")
  # Logical controls are not silently coerced via isTRUE().
  expect_error(rank_weight_tdnn_asymp(s1 = 2L, s2 = 4L, dim = 2L, normalize = NA), "normalize.*TRUE or FALSE")
  expect_error(rank_weight_tdnn_asymp(s1 = 2L, s2 = 4L, dim = 2L, normalize = "yes"), "normalize.*TRUE or FALSE")
  expect_error(
    rank_weight_tdnn_asymp(s1 = 2L, s2 = 4L, dim = 2L, normalize = c(TRUE, FALSE)),
    "normalize.*TRUE or FALSE"
  )
  # valid call must not error
  expect_no_error(rank_weight_tdnn_asymp(s1 = 3L, s2 = 6L, dim = 2L))
})

# ── [spec] ────────────────────────────────────────────────────────────────────

test_that("[spec] rank_weight_tdnn_asymp returns valid rdistnn_rank_weight spec", {
  rw <- rank_weight_tdnn_asymp(s1 = 3L, s2 = 6L, dim = 2L)
  expect_true(is_rank_weight(rw))
  expect_true(inherits(rw, "rdistnn_rank_weight"))
  expect_match(rw$name, "tdnn_asymp")
  expect_true(is.function(rw$eval_batch))
})

test_that("[spec] eval_batch produces correct shape and column sums", {
  .test_rng_scope({
    set.seed(5)
    n <- 14L; m <- 4L; s1 <- 2L; s2 <- 5L; d <- 2L
    R <- .make_rank_matrix(n, m, seed = 5L)
    rw <- rank_weight_tdnn_asymp(s1, s2, d)
    W  <- rank_weight_eval_batch(rw, R)
    expect_equal(dim(W), c(n, m))
    expect_equal(colSums(W), rep(1, m), tolerance = 1e-10)
  })
})

test_that("[spec] eval_batch has signed weights", {
  .test_rng_scope({
    set.seed(6)
    n <- 20L; m <- 3L; s1 <- 3L; s2 <- 8L; d <- 2L
    R <- .make_rank_matrix(n, m, seed = 6L)
    rw <- rank_weight_tdnn_asymp(s1, s2, d)
    W  <- rank_weight_eval_batch(rw, R)
    expect_true(any(W < 0))
  })
})

test_that("[spec] normalize=FALSE spec: colSums may deviate slightly from 1", {
  # With normalize=FALSE the rank table is not renormalized; floating-point
  # drift can cause sum != 1 exactly, but it should still be close.
  n <- 15L; m <- 2L; s1 <- 2L; s2 <- 5L; d <- 2L
  R <- .make_rank_matrix(n, m, seed = 9L)
  rw_norm <- rank_weight_tdnn_asymp(s1, s2, d, normalize = TRUE)
  rw_raw  <- rank_weight_tdnn_asymp(s1, s2, d, normalize = FALSE)
  W_norm <- rank_weight_eval_batch(rw_norm, R)
  W_raw  <- rank_weight_eval_batch(rw_raw,  R)
  # Normalized version sums exactly to 1
  expect_equal(colSums(W_norm), rep(1, m), tolerance = 1e-10)
  # Both have same sign pattern (the normalization only rescales, doesn't change signs)
  expect_equal(sign(W_norm), sign(W_raw))
})

# ── [parity] ──────────────────────────────────────────────────────────────────

test_that("[parity] asymp spec eval_batch agrees with rank tables at the requested ranks", {
  .test_rng_scope({
    set.seed(10)
    n <- 14L; m <- 3L; s1 <- 2L; s2 <- 5L; d <- 2L
    R <- .make_rank_matrix(n, m, seed = 10L)
    rw <- rank_weight_tdnn_asymp(s1, s2, d)
    W_spec <- rank_weight_eval_batch(rw, R)
    W_ref  <- sapply(seq_len(m), function(j)
      tdnn_rank_table_asymp(n, s1, s2, d)[R[, j]])
    expect_equal(W_spec, W_ref, tolerance = 1e-12)
  })
})

test_that("[parity] asymp and exact specs agree for large n (convergence)", {
  .test_rng_scope({
    set.seed(11)
    n <- 500L; m <- 2L; s1 <- 5L; s2 <- 12L; d <- 2L
    R <- .make_rank_matrix(n, m, seed = 11L)
    rw_exact <- rank_weight_tdnn_exact(s1, s2, d)
    rw_asymp <- rank_weight_tdnn_asymp(s1, s2, d)
    W_e <- rank_weight_eval_batch(rw_exact, R)
    W_a <- rank_weight_eval_batch(rw_asymp, R)
    # For large n the max column-wise difference should be small
    expect_lt(max(abs(colSums(W_e) - colSums(W_a))), 1e-8)
    expect_lt(max(abs(W_e - W_a)), 0.01)
  })
})
