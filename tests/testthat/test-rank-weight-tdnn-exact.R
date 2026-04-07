# -------------------------------------------------------------------
# Exact TDNN rank-weight table and from_ranks function — tests
# Files under test: R/wrank-tdnn_exact.R, R/wrank-tdnn-wrapper.R
# Tags:
#   [smoke]    Construction & basic properties
#   [scalars]  Bias-correction scalars w1*, w2*
#   [table]    tdnn_rank_table_exact is w1*·DNN(s1) + w2*·DNN(s2)
#   [signed]   Negative weights exist for reasonable params
#   [golden]   Hand-computed small case
#   [from_ranks] dnn_weights_tdnn_exact_full_from_ranks: mapping & meta
#   [inv]      w_rank invariant to rank permutation
#   [spec]     rank_weight_tdnn_exact spec: class, caching, batch eval
#   [errors]   Input validation
# -------------------------------------------------------------------

# ---- helpers used throughout ----------------------------------------

.tdnn_exact_brute <- function(n, s1, s2, dim) {
  # Compute the TDNN exact rank table directly from the definition:
  #   w1* * w_DNN(s1) + w2* * w_DNN(s2)
  # using .dnn_rank_table_exact_k1 via dnn_weights_exact_full_from_ranks.
  sc <- RdistNN:::.tdnn_scalars(s1, s2, dim)
  t1 <- dnn_weights_exact_full_from_ranks(seq_len(n), k = 1L, s = s1)$w_rank
  t2 <- dnn_weights_exact_full_from_ranks(seq_len(n), k = 1L, s = s2)$w_rank
  sc["w1"] * t1 + sc["w2"] * t2
}

# ---- [smoke] ---------------------------------------------------------

test_that("[smoke] tdnn_rank_table_exact returns length-n numeric summing to 1", {
  n <- 20L; s1 <- 3L; s2 <- 7L; d <- 2L
  w <- tdnn_rank_table_exact(n, s1, s2, d)
  expect_type(w, "double")
  expect_length(w, n)
  expect_true(is.finite(sum(w)))
  expect_equal(sum(w), 1, tolerance = 1e-10)
})

test_that("[smoke] dnn_weights_tdnn_exact_full_from_ranks returns correct structure", {
  set.seed(1)
  n <- 15L; s1 <- 3L; s2 <- 6L; d <- 2L
  ranks <- sample.int(n)
  out <- dnn_weights_tdnn_exact_full_from_ranks(ranks, s1, s2, d)
  expect_type(out, "list")
  expect_named(out, c("w", "w_rank", "meta"))
  expect_length(out$w, n)
  expect_length(out$w_rank, n)
  expect_equal(sum(out$w), 1, tolerance = 1e-10)
  expect_equal(sum(out$w_rank), 1, tolerance = 1e-10)
})

# ---- [scalars] -------------------------------------------------------

test_that("[scalars] .tdnn_scalars produces w1+w2=1 and correct signs", {
  # For s1 < s2 and d > 0: (s1/s2)^(-2/d) > 1 => w1* < 0, w2* > 1
  for (params in list(c(2, 4, 2), c(3, 9, 3), c(5, 10, 1))) {
    s1 <- params[1]; s2 <- params[2]; d <- params[3]
    sc <- RdistNN:::.tdnn_scalars(s1, s2, d)
    expect_equal(sc[["w1"]] + sc[["w2"]], 1, tolerance = 1e-14)
    expect_lt(sc[["w1"]], 0)
    expect_gt(sc[["w2"]], 1)
  }
})

test_that("[scalars] .tdnn_scalars closed form for s2 = 4*s1, dim = 2", {
  # (s1 / 4s1)^(-2/2) = (1/4)^(-1) = 4
  # w1* = 1/(1-4) = -1/3, w2* = 4/3
  sc <- RdistNN:::.tdnn_scalars(s1 = 2L, s2 = 8L, dim = 2L)
  expect_equal(sc[["w1"]], -1 / 3, tolerance = 1e-14)
  expect_equal(sc[["w2"]],  4 / 3, tolerance = 1e-14)
})

# ---- [table] ---------------------------------------------------------

test_that("[table] tdnn_rank_table_exact equals w1*·DNN(s1) + w2*·DNN(s2)", {
  n <- 25L; s1 <- 4L; s2 <- 8L; d <- 3L
  w_tdnn <- tdnn_rank_table_exact(n, s1, s2, d)
  w_ref  <- .tdnn_exact_brute(n, s1, s2, d)
  # re-normalize ref to match any internal normalization guard
  w_ref  <- w_ref / sum(w_ref)
  expect_equal(w_tdnn, w_ref, tolerance = 1e-12)
})

# ---- [signed] --------------------------------------------------------

test_that("[signed] TDNN exact weights can be negative (bias-correction)", {
  n <- 30L; s1 <- 5L; s2 <- 10L; d <- 2L
  w <- tdnn_rank_table_exact(n, s1, s2, d)
  expect_true(any(w < 0))
  expect_equal(sum(w), 1, tolerance = 1e-10)
})

test_that("[signed] dnn_weights_tdnn_exact_full_from_ranks w has negative entries", {
  set.seed(2)
  n <- 20L; s1 <- 4L; s2 <- 8L; d <- 2L
  out <- dnn_weights_tdnn_exact_full_from_ranks(sample.int(n), s1, s2, d)
  expect_true(any(out$w < 0))
  expect_true(any(out$w_rank < 0))
})

# ---- [golden] --------------------------------------------------------

test_that("[golden] tdnn_rank_table_exact matches hand-computed values for small n", {
  # n=6, s1=2, s2=4, dim=2
  # w1* = 1/(1-(2/4)^(-1)) = 1/(1-2) = -1,  w2* = 2
  # DNN(s1=2): w[j] = C(6-j,1)/C(6,2) = (6-j)/15  for j<=5, else 0
  # DNN(s2=4): w[j] = C(6-j,3)/C(6,4) = C(6-j,3)/15 for j<=3, else 0
  # TDNN = -1 * DNN(2) + 2 * DNN(4)
  n <- 6L; s1 <- 2L; s2 <- 4L; d <- 2L
  dnn2 <- c(5, 4, 3, 2, 1, 0) / 15
  dnn4 <- c(choose(5,3), choose(4,3), choose(3,3), 0, 0, 0) / 15
  w_ref <- -1 * dnn2 + 2 * dnn4
  w_ref <- w_ref / sum(w_ref)
  w <- tdnn_rank_table_exact(n, s1, s2, d)
  expect_equal(w, w_ref, tolerance = 1e-12)
})

# ---- [from_ranks] ----------------------------------------------------

test_that("[from_ranks] w = w_rank[ranks] (data-order mapping)", {
  set.seed(3)
  n <- 18L; s1 <- 3L; s2 <- 6L; d <- 2L
  ranks <- sample.int(n)
  out <- dnn_weights_tdnn_exact_full_from_ranks(ranks, s1, s2, d)
  expect_equal(out$w, out$w_rank[ranks], tolerance = 0)
})

test_that("[from_ranks] meta fields are correct", {
  n <- 14L; s1 <- 2L; s2 <- 5L; d <- 3L
  ranks <- sample.int(n)
  out <- dnn_weights_tdnn_exact_full_from_ranks(ranks, s1, s2, d)
  expect_identical(out$meta$n, n)
  expect_identical(out$meta$s1, as.integer(s1))
  expect_identical(out$meta$s2, as.integer(s2))
  expect_identical(out$meta$dim, as.integer(d))
  expect_named(out$meta$scalars, c("w1", "w2"))
  expect_equal(out$meta$scalars[["w1"]] + out$meta$scalars[["w2"]], 1, tolerance = 1e-14)
})

# ---- [inv] -----------------------------------------------------------

test_that("[inv] w_rank depends only on (n, s1, s2, dim), not the permutation", {
  n <- 20L; s1 <- 3L; s2 <- 7L; d <- 2L
  r1 <- seq_len(n)
  r2 <- rev(r1)
  set.seed(4)
  r3 <- sample.int(n)
  w1 <- dnn_weights_tdnn_exact_full_from_ranks(r1, s1, s2, d)$w_rank
  w2 <- dnn_weights_tdnn_exact_full_from_ranks(r2, s1, s2, d)$w_rank
  w3 <- dnn_weights_tdnn_exact_full_from_ranks(r3, s1, s2, d)$w_rank
  expect_equal(w2, w1, tolerance = 0)
  expect_equal(w3, w1, tolerance = 0)
})

# ---- [spec] ----------------------------------------------------------

test_that("[spec] rank_weight_tdnn_exact returns a valid rank_weight spec", {
  rw <- rank_weight_tdnn_exact(s1 = 3L, s2 = 6L, dim = 2L)
  expect_true(is_rank_weight(rw))
  expect_true(inherits(rw, "rank_weight"))
  expect_match(rw$name, "tdnn_exact")
  expect_true(is.function(rw$eval_batch))
})

test_that("[spec] eval_batch produces correct shape and column sums", {
  set.seed(5)
  n <- 12L; m <- 4L; s1 <- 2L; s2 <- 5L; d <- 2L
  R <- apply(replicate(m, sample.int(n)), 2L, identity)
  storage.mode(R) <- "integer"
  rw <- rank_weight_tdnn_exact(s1, s2, d)
  W  <- rank_weight_eval_batch(rw, R)
  expect_equal(dim(W), c(n, m))
  expect_equal(colSums(W), rep(1, m), tolerance = 1e-10)
})

test_that("[spec] table cache: second call with same n reuses cached table", {
  s1 <- 3L; s2 <- 6L; d <- 2L
  rw <- rank_weight_tdnn_exact(s1, s2, d)
  n <- 10L; m <- 2L
  R <- apply(replicate(m, sample.int(n)), 2L, identity)
  storage.mode(R) <- "integer"
  W1 <- rank_weight_eval_batch(rw, R)
  W2 <- rank_weight_eval_batch(rw, R)  # second call, should hit cache
  expect_equal(W1, W2, tolerance = 0)
})

test_that("[spec] eval_batch agrees with dnn_weights_tdnn_exact_full_from_ranks", {
  set.seed(6)
  n <- 14L; m <- 3L; s1 <- 2L; s2 <- 5L; d <- 2L
  R <- apply(replicate(m, sample.int(n)), 2L, identity)
  storage.mode(R) <- "integer"
  rw <- rank_weight_tdnn_exact(s1, s2, d)
  W_spec <- rank_weight_eval_batch(rw, R)
  W_ref  <- sapply(seq_len(m), function(j)
    dnn_weights_tdnn_exact_full_from_ranks(R[, j], s1, s2, d)$w)
  expect_equal(W_spec, W_ref, tolerance = 1e-12)
})

# ---- [spec-errors] rank_weight_tdnn_exact constructor ------------------------

test_that("[spec-errors] rank_weight_tdnn_exact validates s1, s2, dim at construction time", {
  # non-positive arguments
  expect_error(rank_weight_tdnn_exact(s1 = 0L, s2 = 4L, dim = 2L), "positive integer")
  expect_error(rank_weight_tdnn_exact(s1 = 2L, s2 = 0L, dim = 2L), "positive integer")
  expect_error(rank_weight_tdnn_exact(s1 = 2L, s2 = 4L, dim = 0L), "positive integer")
  # s1 >= s2
  expect_error(rank_weight_tdnn_exact(s1 = 4L, s2 = 4L, dim = 2L), "s1 < s2")
  expect_error(rank_weight_tdnn_exact(s1 = 5L, s2 = 4L, dim = 2L), "s1 < s2")
  # non-integer scalars (coercion trap)
  expect_error(rank_weight_tdnn_exact(s1 = 1.5, s2 = 4L, dim = 2L), "integer scalar")
  expect_error(rank_weight_tdnn_exact(s1 = "foo", s2 = 4L, dim = 2L), "integer scalar")
  expect_error(rank_weight_tdnn_exact(s1 = 2L, s2 = 3.9, dim = 2L), "integer scalar")
  expect_error(rank_weight_tdnn_exact(s1 = 2L, s2 = 4L, dim = 1.5), "integer scalar")
  # valid call must not error
  expect_no_error(rank_weight_tdnn_exact(s1 = 3L, s2 = 6L, dim = 2L))
})

# ---- [errors] --------------------------------------------------------

test_that("[errors] tdnn_rank_table_exact validates inputs", {
  expect_error(tdnn_rank_table_exact(0L,  2L, 4L, 2L), "positive integer")
  expect_error(tdnn_rank_table_exact(10L, 4L, 4L, 2L), "s1 < s2")
  expect_error(tdnn_rank_table_exact(10L, 5L, 4L, 2L), "s1 < s2")
  expect_error(tdnn_rank_table_exact(10L, 2L, 11L, 2L), "s2 <= n")
  expect_error(tdnn_rank_table_exact(10L, 2L, 4L, 0L), "positive integer")
})

test_that("[errors] dnn_weights_tdnn_exact_full_from_ranks validates inputs", {
  n <- 8L
  expect_error(dnn_weights_tdnn_exact_full_from_ranks(NULL, 2L, 4L, 2L),
               "`ranks` must be provided")
  expect_error(dnn_weights_tdnn_exact_full_from_ranks(integer(0), 2L, 4L, 2L),
               "positive length")
  expect_error(dnn_weights_tdnn_exact_full_from_ranks(c(1L, NA, 3L), 1L, 2L, 2L),
               "cannot contain NA")
  expect_error(dnn_weights_tdnn_exact_full_from_ranks(c(1L, 1L, 2L, 3L), 1L, 2L, 2L),
               "permutation of 1..n")
  ranks <- seq_len(n)
  expect_error(dnn_weights_tdnn_exact_full_from_ranks(ranks, 4L, 4L, 2L), "s1 < s2")
  expect_error(dnn_weights_tdnn_exact_full_from_ranks(ranks, 2L, n + 1L, 2L), "s2 <= n")
  expect_error(dnn_weights_tdnn_exact_full_from_ranks(ranks, 2L, 5L, 0L), "positive integer")
})
