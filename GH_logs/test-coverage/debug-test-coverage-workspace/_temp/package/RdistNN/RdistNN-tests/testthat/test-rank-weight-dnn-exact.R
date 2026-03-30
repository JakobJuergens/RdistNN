# -------------------------------------------------------------------
# Exact DNN weights (FULL U-statistic, WOR) — tests
# File under test: R/weights_dnn_exact_full.R
# Tags:
#   [smoke]     Construction & basic properties
#   [edge]      Edge cases with closed-form answers
#   [bf]        Brute-force enumeration parity (small n)
#   [inv]       Invariance to rank permutation (w_rank)
#   [errors]    Input validation
#   [meta]      Metadata correctness
# -------------------------------------------------------------------

test_that("[smoke] returns weights of length n that are nonnegative and sum to 1", {
  set.seed(1)
  n <- 12L
  s <- 5L
  k <- 3L
  ranks <- sample.int(n)
  out <- dnn_weights_exact_full_from_ranks(ranks, k = k, s = s)

  expect_type(out$w, "double")
  expect_type(out$w_rank, "double")
  expect_length(out$w, n)
  expect_length(out$w_rank, n)
  expect_true(all(out$w >= 0))
  expect_true(all(out$w_rank >= 0))

  expect_equal(sum(out$w_rank), 1, tolerance = 1e-12)
  expect_equal(sum(out$w), 1, tolerance = 1e-12)
  # mapping: w[i] = w_rank[ ranks[i] ]
  expect_equal(out$w, out$w_rank[ranks], tolerance = 0)
})

test_that("[edge] s = n gives uniform 1/k on top-k ranks; s = 1 gives uniform weights", {
  # Case 1: s = n ⇒ subset is entire sample; i-th NN in subset is the global i-th NN deterministically.
  # Hence w_rank[j] = 1/k for j <= k, else 0.
  for (k in 1:4) {
    n <- 7L
    s <- n
    ranks <- 1:n
    out <- dnn_weights_exact_full_from_ranks(ranks, k = k, s = s)
    expect_equal(out$w_rank, c(rep(1 / k, k), rep(0, n - k)), tolerance = 0)
    expect_equal(out$w, out$w_rank, tolerance = 0)
  }

  # Case 2: s = 1, k = 1 ⇒ subset picks one item uniformly; that single item is the 1st NN in the subset.
  # Hence w_rank is uniform 1/n.
  n <- 9L
  s <- 1L
  k <- 1L
  ranks <- sample.int(n)
  out <- dnn_weights_exact_full_from_ranks(ranks, k = k, s = s)
  expect_equal(out$w_rank, rep(1 / n, n), tolerance = 1e-12) # was tolerance = 0
  expect_equal(out$w, out$w_rank[ranks], tolerance = 1e-12) # was tolerance = 0
})

test_that("[edge] k = 1 matches the simple closed form C(n-j, s-1)/C(n, s)", {
  n <- 10L
  s <- 4L
  k <- 1L
  ranks <- 1:n
  out <- dnn_weights_exact_full_from_ranks(ranks, k = k, s = s)
  j <- seq_len(n)
  ref <- ifelse(j <= n - s + 1L,
    choose(n - j, s - 1L) / choose(n, s),
    0
  )
  expect_equal(out$w_rank, ref, tolerance = 1e-12)
})

test_that("[bf] brute-force enumeration over all subsets matches exact weights (small n)", {
  # Enumerate all subsets of size s; each subset gives 1/k mass to the k smallest global ranks present.
  n <- 8L
  s <- 3L
  k <- 2L
  ranks <- 1:n # using identity simplifies mapping; test invariance separately

  # exact
  ex <- dnn_weights_exact_full_from_ranks(ranks, k = k, s = s)$w_rank

  # brute-force
  subs <- utils::combn(n, s) # each column a subset of indices in data order
  acc <- numeric(n)
  for (col in seq_len(ncol(subs))) {
    S <- subs[, col]
    # find the k smallest global ranks in S
    local_ord <- order(ranks[S]) # positions within S
    S_topk <- S[local_ord[seq_len(k)]]
    acc[S_topk] <- acc[S_topk] + 1 / k
  }
  bf <- acc / ncol(subs)

  expect_equal(ex, bf, tolerance = 1e-12)
  expect_equal(sum(bf), 1, tolerance = 1e-12)
})

test_that("[inv] w_rank depends only on (n,k,s) and not on the particular permutation", {
  n <- 12L
  s <- 5L
  k <- 3L
  r1 <- 1:n
  r2 <- rev(r1)
  r3 <- sample.int(n)

  w1 <- dnn_weights_exact_full_from_ranks(r1, k, s)$w_rank
  w2 <- dnn_weights_exact_full_from_ranks(r2, k, s)$w_rank
  w3 <- dnn_weights_exact_full_from_ranks(r3, k, s)$w_rank

  expect_equal(w2, w1, tolerance = 0)
  expect_equal(w3, w1, tolerance = 0)

  # mapping back to data order changes accordingly
  out2 <- dnn_weights_exact_full_from_ranks(r2, k, s)
  out3 <- dnn_weights_exact_full_from_ranks(r3, k, s)
  expect_equal(out2$w, w1[r2], tolerance = 0)
  expect_equal(out3$w, w1[r3], tolerance = 0)
})

test_that("[errors] input validation triggers appropriate errors", {
  n <- 6L
  expect_error(dnn_weights_exact_full_from_ranks(NULL, k = 1, s = 1), "`ranks` must be provided")
  expect_error(dnn_weights_exact_full_from_ranks(integer(0), k = 1, s = 1), "positive length")
  expect_error(dnn_weights_exact_full_from_ranks(c(1L, NA, 3L), k = 1, s = 1), "cannot contain NA")
  expect_error(dnn_weights_exact_full_from_ranks(c(1L, 1L, 2L, 3L, 4L, 5L), k = 1, s = 1), "permutation of 1..n")

  ranks <- 1:n
  expect_error(dnn_weights_exact_full_from_ranks(ranks, k = 0, s = 1), "positive integer")
  expect_error(dnn_weights_exact_full_from_ranks(ranks, k = 2, s = 1), "k <= s")
  expect_error(dnn_weights_exact_full_from_ranks(ranks, k = 1, s = 0), "1 <= s <= n")
  expect_error(dnn_weights_exact_full_from_ranks(ranks, k = 1, s = n + 1), "1 <= s <= n")
})

test_that("[meta] meta fields are correct", {
  n <- 9L
  s <- 4L
  k <- 2L
  ranks <- sample.int(n)
  out <- dnn_weights_exact_full_from_ranks(ranks, k, s)
  expect_identical(out$meta$n, n)
  expect_identical(out$meta$k, as.integer(k))
  expect_identical(out$meta$s, as.integer(s))
})
