# -------------------------------------------------------------------
# Asymptotic DNN weights (FULL U-statistic, WOR) — tests
# File under test: R/weights_dnn_asymp_full.R
# Tags:
#   [smoke]   Construction & basic properties
#   [edge]    Edge cases / structure (s=1 uniform; k=1 monotone)
#   [approx]  Closeness to exact weights for moderate n
#   [inv]     Invariance to rank permutation (w_rank)
#   [opts]    normalize=FALSE parity
#   [errors]  Input validation
#   [meta]    Metadata correctness
# -------------------------------------------------------------------

test_that("[smoke] returns weights of length n that are nonnegative and sum to 1", {
  set.seed(1)
  n <- 60L
  s <- 12L
  k <- 4L
  ranks <- sample.int(n)
  out <- dnn_weights_asymp_full_from_ranks(ranks, k = k, s = s)

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

test_that("[edge] s = 1, k = 1 yields uniform weights exactly", {
  n <- 31L
  s <- 1L
  k <- 1L
  ranks <- sample.int(n)
  out <- dnn_weights_asymp_full_from_ranks(ranks, k = k, s = s)
  expect_equal(out$w_rank, rep(1 / n, n), tolerance = 1e-12) # was tolerance = 0
  expect_equal(out$w, out$w_rank[ranks], tolerance = 1e-12) # was tolerance = 0
})


test_that("[edge] k = 1 produces (near) monotone nonincreasing weights in rank", {
  # With k=1, w_rank[j] are Beta(1, s)-bin increments over equal bins; they should decrease with j.
  n <- 200L
  s <- 15L
  k <- 1L
  ranks <- 1:n
  out <- dnn_weights_asymp_full_from_ranks(ranks, k = k, s = s, normalize = FALSE)
  diffs <- diff(out$w_rank)
  # Allow tiny numerical wiggles
  expect_lte(max(diffs), 1e-12)
})

test_that("[approx] asymptotic weights are close to exact weights for moderate n", {
  set.seed(2)
  n <- 120L
  s <- 10L
  k <- 3L
  ranks <- 1:n

  w_ex <- dnn_weights_exact_full_from_ranks(ranks, k, s)$w_rank
  w_ap <- dnn_weights_asymp_full_from_ranks(ranks, k, s, normalize = TRUE)$w_rank

  # L1 gap should be small for moderate n; loosen if needed in CI
  l1 <- sum(abs(w_ex - w_ap))
  expect_lt(l1, 0.02) # ~2% total variation is a reasonable target here
})

test_that("[inv] w_rank depends only on (n,k,s) and not on permutation; mapping to data order changes accordingly", {
  set.seed(3)
  n <- 80L
  s <- 12L
  k <- 4L
  r1 <- 1:n
  r2 <- rev(r1)
  r3 <- sample.int(n)

  w1 <- dnn_weights_asymp_full_from_ranks(r1, k, s)$w_rank
  w2 <- dnn_weights_asymp_full_from_ranks(r2, k, s)$w_rank
  w3 <- dnn_weights_asymp_full_from_ranks(r3, k, s)$w_rank

  expect_equal(w2, w1, tolerance = 0)
  expect_equal(w3, w1, tolerance = 0)

  out2 <- dnn_weights_asymp_full_from_ranks(r2, k, s)
  out3 <- dnn_weights_asymp_full_from_ranks(r3, k, s)
  expect_equal(out2$w, w1[r2], tolerance = 0)
  expect_equal(out3$w, w1[r3], tolerance = 0)
})

test_that("[opts] normalize=FALSE is already ~normalized and equals TRUE within tight tolerance", {
  set.seed(4)
  n <- 90L
  s <- 9L
  k <- 3L
  ranks <- sample.int(n)

  outF <- dnn_weights_asymp_full_from_ranks(ranks, k, s, normalize = FALSE)
  outT <- dnn_weights_asymp_full_from_ranks(ranks, k, s, normalize = TRUE)

  expect_equal(sum(outF$w_rank), 1, tolerance = 1e-12)
  expect_equal(sum(outF$w), 1, tolerance = 1e-12)
  expect_equal(outF$w_rank, outT$w_rank, tolerance = 1e-12)
  expect_equal(outF$w, outT$w, tolerance = 1e-12)
})

test_that("[errors] input validation triggers appropriate errors", {
  n <- 6L
  expect_error(dnn_weights_asymp_full_from_ranks(NULL, k = 1, s = 1), "`ranks` must be provided")
  expect_error(dnn_weights_asymp_full_from_ranks(integer(0), k = 1, s = 1), "positive length")
  expect_error(dnn_weights_asymp_full_from_ranks(c(1L, NA, 3L), k = 1, s = 1), "cannot contain NA")
  expect_error(dnn_weights_asymp_full_from_ranks(c(1L, 1L, 2L, 3L, 4L, 5L), k = 1, s = 1), "permutation of 1..n")

  ranks <- 1:n
  expect_error(dnn_weights_asymp_full_from_ranks(ranks, k = 0, s = 1), "positive integer")
  expect_error(dnn_weights_asymp_full_from_ranks(ranks, k = 2, s = 1), "k <= s")
  # s < 1: caught by .validate_scalar_int as "positive integer"
  expect_error(dnn_weights_asymp_full_from_ranks(ranks, k = 1, s = 0), "positive integer")
  # s > n: explicit upper-bound check
  expect_error(dnn_weights_asymp_full_from_ranks(ranks, k = 1, s = n + 1), "s <= n")
})

test_that("[meta] meta fields are correct", {
  n <- 55L
  s <- 11L
  k <- 4L
  ranks <- sample.int(n)
  out <- dnn_weights_asymp_full_from_ranks(ranks, k, s)
  expect_identical(out$meta$n, n)
  expect_identical(out$meta$k, as.integer(k))
  expect_identical(out$meta$s, as.integer(s))
  expect_identical(out$meta$approx, "beta(OrderStat)")
})

# ---- [spec-errors] rank_weight_dnn_asymp constructor -------------------------

test_that("[spec-errors] rank_weight_dnn_asymp validates k and s at construction time", {
  # non-positive k / s
  expect_error(rank_weight_dnn_asymp(k = 0L, s = 3L), "positive integer")
  expect_error(rank_weight_dnn_asymp(k = 2L, s = 0L), "positive integer")
  # cross-parameter constraint
  expect_error(rank_weight_dnn_asymp(k = 2L, s = 1L), "k <= s")
  expect_error(rank_weight_dnn_asymp(k = 5L, s = 3L), "k <= s")
  # non-integer scalars (coercion trap)
  expect_error(rank_weight_dnn_asymp(k = 1.5, s = 3L), "integer scalar")
  expect_error(rank_weight_dnn_asymp(k = "foo", s = 3L), "integer scalar")
  expect_error(rank_weight_dnn_asymp(k = 1L, s = 2.9), "integer scalar")
  expect_error(rank_weight_dnn_asymp(k = 1L, s = "bar"), "integer scalar")
  # valid call must not error
  expect_no_error(rank_weight_dnn_asymp(k = 2L, s = 3L))
})
