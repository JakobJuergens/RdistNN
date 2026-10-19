test_that("exact rank tables use the same represented law in both engines", {
  old <- options(RdistNN.use_compiled = FALSE)
  on.exit(options(old))
  for (n in c(2L, 8L, 20L, 257L)) for (s in unique(c(2L, n - 1L, n))) {
    for (k in unique(c(1L, min(3L, s), s))) {
      options(RdistNN.use_compiled = FALSE)
      a <- rank_weight_eval(rank_weight_dnn_exact(k, s), seq_len(n))
      options(RdistNN.use_compiled = TRUE)
      cpp <- .rdistnn_compiled_fun("rustat_exact_rank_table")
      expect_true(is.function(cpp))
      if (!is.function(cpp)) next
      options(RdistNN.use_compiled = TRUE)
      b <- rank_weight_eval(rank_weight_dnn_exact(k, s), seq_len(n))
      expect_identical(a, b)
      expect_identical(a, cpp(n, k, s))
      expect_equal(.compensated_sum_R(a), 1, tolerance = 2e-15)
      if (k == s) expect_identical(a, rep(1 / n, n))
    }
  }
})

test_that("exact tables agree with independent subset enumeration", {
  for (n in 2:8) for (s in 1:n) for (k in unique(c(1L, s))) {
    subsets <- utils::combn(n, s)
    counts <- tabulate(as.vector(subsets[seq_len(k), , drop = FALSE]), nbins = n)
    w <- rank_weight_eval(rank_weight_dnn_exact(k, s), seq_len(n))
    expect_equal(w, counts / (ncol(subsets) * k), tolerance = 3e-15)
  }
})

test_that("scaled recurrences terminate tails without a subnormal plateau", {
  w <- rank_weight_eval(rank_weight_dnn_exact(1, 500), seq_len(2000))
  expect_true(all(w >= 0))
  expect_true(all(diff(w) <= 0))
  expect_true(all(w[1502:2000] == 0))
  # lchoose supplies an independent approximate tail oracle away from underflow.
  at <- c(1L, 50L, 500L, 1000L)
  oracle <- exp(lchoose(2000 - at, 499) - lchoose(2000, 500))
  expect_equal(w[at] / oracle, rep(1, length(at)), tolerance = 5e-12)
})

test_that("the public exact quantile pipeline keeps the represented boundary", {
  # The rational prefix is 9/10, but this computed table's rounded prefix is
  # one ULP below .9. Both hosts must select 19 for this represented table.
  outcome_order <- c(14, 2, 9, 7, 19, 5, 11, 20, 4, 6, 16, 10, 3, 17, 18, 15, 1, 12, 8, 13)
  old <- options(RdistNN.use_compiled = TRUE)
  on.exit(options(old))
  for (native in c(TRUE, FALSE)) {
    options(RdistNN.use_compiled = native)
    fit <- np_dnn(X = matrix(1:20), Y = order(outcome_order), x0 = 0, s = 2, se = FALSE)
    expect_identical(as.data.frame(diagnostics(fit))$weighted_y_q90, 19)
  }
})

test_that("near-full subsets retain their analytic near-uniform law", {
  n <- 2000L; k <- n - 2L; s <- n - 1L
  oracle <- c(rep((n - 1) / n / k, k), 1 / n, 0)
  old <- options(RdistNN.use_compiled = FALSE); on.exit(options(old))
  a <- rank_weight_eval(rank_weight_dnn_exact(k, s), seq_len(n))
  options(RdistNN.use_compiled = TRUE)
  b <- rank_weight_eval(rank_weight_dnn_exact(k, s), seq_len(n))
  expect_identical(a, b)
  expect_equal(a, oracle, tolerance = 3e-15)
})
