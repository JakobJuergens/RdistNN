test_that("zero-weight opposite-sign outliers do not create jackknife variation", {
  old <- options(RdistNN.use_compiled = FALSE)
  on.exit(options(old), add = TRUE)
  X <- matrix(seq_len(30L), ncol = 1L)
  for (compiled in c(FALSE, TRUE)) for (jackknife in c("fast", "loo")) {
    options(RdistNN.use_compiled = compiled)
    for (k in c(1L, 3L)) for (sign in c(-1, 1)) {
      y <- sign * c(rep(1e16, 29L), -1)
      np <- np_dnn(X = X, x0 = 0, Y = y, k = k, s = 10L, jackknife = jackknife)
      cate <- cate_dnn(X = X, x0 = 0, pseudo_outcome = y, k = k, s = 10L, jackknife = jackknife)
      ratio <- affine_ratio_dnn(X = X, x0 = 0, numerator = y,
                               denominator = rep(1, 30L), k = k, s = 10L, jackknife = jackknife,
                               return = c("deletions", "pseudo_values"))
      expect_equal(as.numeric(np$se), 0)
      expect_equal(as.numeric(cate$se), 0)
      expect_equal(ratio$se, 0)
      expect_identical(ratio$status$covariance_status, "zero_score_variance")
      expect_equal(as.numeric(ratio$deletions[, 1L, 1L]), rep(sign * 1e16, 30L))
      plqr <- plqr_dnn(X = X, x0 = 0, Y = rep(1, 30L), D = matrix(2 * y, ncol = 1L),
                       mD = matrix(0, 30L, 1L), g = rep(0, 30L), theta_eval = 0,
                       k = k, s = 10L, jackknife = jackknife)
      expect_equal(as.numeric(plqr$moment_se), 0)
    }
  }
})

test_that("centering uses each target's support and preserves cross-target covariance", {
  old <- options(RdistNN.use_compiled = FALSE)
  on.exit(options(old), add = TRUE)
  i <- seq_len(60L)
  ranks <- cbind(i, rev(i))
  storage.mode(ranks) <- "integer"
  base <- c(4 * seq_len(30L), -8 * seq_len(30L))
  shifted <- base + c(rep(1e16, 30L), rep(-1e16, 30L))
  rw <- rank_weight_dnn_exact(k = 3L, s = 40L)
  table <- RdistNN:::.rank_weight_table(rw, 59L)
  full <- RdistNN:::.rank_weight_table(rw, 60L)
  for (compiled in c(FALSE, TRUE)) for (naive in c(FALSE, TRUE)) {
    options(RdistNN.use_compiled = compiled)
    engine <- if (naive) rw else NULL
    a <- RdistNN:::.dnn_loo_centered(base, ranks, table, full, rw = engine)$loo
    b <- RdistNN:::.dnn_loo_centered(shifted, ranks, table, full, rw = engine)$loo
    expect_equal(stats::cov(a), stats::cov(b), tolerance = 1e-11)
  }
})

test_that("signed TDNN and exposed raw rank engines ignore zero-weight outliers", {
  old <- options(RdistNN.use_compiled = FALSE)
  on.exit(options(old), add = TRUE)
  X <- matrix(seq_len(30L), ncol = 1L)
  y <- c(rep(1e16, 29L), -1e307)
  ranks <- matrix(seq_len(30L), ncol = 1L)
  rw <- rank_weight_dnn_exact(k = 1L, s = 10L)
  table <- RdistNN:::.rank_weight_table(rw, 29L)
  full <- RdistNN:::.rank_weight_table(rw, 30L)
  for (compiled in c(FALSE, TRUE)) {
    options(RdistNN.use_compiled = compiled)
    raw <- jackknife_rank_full(ranks, matrix(y), full, table)$theta_loo
    expect_equal(diff(range(raw)), 0)
    raw <- RdistNN:::.dnn_loo_prefix(y, ranks, table)
    expect_equal(diff(range(raw)), 0)
    raw <- RdistNN:::.dnn_loo_loop(y, ranks, rw)
    expect_equal(diff(range(raw)), 0)
    for (jackknife in c("fast", "loo")) {
      tdnn <- np_tdnn(X = X, x0 = 0, Y = y, s1 = 10L, s2 = 15L,
                      jackknife = jackknife)
      expect_equal(as.numeric(tdnn$se), 0)
    }
  }
})

test_that("large-level variance matches the independent R/Stata golden fixture", {
  old <- options(RdistNN.use_compiled = FALSE)
  on.exit(options(old), add = TRUE)
  # Shifted deletions are (4,0,0,0): (3/4) * (9+1+1+1) = 9.
  y <- c(1e16 + c(4, 8, 12), -1)
  for (compiled in c(FALSE, TRUE)) for (jackknife in c("fast", "loo")) {
    options(RdistNN.use_compiled = compiled)
    fit <- np_dnn(X = matrix(1:4), x0 = 0, Y = y, s = 3L, jackknife = jackknife)
    expect_equal(as.numeric(fit$se)^2, 9, tolerance = 1e-12)
  }
})

test_that("support-aware shifts preserve local variation and deletion boundary ranks", {
  old <- options(RdistNN.use_compiled = FALSE)
  on.exit(options(old), add = TRUE)
  X <- matrix(seq_len(30L), ncol = 1L)
  for (compiled in c(FALSE, TRUE)) for (jackknife in c("fast", "loo")) {
    options(RdistNN.use_compiled = compiled)
    for (k in c(1L, 3L)) {
      y <- 4 * seq_len(30L)
      base <- np_dnn(X = X, x0 = 0, Y = y, k = k, s = 10L, jackknife = jackknife)
      shifted <- 1e16 + y
      shifted[30L] <- -1e308
      fit <- np_dnn(X = X, x0 = 0, Y = shifted, k = k, s = 10L, jackknife = jackknife)
      expect_equal(fit$se, base$se, tolerance = 1e-11)
    }
    y <- c(rep(2, 20L), 7, rep(-1, 9L))
    fast <- np_dnn(X = X, x0 = 0, Y = y, s = 10L)
    slow <- np_dnn(X = X, x0 = 0, Y = y, s = 10L, jackknife = "loo")
    expect_equal(fast$se, slow$se, tolerance = 1e-11)
    for (normalize in c(FALSE, TRUE)) {
      fast <- np_dnn(X = X, x0 = 0, Y = y, s = 10L, method = "asymp", normalize = normalize)
      slow <- np_dnn(X = X, x0 = 0, Y = y, s = 10L, method = "asymp", normalize = normalize,
                     jackknife = "loo")
      expect_equal(fast$se, slow$se, tolerance = 1e-11)
    }
  }
})

test_that("asymptotic tails remain active beyond the exact support cutoff", {
  old <- options(RdistNN.use_compiled = FALSE)
  on.exit(options(old), add = TRUE)
  X <- matrix(seq_len(10L), ncol = 1L)
  y <- c(0, 0, rep(1, 8L))
  # Deleting either zero gives (8/9)^9; deleting any one gives (7/9)^9.
  # The variance of these two/eight deletion groups is (36/25)*(A-B)^2.
  expected_var <- (36 / 25) * ((8 / 9)^9 - (7 / 9)^9)^2
  for (compiled in c(FALSE, TRUE)) for (jackknife in c("fast", "loo")) {
    options(RdistNN.use_compiled = compiled)
    exact <- np_dnn(X = X, x0 = 0, Y = y, s = 9L, jackknife = jackknife)
    expect_equal(as.numeric(exact$se), 0)
    for (normalize in c(FALSE, TRUE)) {
      fit <- np_dnn(X = X, x0 = 0, Y = y, s = 9L, method = "asymp",
                     normalize = normalize, jackknife = jackknife)
      expect_equal(as.numeric(fit$mu), (4 / 5)^9, tolerance = 1e-12)
      expect_equal(as.numeric(fit$se)^2, expected_var, tolerance = 1e-12)
    }
  }
})

test_that("unused observations are excluded before an overflowing subtraction", {
  old <- options(RdistNN.use_compiled = FALSE)
  on.exit(options(old), add = TRUE)
  X <- matrix(seq_len(30L), ncol = 1L)
  for (compiled in c(FALSE, TRUE)) for (jackknife in c("fast", "loo")) {
    options(RdistNN.use_compiled = compiled)
    for (sign in c(-1, 1)) {
      # Subtracting the active reference from the last row would overflow.
      y <- sign * c(rep(1e308, 29L), -1e308)
      fit <- np_dnn(X = X, x0 = 0, Y = y, s = 10L, jackknife = jackknife)
      expect_true(is.finite(fit$mu))
      expect_equal(as.numeric(fit$se), 0)
    }
  }
})
