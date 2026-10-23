test_that("frozen-rank loop and R fast paths agree for paired signals", {
  .with_localization_policy_warnings({
  old <- options(RdistNN.use_compiled = FALSE)
  on.exit(options(old), add = TRUE)
  for (n in c(2L, 4L, 17L)) for (m in c(1L, 3L)) for (s in unique(c(1L, min(2L, n - 1L), n - 1L))) {
    X <- cbind(seq_len(n) %% 3L, seq_len(n) %% 2L)
    X0 <- rbind(c(0, 0), c(.2, .5), c(2, 1))[seq_len(m), , drop = FALSE]
    a <- list(X = X, X0 = X0, numerator = sin(seq_len(n)),
              denominator = 1 + seq_len(n) / n, s = s,
              return = c("deletions", "pseudo_values"))
    if (s == 1L) {
      expect_warning(fast <- do.call(affine_ratio_dnn, c(a, list(jackknife = "fast"))), "uniform")
      expect_warning(loop <- do.call(affine_ratio_dnn, c(a, list(jackknife = "loo"))), "uniform")
    } else {
      fast <- do.call(affine_ratio_dnn, c(a, list(jackknife = "fast")))
      loop <- do.call(affine_ratio_dnn, c(a, list(jackknife = "loo")))
    }
    for (nm in c("theta", "deletions", "pseudo_values", "component_vcov", "theta_vcov", "status")) {
      expect_equal(fast[[nm]], loop[[nm]], tolerance = 1e-12)
    }
  }
  })
})

test_that("available compiled paired engine agrees with the R implementation", {
  .with_distance_tie_warning({
  old <- options(RdistNN.use_compiled = TRUE)
  on.exit(options(old), add = TRUE)
  skip_if(is.null(.rdistnn_compiled_fun("rustat_jackknife_rank_full")))
  for (n in c(2L, 4L, 17L, 37L)) for (m in c(1L, 3L)) for (s in unique(c(1L, min(11L, n - 1L), n - 1L))) {
    a <- list(
      X = cbind(seq_len(n) %% 5L, seq_len(n) / n),
      X0 = rbind(c(0, 0), c(2, .5), c(4, 1))[seq_len(m), , drop = FALSE],
      numerator = cos(seq_len(n)), denominator = 1 + sin(seq_len(n)) / 4,
      s = s, return = c("deletions", "pseudo_values")
    )
    options(RdistNN.use_compiled = TRUE)
    if (s == 1L) expect_warning(cpp <- do.call(affine_ratio_dnn, a), "uniform") else cpp <- do.call(affine_ratio_dnn, a)
    options(RdistNN.use_compiled = FALSE)
    if (s == 1L) expect_warning(r <- do.call(affine_ratio_dnn, a), "uniform") else r <- do.call(affine_ratio_dnn, a)
    for (nm in c("theta", "deletions", "pseudo_values", "component_vcov", "theta_vcov", "status")) {
      expect_equal(cpp[[nm]], r[[nm]], tolerance = 1e-11)
    }
  }
  })
})

test_that("custom geometry is frozen even for the explicit deletion loop", {
  calls <- 0L
  distance <- new_distance("sample_dependent", eval = function(X, x0, ctx = NULL) {
    calls <<- calls + 1L
    abs(X[, 1L] - x0[1L]) + seq_len(nrow(X)) / nrow(X)^2
  })
  args <- list(X = matrix(1:8, ncol = 1L), X0 = matrix(c(0, 5), ncol = 1L),
               numerator = sin(1:8), denominator = rep(1, 8), s = 3L,
               distance = distance, return = "deletions")
  fast <- do.call(affine_ratio_dnn, c(args, list(jackknife = "fast")))
  expect_equal(calls, 2L)
  calls <- 0L
  loop <- do.call(affine_ratio_dnn, c(args, list(jackknife = "loo")))
  expect_equal(calls, 2L)
  expect_equal(fast$deletions, loop$deletions, tolerance = 1e-12)
})

test_that("distinct-distance row permutations preserve pairing and IDs", {
  X <- matrix(c(.1, .7, .3, 1.2, -.4, .9), ncol = 1L)
  num <- c(2, 1, 4, 3, 5, 8)
  den <- c(1, -1, 2, 1, 1, 3)
  a <- affine_ratio_dnn(X = X, x0 = 0, numerator = num, denominator = den, s = 3L)
  p <- c(6, 2, 4, 3, 5, 1)
  b <- affine_ratio_dnn(X = X[p, , drop = FALSE], x0 = 0,
                        numerator = num[p], denominator = den[p], s = 3L, row_ids = p)
  expect_equal(a$theta, b$theta)
  expect_equal(a$theta_vcov, b$theta_vcov)
  expect_identical(b$row_ids, p)
})

test_that("large signal scales and per-target overflow are explicit", {
  .with_partial_result_warning({
  base <- affine_ratio_dnn(X = matrix(1:4, ncol = 1L), x0 = 0,
                           numerator = c(2, 1, 4, 3), denominator = c(1, -1, 2, 1), s = 2L)
  large <- affine_ratio_dnn(X = matrix(1:4, ncol = 1L), x0 = 0,
                            numerator = 1e150 * c(2, 1, 4, 3),
                            denominator = 1e150 * c(1, -1, 2, 1), s = 2L)
  expect_equal(large$theta_vcov, base$theta_vcov, tolerance = 1e-11)
  overflow <- affine_ratio_dnn(X = matrix(1:4, ncol = 1L), x0 = 0,
                               numerator = rep(1e308, 4), denominator = rep(1e-308, 4), s = 2L)
  expect_identical(overflow$status$ratio_status, "ratio_overflow")
  expect_true(is.na(overflow$theta))
  separate <- affine_ratio_dnn(
    X = matrix(1:4, ncol = 1L), X0 = matrix(c(0, 5), ncol = 1L),
    numerator = c(1e200, 0, 0, 0), denominator = rep(1, 4), s = 3L
  )
  expect_identical(separate$status$covariance_status, c("numerical_failure", "zero_score_variance"))
  expect_true(is.na(separate$se[1L]))
  expect_equal(separate$se[2L], 0)
  })
})
