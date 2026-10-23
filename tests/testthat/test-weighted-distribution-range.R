test_that("weighted SD retains square roots of unrepresentable normalized weights", {
  # For two observations, SD = |x2-x1|*sqrt(w1*w2)/(w1+w2).
  # Here it rounds to one, although w2/sum(w) = 2^-2000 underflows.
  # Every input, both square-root factors and the SD are representable.
  for (at in list(1:2, 2:1)) {
    for (sign in c(-1, 1)) {
      expect_equal(.weighted_sd((sign * c(0, 2^1000))[at],
        c(2^1000, 2^-1000)[at]), 1, tolerance = 2e-15)
    }
  }
})

test_that("weighted outcome SD avoids squaring physical units", {
  old <- options(RdistNN.use_compiled = TRUE)
  on.exit(options(old), add = TRUE)
  for (native in c(TRUE, FALSE)) {
    options(RdistNN.use_compiled = native)
    expect_equal(.np_outcome_diagnostics(matrix(rep(1/8, 8)),
      rep(c(-1e155, 1e155), 4))$weighted_sd_y / 1e155, 1, tolerance = 1e-15)
    expect_equal(.np_outcome_diagnostics(matrix(rep(.25, 4)),
      1e16 + c(0, 2, 0, 2))$weighted_sd_y, 1, tolerance = 0)
    expect_equal(.np_outcome_diagnostics(matrix(c(.5, .5, 0)),
      c(-1e-200, 1e-200, 1e300))$weighted_sd_y / 1e-200, 1, tolerance = 1e-15)
    expect_warning(fit <- np_dnn(X = matrix(1:8), Y = rep(c(-1e155, 1e155), 4),
      x0 = 0, s = 1, se = FALSE), "uniform")
    expect_equal(diagnostics(unserialize(serialize(fit, NULL)))$weighted_sd_y / 1e155,
      1, tolerance = 1e-15)
  }
})

test_that("positive tiny-tail outcomes preserve local spread under row permutations", {
  old <- options(RdistNN.use_compiled = TRUE)
  on.exit(options(old), add = TRUE)
  orders <- list(1:3, 3:1, c(2L, 1L, 3L))
  for (native in c(TRUE, FALSE)) {
    options(RdistNN.use_compiled = native)
    if (native) expect_false(is.null(.rdistnn_compiled_fun("rustat_compensated_sum")))
    for (sign in c(1, -1)) {
      x <- c(sign * 2^100, 1, 2)
      w <- c(2^-250, .5, .5)
      # Centering at the known local midpoint keeps this independent reference
      # well conditioned. The positive tail contributes the represented 2^-50.
      centered <- x - 1.5
      expected <- sqrt(sum(w * centered^2) - sum(w * centered)^2)
      for (at in orders) {
        expect_equal(.weighted_sd(x[at], w[at]), expected, tolerance = 1e-15)
      }
    }
    X <- matrix(c(129, 0:128, 130:256), ncol = 1L)
    y <- rep(c(1, 2), length.out = 257L)
    y[1L] <- 2^100
    fit <- np_dnn(X = X, x0 = -1, Y = y, k = 1L, s = 128L,
      se = FALSE, return = "W")
    w <- fit$W[, 1L]
    expect_gt(w[1L], 0)
    z <- y - 1
    expected <- sqrt(sum(w * z^2) - sum(w * z)^2)
    expect_equal(diagnostics(fit)$weighted_sd_y, expected, tolerance = 1e-15)
    at <- c(2:257, 1L)
    reordered <- np_dnn(X = X[at, , drop = FALSE], x0 = -1,
      Y = y[at], k = 1L, s = 128L, se = FALSE)
    expect_equal(diagnostics(reordered)$weighted_sd_y, expected, tolerance = 1e-15)
    expect_equal(diagnostics(unserialize(serialize(fit, NULL)))$weighted_sd_y,
      expected, tolerance = 1e-15)
    negative <- y * 0 + rep(c(1e-300, 2e-300), length.out = length(y))
    negative[1L] <- (-5e307) * 2
    fit_negative <- np_dnn(X = X, x0 = -1, Y = negative, k = 1L, s = 128L,
      se = FALSE)
    expected_negative <- -(negative[1L] * sqrt(w[1L]))
    expect_equal(diagnostics(fit_negative)$weighted_sd_y / expected_negative,
      1, tolerance = 2e-14)
  }
})

test_that("physical radii reuse stable Euclidean distances in original coordinates", {
  old <- options(RdistNN.use_compiled = TRUE)
  on.exit(options(old), add = TRUE)
  for (native in c(TRUE, FALSE)) {
    options(RdistNN.use_compiled = native)
    for (scale in c(1e200, 1e-200)) {
      X <- matrix((1:8) * scale)
      fit <- cate_dnn(X = X, pseudo_outcome = 1:8, x0 = 0, s = 2,
        distance = "l1", se = FALSE)
      expect_equal(diagnostics(fit)$physical_radius_90 / scale, 6, tolerance = 1e-15)
      saved <- unserialize(serialize(fit, NULL))
      expect_equal(diagnostics(saved)$physical_radius_90 / scale, 6, tolerance = 1e-15)
      ratio <- affine_ratio_dnn(X = X, numerator = 1:8, denominator = rep(1, 8),
        x0 = 0, s = 2, distance = "l1", se = FALSE)
      expect_equal(diagnostics(unserialize(serialize(ratio, NULL)))$physical_radius_90 / scale,
        6, tolerance = 1e-15)
    }
  }
})

test_that("weighted quantile decisions use the canonical compensated CDF", {
  # All masses and prefixes below are exact binary fractions. The selected
  # observation at .75 is 3, including with zero mass at the first observation.
  old <- options(RdistNN.use_compiled = TRUE)
  on.exit(options(old), add = TRUE)
  for (native in c(TRUE, FALSE)) {
    options(RdistNN.use_compiled = native)
    if (native) expect_false(is.null(.rdistnn_compiled_fun("rustat_weighted_cdf_index")))
    w <- c(0, .25, .5, .25)
    for (p in c(0, .25, .75, 1)) {
      expect_identical(.dnn_weighted_quantile(1:4, w, p),
        c(1, 2, 3, 4)[match(p, c(0, .25, .75, 1))])
    }
    expect_identical(.dnn_weighted_quantile(c(3, 2, 2, 1), c(.25, .5, .25, 0), .75), 2)
    expect_true(is.na(.dnn_weighted_quantile(1:2, c(0, 0))))
    fixtures <- read.csv(test_path("fixtures", "weighted-cdf.csv"), stringsAsFactors = FALSE)
    for (case in split(fixtures, fixtures$case)) {
      expect_identical(.dnn_weighted_quantile(case$position,
        as.numeric(case$weight_hex), as.numeric(case$prob_hex[1L])),
        as.numeric(case$expected[1L]))
    }
  }
})

test_that("weighted empirical quantiles respect exact and machine-scale boundaries", {
  expect_equal(.dnn_weighted_quantile(1:10, rep(0.1, 10), 0.1), 1)
  expect_equal(.dnn_weighted_quantile(1:10, rep(0.1, 10), 0.5), 5)
  expect_equal(.dnn_weighted_quantile(1:10, rep(0.1, 10), 0.9), 9)

  ulp_below_half <- .Machine$double.eps / 4
  expect_equal(
    .dnn_weighted_quantile(1:2, c(0.5 - ulp_below_half, 0.5 + ulp_below_half), 0.5),
    2
  )
  expect_equal(
    .dnn_weighted_quantile(1:2, c(0.5 + ulp_below_half, 0.5 - ulp_below_half), 0.5),
    1
  )
  expect_true(is.na(.dnn_weighted_quantile(1:2, c(0.5, 0.5), Inf)))
  expect_true(is.na(.dnn_weighted_quantile(1:2, c(-0.5, 1.5), 0.5)))
  expect_true(is.na(.dnn_weighted_quantile(1:2, c(0, 0), 0.5)))
})
