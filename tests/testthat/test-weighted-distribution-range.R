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
