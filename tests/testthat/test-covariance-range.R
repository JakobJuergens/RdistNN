test_that("shared covariance recovers finite large and subnormal cross-products", {
  a <- 1.25 * 2^511
  X <- cbind(c(-a, a), c(a / 2, -a / 2), c(7, 7))
  expected <- outer(c(a, -a / 2, 0), c(a, -a / 2, 0))
  V <- .covariance_matrix(X, factor = .5)
  expect_equal(V / a^2, expected / a^2, tolerance = 2e-15)
  expect_identical(V, t(V))
  expect_identical(.covariance_matrix(X, factor = .5, diagonal = TRUE), diag(V))
  # The unscaled sum overflows IEEE double, although its half is finite.
  large <- 1.5 * 2^511
  expect_equal(.covariance_matrix(matrix(c(-large, large)), .5)[1, 1] / large^2, 1,
    tolerance = 2e-15)
  tiny <- matrix(rep(c(-1, 1), 4) * 2^-538)
  expect_identical(.covariance_matrix(tiny, 7 / 8)[1, 1], 2^-1073)
})

test_that("shared covariance preserves ordinary identities and adapter scaling", {
  .with_distance_tie_warning({
  X <- cbind(c(-2, -1, 0, 1, 2, 3), c(1, 2, 1, -1, 0, 2))
  Y <- c(1, -2, 0, 3, -1, 1)
  V <- .covariance_matrix(X, 5 / 6)
  centered <- sweep(X, 2, colMeans(X), "-")
  expect_equal(V, 5 / 6 * crossprod(centered), tolerance = 2e-15)
  expect_identical(diag(V), .covariance_matrix(X, 5 / 6, diagonal = TRUE))
  old <- options(RdistNN.use_compiled = TRUE)
  on.exit(options(old), add = TRUE)
  for (compiled in c(TRUE, FALSE)) {
    options(RdistNN.use_compiled = compiled)
    for (engine in c("np", "tdnn", "cate", "ratio")) {
      run <- function(y) switch(engine,
        np = np_dnn(X = X, Y = y, X0 = X[c(1, 4), ], s = 3),
        tdnn = np_tdnn(X = X, Y = y, X0 = X[c(1, 4), ], s1 = 2, s2 = 3),
        cate = cate_dnn(X = X, pseudo_outcome = y, X0 = X[c(1, 4), ], s = 3),
        ratio = affine_ratio_dnn(X = X, numerator = y, denominator = rep(2, 6),
          X0 = X[c(1, 4), ], s = 3))
      ordinary <- run(Y)
      for (h in c(2^510, 2^-510)) {
        fit <- run(Y * h)
        expect_equal(as.numeric(fit$se) / h, as.numeric(ordinary$se), tolerance = 3e-14)
        if (engine %in% c("np", "tdnn")) {
          expect_equal(vcov(fit) / h^2, vcov(ordinary), tolerance = 4e-14)
          expected <- predict(ordinary, X0 = X[c(3, 1, 3), ], se.fit = TRUE)
          actual <- predict(unserialize(serialize(fit, NULL)), X0 = X[c(3, 1, 3), ], se.fit = TRUE)
          expect_equal(actual$se.fit / h, expected$se.fit, tolerance = 3e-14)
        }
      }
    }
  }
  })
})

test_that("shared covariance distinguishes unavailable variance from determinism", {
  X <- cbind(c(-1, 1) * 2^-540, c(3, 3), c(-1, 1))
  expect_error(.covariance_matrix(X, .5), "[Cc]ovariance.*rescale")
  V <- .covariance_matrix(X, .5, allow_failure = TRUE)
  expect_true(all(is.na(V[1, ])))
  expect_true(all(is.na(V[, 1])))
  expect_identical(V[2:3, 2:3], diag(c(0, 1)))
  expect_identical(.covariance_matrix(matrix(rep(2^900, 3))), matrix(0))
  V <- matrix(c(2^1022, 2^1021, 2^1021, 2^1022), 2)
  expect_identical(.covariance_symmetrize(V), V)
})

test_that("NP fitted, joint and marginal covariance share range-safe arithmetic", {
  X <- matrix(c(0, 1))
  old <- options(RdistNN.use_compiled = TRUE)
  on.exit(options(old), add = TRUE)
  for (compiled in c(TRUE, FALSE)) {
    options(RdistNN.use_compiled = compiled)
    expect_identical(is.null(.rdistnn_compiled_fun("rustat_jackknife_rank_full")), !compiled)
    for (a in c(1.25, 1.5) * 2^511) {
      fit <- suppressWarnings(np_dnn(X = X, Y = c(-a, a), X0 = matrix(c(0, 1)), s = 1))
      expect_equal(fit$se / a, c(1, 1), tolerance = 2e-15)
      expect_equal(unname(vcov(fit)) / a^2, matrix(1, 2, 2), tolerance = 2e-15)
      expect_equal(predict(fit, X0 = matrix(c(1, 0, 1)), se.fit = TRUE)$se.fit / a,
        rep(1, 3), tolerance = 2e-15)
    }
    expect_error(suppressWarnings(np_dnn(X = X, Y = c(-1, 1) * 2^-540,
      x0 = 0, s = 1)), "[Cc]ovariance.*rescale")
  }
})

test_that("covariance range failure has a numeric identity and useful action", {
  error <- tryCatch(.covariance_matrix(matrix(c(-1, 1) * 2^-540), .5), error = identity)
  expect_s3_class(error, "rdistnn_numerical_error")
  expect_identical(error$code, "numerical_range")
  expect_match(conditionMessage(error), "rescale.*point estimates")
})

test_that("range normalization preserves finite small cross-covariances", {
  old <- options(RdistNN.use_compiled = TRUE)
  on.exit(options(old), add = TRUE)
  large <- 2^510
  cases <- list(
    list(x = cbind(c(large, -large, 2^-510, -2^-510, 0, 0),
                  c(0, 0, 2^-510, -2^-510, large, -large)), value = 2^-1019),
    list(x = cbind(c(large, -large, 2^-1000, -2^-1000, 0, 0),
                  c(0, 0, 2^500, -2^500, large, -large)), value = 2^-499),
    # A nonzero mean leaves two-difference centering residuals at the huge rows.
    list(x = cbind(c(large, -large, 2^-1000, 0, 0, 0),
                  c(0, 0, 2^500, -2^500, large, -large)), value = 2^-500)
  )
  for (compiled in c(TRUE, FALSE)) {
    options(RdistNN.use_compiled = compiled)
    for (case in cases) for (factor in c(1, .5, 5 / 6)) {
      V <- .covariance_matrix(case$x, factor)
      expect_equal(V[1, 2] / case$value, factor, tolerance = 3e-15)
      expect_identical(V, t(V))
      expect_identical(diag(V), .covariance_matrix(case$x, factor, diagonal = TRUE))
      expect_identical(.covariance_matrix(case$x[, 2:1], factor), V[2:1, 2:1])
    }
    # The ratio is undefined, but its documented paired component covariance
    # remains finite and must retain the small numerator/denominator entry.
    signal <- cbind(c(2^513, -2^513, 2^-507, -2^-507, 0, 0),
                    c(0, 0, 2^-507, -2^-507, 2^513, -2^513))
    fit <- NULL
    expect_warning(fit <- .with_localization_policy_warnings(affine_ratio_dnn(
      X = matrix(1:6), x0 = 0, s = 1L,
      numerator = signal[, 1], denominator = signal[, 2])),
      class = "rdistnn_partial_result_warning")
    expect_equal(fit$component_vcov[[1]][1, 2] / (2^-1014 / 15), 1,
      tolerance = 2e-14)
  }
})

test_that("original-scale covariance recovery preserves protected and ordinary cases", {
  # These products overflow before summation, while the normalized covariance
  # and its zero cross-entry are finite. Direct recovery must not replace them.
  large <- 1.5 * 2^512
  X <- cbind(c(large, -large, large, -large),
             c(large, large, -large, -large))
  V <- .covariance_matrix(X, 1 / 64)
  expect_identical(V[1, 2], 0)
  expect_equal(diag(V) / (large / 4)^2, c(1, 1), tolerance = 2e-15)
  tiny <- matrix(rep(c(-1, 1), 4) * 2^-538)
  expect_identical(.covariance_matrix(cbind(tiny, tiny), 7 / 8)[1, 2], 2^-1073)

  large <- 2^510
  X <- cbind(c(large, -large, 2^-1000, -2^-1000, 0, 0),
             c(0, 0, 2^500, 2^500, large, -large))
  expect_identical(.covariance_matrix(X)[1, 2], 0)
  units <- cbind(c(-2, -1, 0, 1, 2, 3), c(3, -1, 2, 0, -2, 1))
  expected <- .covariance_matrix(units, 5 / 6)
  offset <- 2^511 + units * 2^459
  expect_equal(.covariance_matrix(offset, 5 / 6) / 2^918, expected,
    tolerance = 3e-15)

  # Ordinary columns retain their single centering pass and do not prepare
  # the original-scale recovery data or perform a second dot-product scan.
  calls <- 0L
  original <- .numeric_signal_level
  local_mocked_bindings(.numeric_signal_level = function(...) {
    calls <<- calls + 1L
    original(...)
  }, .package = "RdistNN")
  .covariance_matrix(units)
  expect_identical(calls, 2L)
})

test_that("original-scale recovery preserves accurately normalized subnormal products", {
  withr::local_options(RdistNN.use_compiled = FALSE)
  smallest <- 2^-1074
  unit <- 2^-524
  a <- 1.1 * unit
  b <- unit * (1 - 2^-27)
  uncentered <- cbind(c(a, a), c(unit, -b))
  centered <- cbind(c(a, a, -a, -a), c(unit, -b, -unit, b))
  # Exact represented-input products sum to .55... and 1.1... times the
  # smallest double. Both round to that positive value. Direct subnormal
  # products round separately to equal opposite values and lose the sum.
  for (compiled in c(FALSE, TRUE)) {
    options(RdistNN.use_compiled = compiled)
    .expect_identical_fast(.covariance_matrix(uncentered, center = FALSE)[1L, 2L], smallest)
    .expect_identical_fast(.covariance_matrix(centered)[1L, 2L], smallest)
    expect_warning(fit <- .with_localization_policy_warnings(affine_ratio_dnn(
      X = matrix(1:4), x0 = 0, s = 1L,
      numerator = 3 * centered[, 1L], denominator = 3 * centered[, 2L])),
      class = "rdistnn_partial_result_warning")
    # For these actual represented signals the exact covariance is
    # 7430939385161319 / 9007199254740992 times the smallest double.
    .expect_identical_fast(fit$component_vcov[[1L]][1L, 2L], smallest)
  }
})
