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
