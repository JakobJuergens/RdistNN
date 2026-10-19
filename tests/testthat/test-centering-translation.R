test_that("nonzero remote tails do not destroy local jackknife variation", {
  old <- options(RdistNN.use_compiled = FALSE)
  on.exit(options(old))
  n <- 50L
  y <- 1e16 + rep(c(0, 2), 25)
  y[n] <- -1e16
  oracle <- 0.98398589833792965115112634738713
  for (native in c(FALSE, TRUE)) {
    options(RdistNN.use_compiled = native)
    if (native) expect_true(is.function(.rdistnn_compiled_fun("rustat_jackknife_rank_full")))
    for (jackknife in c("fast", "loo")) for (shift in c(0, -1e16)) {
      z <- y + shift
      a <- list(X = matrix(seq_len(n)), x0 = 0, s = 49L,
        method = "asymp", jackknife = jackknife)
      fits <- list(do.call(np_dnn, c(a, list(Y = z))),
        do.call(cate_dnn, c(a, list(pseudo_outcome = z))),
        do.call(affine_ratio_dnn, c(a, list(numerator = z, denominator = rep(1, n)))))
      for (fit in fits) expect_equal(as.numeric(fit$se), oracle, tolerance = 2e-14)
    }
  }
})

test_that("ordinary SD and covariance retain representable local differences", {
  old <- options(RdistNN.use_compiled = TRUE)
  on.exit(options(old), add = TRUE)
  x <- c(0, 2, 0, 2)
  expect_equal(.sample_sd(1e16 + x), sqrt(4 / 3), tolerance = 1e-15)
  expect_equal(.covariance_matrix(cbind(1e16 + x, 7 + 3 * x), factor = 1 / 3),
    unname(stats::cov(cbind(x, 3 * x))), tolerance = 1e-14)
  X <- cbind(x, c(0, 0, 1, 1))
  for (native in c(FALSE, TRUE)) {
    options(RdistNN.use_compiled = native)
    a <- np_dnn(X = X, Y = c(10, 20, 30, 40), x0 = c(0, .9), s = 3,
      distance = distance_l2_scaled("sd"), se = FALSE, return = "ranks")
    b <- np_dnn(X = sweep(X, 2, c(1e16, 0), "+"), Y = c(10, 20, 30, 40),
      x0 = c(1e16, .9), s = 3, distance = distance_l2_scaled("sd"),
      se = FALSE, return = "ranks")
    expect_identical(a$ranks, b$ranks)
    expect_equal(a$mu, b$mu, tolerance = 0)
  }
})

test_that("local reference and inference survive a remote tail moved to the first row", {
  n <- 50L
  y <- 1e16 + rep(c(0, 2), 25); y[n] <- -1e16
  old <- options(RdistNN.use_compiled = FALSE); on.exit(options(old))
  for (native in c(FALSE, TRUE)) for (jk in c("fast", "loo")) {
    options(RdistNN.use_compiled = native)
    for (ix in list(c(50L, 1:49), 50:1, c(2:50, 1L))) {
      f <- np_dnn(X = matrix((1:n)[ix]), Y = y[ix], x0 = 0, s = 49,
        method = "asymp", jackknife = jk)
      expect_equal(f$se, .983985898337929651, tolerance = 2e-14)
    }
  }
})

test_that("covariance centering does not erase separately represented components", {
  X <- cbind(c(2^55, 3, -2^55, 17), c(1, -1, 1, -1))
  old <- options(RdistNN.use_compiled = TRUE)
  on.exit(options(old), add = TRUE)
  for (native in c(TRUE, FALSE)) {
    options(RdistNN.use_compiled = native)
    for (ix in list(1:4, 4:1, c(2, 4, 1, 3))) {
      expect_identical(.covariance_matrix(X[ix, ])[1, 2], -20)
    }
    # The rounded centered leading values alone sum to -22. The residual
    # from mean subtraction must be included in the covariance correction.
    centered <- cbind(c(2^55 - 4, -2, -2^55 - 8, 12), X[, 2])
    residual <- cbind(c(-1, 0, 3, 0), rep(0, 4))
    expect_identical(.covariance_crossprod(centered, residual)[1, 2], -20)
    expect_identical(.covariance_crossprod(X, NULL)[1, 2], -20)
    scaled <- sweep(X, 2L, c(2^454, 2^-450), "*")
    expect_identical(.covariance_matrix(scaled)[1, 2], -320)
  }
})
