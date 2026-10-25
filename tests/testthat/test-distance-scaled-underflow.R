test_that("scaled norms retain collectively representable subnormal products", {
  eta <- .Machine$double.xmin * .Machine$double.eps
  X <- rbind(c(0, 0), c(eta, eta), c(2 * eta, 2 * eta), c(3 * eta, 3 * eta))
  old <- options(RdistNN.use_compiled = FALSE)
  on.exit(options(old), add = TRUE)
  for (compiled in c(FALSE, TRUE)) {
    options(RdistNN.use_compiled = compiled)
    if (compiled) {
      expect_true(is.function(.rdistnn_compiled_fun("rustat_dist_l1_scaled")))
      expect_true(is.function(.rdistnn_compiled_fun("rustat_dist_l2_scaled")))
    }
    for (metric in c("l1", "l2")) {
      spec <- if (metric == "l1") {
        distance_l1_scaled("custom", scale_vec = c(2, 2))
      } else {
        distance_l2_scaled("custom", scale_vec = c(2, 2))
      }
      # Exact binary inputs give coordinate products eta/2, eta and 3*eta/2.
      # Round the complete L1/L2 norm, retaining collective half-eta mass.
      want <- if (metric == "l1") c(0, 1, 2, 3) else c(0, 1, 1, 2)
      ctx <- spec$setup(X)
      expect_identical(spec$eval(X, c(0, 0), ctx) / eta, want)
      expect_identical(spec$eval_batch(X, matrix(0, 2L, 2L), ctx) / eta,
        matrix(want, nrow = 4L, ncol = 2L))
    }
    squared <- distance_l2_scaled("custom", scale_vec = c(2, 2), squared = TRUE)
    expect_identical(squared$eval(X, c(0, 0), squared$setup(X)), numeric(4L))
    # Restoration must round only once at the subnormal boundary. Sequential
    # large power chunks can erase the residual just above half of eta.
    above_half <- .5 + .Machine$double.eps / 2
    expect_identical(distance_l1_scaled()$eval(matrix(eta), 0, list(inv_s = above_half)), eta)
    expect_identical(distance_l2_scaled()$eval(matrix(eta), 0, list(inv_s = above_half)), eta)
    # Correct rounded L2 distances are (0,eta,eta), so rows two and three
    # share the second rank's 1/3 mass rather than forming a zero-distance block.
    fit <- .with_localization_policy_warnings(np_dnn(X = X[1:3, , drop = FALSE], Y = c(0, 0, 1),
      x0 = c(0, 0), distance = distance_l2_scaled("custom", c(2, 2)),
      k = 1L, s = 2L, se = FALSE, return = c("W", "dx")))
    expect_equal(as.numeric(fit$W), c(2 / 3, 1 / 6, 1 / 6), tolerance = 0)
    expect_equal(fit$mu, 1 / 6, tolerance = 0)
  }
})

test_that("scaled norm recovery preserves ordinary coordinate arithmetic", {
  X <- rbind(c(3, 4, -2), c(.7, -1.3, 2.1), c(-2.2, .1, .9))
  x0 <- c(.1, -.2, .3)
  inverse <- c(.3, .7, 1.1)
  l1 <- l2_squared <- numeric(nrow(X))
  for (j in seq_len(ncol(X))) {
    d <- (X[, j] - x0[j]) * inverse[j]
    l1 <- l1 + abs(d)
    square <- d * d
    l2_squared <- l2_squared + square
  }
  old <- options(RdistNN.use_compiled = FALSE)
  on.exit(options(old), add = TRUE)
  for (compiled in c(FALSE, TRUE)) {
    options(RdistNN.use_compiled = compiled)
    expect_identical(distance_l1_scaled()$eval(X, x0, list(inv_s = inverse)), l1)
    expect_identical(distance_l2_scaled()$eval(X, x0, list(inv_s = inverse)), sqrt(l2_squared))
    expect_identical(distance_l2_scaled(squared=TRUE)$eval(X, x0, list(inv_s = inverse)), l2_squared)
  }
})

test_that("negative raw differences select the common exceptional range", {
  # Decimal conversion near DBL_MAX depends on R's long-double precision.
  # Preserve the exact input bits used by the independent norm oracle.
  X <- matrix(c(-0x1.59ead5d73161dp+1023, -0x1.38424eecf014bp+1023,
    -0x1.834c34d2acddfp+1023), nrow = 1L)
  inverse <- rep(2^-700, 3L)
  # Independent exact-binary products with a 2000-digit square-root oracle.
  want <- as.numeric("0x1.2ef864d71ee05p+324")
  old <- options(RdistNN.use_compiled = FALSE)
  on.exit(options(old), add = TRUE)
  for (compiled in c(FALSE, TRUE)) {
    options(RdistNN.use_compiled = compiled)
    expect_identical(distance_l2_scaled()$eval(X, rep(0, 3L), list(inv_s = inverse)), want)
  }
})

test_that("one active scaled coordinate keeps direct product and square rounding", {
  eta <- .Machine$double.xmin * .Machine$double.eps
  cases <- list(
    list(raw = 3 * eta, inverse = as.numeric("0x1.5555555555556p-3")),
    list(raw = 5 * eta, inverse = .3),
    list(raw = 2^-537, inverse = as.numeric("0x1.6a09e667f3bccp-1")),
    list(raw = 2^-537, inverse = as.numeric("0x1.6a09e667f3bcdp-1")),
    list(raw = 2^-537, inverse = as.numeric("0x1.6a09e667f3bcep-1"))
  )
  old <- options(RdistNN.use_compiled = FALSE)
  on.exit(options(old), add = TRUE)
  for (compiled in c(FALSE, TRUE)) {
    options(RdistNN.use_compiled = compiled)
    for (case in cases) for (padded in c(FALSE, TRUE)) {
      x <- if (padded) matrix(c(0, case$raw, 0), nrow = 1L) else matrix(case$raw)
      inverse <- if (padded) c(1, case$inverse, 1) else case$inverse
      x0 <- rep(0, ncol(x))
      product <- abs(case$raw * case$inverse)
      expect_identical(distance_l1_scaled()$eval(x, x0, list(inv_s = inverse)), product)
      expect_identical(distance_l2_scaled()$eval(x, x0, list(inv_s = inverse)), product)
      expect_identical(distance_l2_scaled(squared=TRUE)$eval(x, x0, list(inv_s = inverse)), product * product)
    }
    huge <- matrix(c(0, 1.6e308, 0), nrow = 1L)
    target <- c(0, -1.6e308, 0)
    expect_identical(distance_l1_scaled()$eval(huge, target, list(inv_s = c(1, .25, 1))), 8e307)
    expect_identical(distance_l2_scaled()$eval(huge, target, list(inv_s = c(1, .25, 1))), 8e307)
  }
})

test_that("exceptional norms stay above each rounded coordinate magnitude", {
  eta <- .Machine$double.xmin * .Machine$double.eps
  X <- matrix(c(3 * eta, eta), nrow = 1L)
  inverse <- c(as.numeric("0x1.5555555555556p-3"), eta)
  old <- options(RdistNN.use_compiled = FALSE)
  on.exit(options(old), add = TRUE)
  for (compiled in c(FALSE, TRUE)) {
    options(RdistNN.use_compiled = compiled)
    expect_identical(distance_l1_scaled()$eval(X, c(0, 0), list(inv_s = inverse)), eta)
    expect_identical(distance_l2_scaled()$eval(X, c(0, 0), list(inv_s = inverse)), eta)
  }
})
