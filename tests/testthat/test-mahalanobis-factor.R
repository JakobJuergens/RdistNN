test_that("computed negative extremes retain finite rank-one geometry", {
  old <- options(RdistNN.use_compiled = TRUE)
  on.exit(options(old), add = TRUE)
  X <- cbind(c((-5e307) * 2, 2e-300, 3e-300), 0)
  targets <- rbind(c(0, 0), X[1L, ])
  # The tiny local coordinates underflow only after division by the finite
  # sample SD. The distant row still has the independent distance sqrt(3).
  expected <- rbind(c(sqrt(3), 0), c(0, sqrt(3)), c(0, sqrt(3)))
  expected_ranks <- cbind(c(3L, 1L, 2L), c(1L, 2L, 3L))
  for (native in c(TRUE, FALSE)) {
    options(RdistNN.use_compiled = native)
    if (native) expect_true(is.function(.rdistnn_compiled_fun("rustat_dist_l2_batch")))
    for (metric in list(distance_l2_scaled(), distance_mahalanobis())) {
      ctx <- metric$setup(X)
      if (metric$name == "mahalanobis") expect_identical(ctx$rank, 1L)
      expect_equal(metric$eval_batch(X, targets, ctx), expected, tolerance = 2e-14)
      weights <- weights_dnn(X, X0 = targets, k = 1L, s = 3L,
        distance = metric, return = c("W", "ranks"))
      expect_identical(weights$ranks, expected_ranks)
      fit <- np_dnn(X = X, X0 = targets, Y = 1:3, k = 1L, s = 3L,
        distance = metric, se = FALSE)
      expect_identical(as.numeric(fit$mu), c(2, 1))
      restored <- unserialize(serialize(fit, NULL))
      expect_identical(as.numeric(predict(restored, X0 = targets)), c(2, 1))
      at <- c(3L, 1L, 2L)
      expect_equal(metric$eval_batch(X[at, , drop = FALSE], targets,
        metric$setup(X[at, , drop = FALSE])), expected[at, ], tolerance = 2e-14)
    }
  }
})

test_that("whitening projections retain small contributions across row shapes", {
  X <- rbind(c(2^54, 1, -2^54), c(1, 2, 3))
  # A rank-one whitening map with three equal loadings has an exact unit
  # projection for the cancellation row, independently of its large null part.
  ctx <- list(version = 2L, scale = 1, origin = c(0, 0, 0),
    factor = matrix(c(1, 1, 1, rep(0, 6)), 3, 3), rank = 1L)
  expected <- cbind(c(1, 6), 0, 0)
  old <- options(RdistNN.use_compiled = TRUE)
  on.exit(options(old), add = TRUE)
  for (native in c(TRUE, FALSE)) {
    options(RdistNN.use_compiled = native)
    expect_identical(.mahal_project(X, ctx), expected)
    for (i in 1:2) expect_identical(.mahal_project(X[i, , drop = FALSE], ctx),
      expected[i, , drop = FALSE])
  }
})

test_that("fresh Mahalanobis origins preserve small neighbors beside a distant observation", {
  old <- options(RdistNN.use_compiled = TRUE)
  on.exit(options(old), add = TRUE)
  X <- matrix(c(2^100, 1, 2, 3))
  targets <- matrix(c(2.1, 1.1, 2.1))
  y <- c(100, 10, 20, 30)
  # In one dimension the covariance metric has this direct scalar formula.
  expected <- abs(outer(X[, 1L], targets[, 1L], "-")) / stats::sd(X[, 1L])
  expected_ranks <- apply(expected, 2L, function(d) order(order(d)))
  expected_fit <- c(20, 10, 20)
  for (native in c(TRUE, FALSE)) {
    options(RdistNN.use_compiled = native)
    if (native) expect_false(is.null(.rdistnn_compiled_fun("rustat_dist_l2_batch")))
    metric <- distance_mahalanobis()
    ctx <- metric$setup(X)
    value <- metric$eval_batch(X, targets, ctx)
    # Relative comparisons preserve the tiny, positive distances in the oracle.
    expect_equal(value / expected, matrix(1, 4, 3), tolerance = 1e-14)
    for (j in seq_len(nrow(targets))) {
      expect_identical(metric$eval(X, targets[j, ], ctx), value[, j])
    }
    weights <- weights_dnn(X, X0 = targets, k = 1L, s = 4L,
      distance = metric, return = c("W", "ranks"))
    expect_identical(weights$ranks, expected_ranks)
    fit <- np_dnn(X = X, X0 = targets, Y = y, k = 1L, s = 4L,
      distance = metric, se = FALSE)
    expect_identical(as.numeric(fit$mu), expected_fit)
    restored <- unserialize(serialize(fit, NULL))
    expect_identical(as.numeric(predict(restored, X0 = targets[c(2, 1, 2), , drop = FALSE])),
      expected_fit[c(2, 1, 2)])
    at <- c(2:4, 1L)
    reordered <- metric$eval_batch(X[at, , drop = FALSE], targets,
      metric$setup(X[at, , drop = FALSE]))
    expect_equal(reordered[order(at), ] / expected, matrix(1, 4, 3), tolerance = 1e-14)
  }
})

test_that("fresh Mahalanobis origins retain ordinary variation around large offsets", {
  metric <- distance_mahalanobis()
  X <- matrix(2^100 + (0:3) * 2^49)
  target <- 2^100 + 5 * 2^48
  expected <- abs(X[, 1L] - target) / stats::sd(X[, 1L] - X[1L, 1L])
  expect_equal(metric$eval(X, target, metric$setup(X)), expected, tolerance = 1e-14)
})

test_that("Mahalanobis rank cutoff follows original covariate scales", {
  a <- c(-1, -1, 1, 1)
  b <- c(-1, 1, -1, 1)
  spec <- distance_mahalanobis()
  old <- options(RdistNN.use_compiled = TRUE)
  on.exit(options(old), add = TRUE)
  for (compiled in c(TRUE, FALSE)) {
    options(RdistNN.use_compiled = compiled)
    expect_identical(is.null(.rdistnn_compiled_fun("rustat_dist_l2_batch")), !compiled)
    for (gap in c(48, 1200)) {
      X <- cbind(a * 2^600, b * 2^(600 - gap))
      ctx <- spec$setup(X)
      # For n=4 the raw-design singular-value cutoff is 2^-50 relative
      # to the largest value. Gap 48 is retained; gap 1200 is discarded.
      expect_identical(ctx$rank, if (gap == 48) 2L else 1L)
      expected <- if (gap == 48) c(0, sqrt(3), sqrt(3), sqrt(6)) else
        c(0, 0, sqrt(3), sqrt(3))
      expect_equal(spec$eval(X, X[1, ], ctx), expected, tolerance = 2e-14)
    }
  }
})

test_that("Mahalanobis factor preserves resolvable near-collinear geometry", {
  x <- c(-2, -1, 0, 1, 2)
  z <- c(1, -1, 0, -1, 1)
  delta <- 2^-22
  X <- cbind(x, x + delta * z)
  a <- -17 / 8 + c(-1, 1) * 2^-10
  b <- -1 / 8
  X0 <- cbind(a, a + delta * b)
  expected <- vapply(a, function(target) sqrt((2 / 5) * (x - target)^2 + (z - b)^2), numeric(5))
  spec <- distance_mahalanobis()
  ctx <- spec$setup(X)
  expect_identical(ctx$rank, 2L)
  old <- options(RdistNN.use_compiled = TRUE)
  on.exit(options(old), add = TRUE)
  previous <- NULL
  for (compiled in c(TRUE, FALSE)) {
    options(RdistNN.use_compiled = compiled)
    expect_identical(is.null(.rdistnn_compiled_fun("rustat_dist_l2_batch")), !compiled)
    D <- spec$eval_batch(X, X0, ctx)
    if (!is.null(previous)) expect_identical(D, previous)
    previous <- D
    expect_equal(D, expected, tolerance = 2e-8)
    expect_identical(apply(D, 2, which.min), 1:2)
    expect_identical(spec$eval(X, X0[1, ], ctx), D[, 1])
    expect_identical(spec$eval_batch(X[c(3, 1), ], X0[c(2, 1, 2), ], ctx),
      D[c(3, 1), c(2, 1, 2)])
    fit <- np_dnn(X = X, Y = c(0, 100, 0, 0, 0), X0 = X0, s = 5,
      distance = spec, se = FALSE)
    expect_identical(as.numeric(fit$mu), c(0, 100))
    restored <- unserialize(serialize(fit, NULL))
    expect_identical(as.numeric(predict(restored, X0 = X0[c(2, 1, 2), ])), c(100, 0, 100))
  }
})

test_that("exceptional projection preserves tiny target differences beside offsets", {
  X <- cbind(2^50 + c(0, 1, 0, -1, 0, 0),
    c(0, 0, 2^-40, 0, -2^-40, 2^-1074))
  fit <- np_dnn(X = X, Y = c(10, 0, 0, 0, 0, 99), x0 = X[6, ],
    s = 6, se = FALSE, distance = distance_mahalanobis())
  expect_identical(as.numeric(fit$mu), 99)
  spec <- distance_mahalanobis()
  D <- spec$eval(X, X[6, ], spec$setup(X))
  expect_identical(D[6], 0)
  expect_true(D[1] > 0)
})

test_that("Mahalanobis factor handles extreme scales, zero rank and current contexts", {
  base <- rbind(c(2, 1), c(2, -1), c(-2, 1), c(-2, -1), c(0, 0))
  spec <- distance_mahalanobis()
  expected <- sqrt(base[, 1]^2 / 4 + base[, 2]^2)
  for (scale in c(2^600, 2^-600)) {
    X <- base * scale
    ctx <- spec$setup(X)
    expect_equal(spec$eval(X, c(0, 0), ctx), expected, tolerance = 2e-14)
    expect_identical(spec$eval(X, X[2, ], ctx)[2], 0)
  }
  X <- cbind(-2:2, -2:2)
  ctx <- spec$setup(X)
  expect_identical(ctx$rank, 1L)
  expect_equal(spec$eval(X, c(1, -1), ctx), abs(-2:2) / sqrt(2.5), tolerance = 2e-14)
  expect_identical(spec$eval(matrix(7, 3, 2), c(8, -1), spec$setup(matrix(7, 3, 2))), rep(0, 3))
  expect_error(spec$eval(base, c(0, 0), list(Si = diag(c(.25, 1)))),
    "current.*whitening.*refit")
  expect_error(spec$setup(matrix(1, 1, 2)), "at least two")
})

test_that("Mahalanobis factor handles wide designs and rejects incomplete contexts", {
  X <- rbind(c(1, 0, 1, 0), c(0, 1, 0, 1), c(-1, -1, -1, -1))
  X0 <- matrix(c(.25, .125, .25, .125), 1)
  d <- sweep(X[, 1:2], 2L, X0[1, 1:2], "-")
  expected <- sqrt(4 / 3 * (d[, 1]^2 - d[, 1] * d[, 2] + d[, 2]^2))
  spec <- distance_mahalanobis()
  ctx <- spec$setup(X)
  expect_identical(ctx$rank, 2L)
  expect_equal(spec$eval_batch(X, X0, ctx)[, 1], expected, tolerance = 2e-14)
  expect_error(spec$eval_batch(X, X0, list()), "current.*whitening.*refit")
  constant <- matrix(2^-600, 3, 2)
  expect_identical(spec$eval(constant, c(2^1000, -2^1000), spec$setup(constant)), rep(0, 3))
})

test_that("Mahalanobis normalization preserves small variation beside huge constants", {
  spec <- distance_mahalanobis()
  X <- cbind(rep(2^600, 3), c(-1, 0, 1) * 2^-600)
  ctx <- spec$setup(X)
  expect_identical(ctx$rank, 1L)
  expect_equal(spec$eval(X, c(2^600, 0), ctx), c(1, 0, 1), tolerance = 2e-14)
  X <- cbind(c(-1, 0, 1) * 2^-600, 0)
  expect_equal(spec$eval(X, c(0, 2^1000), spec$setup(X)), c(1, 0, 1), tolerance = 2e-14)
  # Subtracting opposite huge operands may overflow even though their scaled
  # difference and the final distance are both finite.
  X <- cbind(c(-1, 0, 1) * 2^1023, 0)
  expect_equal(spec$eval(X, c(2^1023, 0), spec$setup(X)), c(2, 1, 0), tolerance = 2e-14)
  X <- cbind(c(-1, 1) * 2^-1000, 0)
  expect_equal(spec$eval(X, c(2^24, 0), spec$setup(X)) / 2^1023.5,
    c(1, 1), tolerance = 2e-14)
  X <- cbind(c(-1, 0, 1) * 2^-600, c(-1, 0, 1) * 2^-1074)
  expect_equal(spec$eval(X, c(0, 2^600), spec$setup(X)) / 2^726,
    c(1, 1, 1), tolerance = 2e-14)
})
