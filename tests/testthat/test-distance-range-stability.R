# Analytic 3-4-5 triangles and binary powers avoid overflow in the oracle itself.
test_that("L2 retains finite huge and tiny norms through scalar, batch and order", {
  old <- options(RdistNN.use_compiled = TRUE)
  on.exit(options(old), add = TRUE)
  expect_true(is.function(.rdistnn_compiled_fun("rustat_dist_l2_batch")))
  for (h in c(1e200, 1e-200, 2^510, 2^-540)) {
    unit <- rbind(c(3, 4), c(4, 3), c(0, 1), c(0, 0), c(0, 5))
    X <- unit * h
    X0 <- matrix(0, 3, 2)
    for (scaled in c(FALSE, TRUE)) {
      XX <- if (scaled) sweep(X, 2, c(2, 4), "*") else X
      spec <- if (scaled) distance_l2_scaled("custom", c(2, 4)) else distance_l2()
      ctx <- spec$setup(XX)
      reference <- c(5, 5, 1, 0, 5)
      previous <- NULL
      for (compiled in c(TRUE, FALSE)) {
        options(RdistNN.use_compiled = compiled)
        expect_identical(is.null(.rdistnn_compiled_fun("rustat_dist_l2_batch")), !compiled)
        D <- spec$eval_batch(XX, X0, ctx)
        expect_equal(D / h, matrix(reference, 5, 3), tolerance = 4e-15)
        expect_identical(order(D[, 1], seq_len(5)), c(4L, 3L, 1L, 2L, 5L))
        .expect_identical_fast(spec$eval(XX, X0[1, ], ctx), D[, 1])
        .expect_identical_fast(spec$eval_batch(XX[c(3, 1), , drop = FALSE], X0, ctx),
                         D[c(3, 1), , drop = FALSE])
        ordered <- if (scaled) order_by_l2_scaled(XX, X0[1, ], "custom", c(2, 4)) else order_by_l2(XX, X0[1, ])
        .expect_identical_fast(ordered$dx, D[, 1])
        # Every represented 5h distance tie preserves the input-row index.
        expect_identical(ordered$ord, c(4L, 3L, 1L, 2L, 5L))
        if (!is.null(previous)) .expect_identical_fast(D, previous)
        previous <- D
      }
    }
  }
})

test_that("L2 handles mixed scales, subnormal squared sums and true overflow", {
  old <- options(RdistNN.use_compiled = TRUE)
  on.exit(options(old), add = TRUE)
  huge <- 0x1.1ccf385ebc8ap+1023 # 1e308, with platform-independent input bits.
  for (compiled in c(TRUE, FALSE)) {
    options(RdistNN.use_compiled = compiled)
    X <- rbind(c(1e200, 1e-200, 0), c(1e-200, 0, 0), c(0, 0, 0))
    expect_equal(dist_l2(X, c(0, 0, 0)) / c(1e200, 1e-200, 1), c(1, 1, 0), tolerance = 0)
    # Each coordinate's square rounds to zero, but their exact sum is 2^-1073.
    tiny <- matrix(2^-538, 1, 8)
    expect_identical(dist_l2(tiny, rep(0, 8), squared = TRUE), 2^-1073)
    expect_equal(dist_l2(tiny, rep(0, 8)) / (sqrt(8) * 2^-538), 1, tolerance = 2e-15)
    expect_identical(dist_l2(matrix(1e200), 0, squared = TRUE), Inf)
    expect_identical(dist_l2(matrix(c(1.5e308, 1.5e308), 1), c(0, 0)), Inf)
    # Opposite signed finite coordinates may overflow subtraction; a small
    # inverse scale makes the distance finite and must avoid Inf * 0.
    expect_equal(dist_l2_scaled(matrix(huge), -huge, 0.25) / (huge / 2), 1, tolerance = 0)
    expect_identical(dist_l2_scaled(matrix(huge), -huge, 0), 0)
    expect_equal(dist_l2_scaled(matrix(huge), -huge, 1e-308), 2, tolerance = 4e-15)
    expect_identical(dist_l2_scaled(matrix(.Machine$double.xmax),
                     -.Machine$double.xmax, 2^-1074), 0x1.fffffffffffffp-50)
  }
})

test_that("extreme L2 fitted geometry survives serialization and target reordering", {
  old <- options(RdistNN.use_compiled = TRUE)
  on.exit(options(old), add = TRUE)
  # With n=4,s=2,k=1 the ordered weights are (1/2,1/3,1/6,0),
  # obtained by counting the six pairs; the nearest member supplies each value.
  for (h in c(1e200, 1e-200)) for (compiled in c(TRUE, FALSE)) {
    options(RdistNN.use_compiled = compiled)
    f <- np_dnn(X = matrix(c(3, 1, 0, 2) * h), Y = c(30, 10, 0, 20),
                x0 = 0, s = 2, k = 1, se = FALSE, return = character(),
                distance = distance_l2_scaled("custom", 2))
    expect_equal(f$mu, 20 / 3, tolerance = 2e-15)
    restored <- unserialize(serialize(f, NULL))
    # Predict with the other engine, preserving the original fitted scale.
    options(RdistNN.use_compiled = !compiled)
    expect_equal(predict(restored, X0 = matrix(c(4, 0, 4) * h)),
                 c(70 / 3, 20 / 3, 70 / 3), tolerance = 2e-15)
  }
})

test_that("mixed ordinary and underflowing coordinates retain cross-engine arithmetic", {
  old <- options(RdistNN.use_compiled = TRUE)
  on.exit(options(old), add = TRUE)
  # A tiny coordinate must not change how exceptional rows are selected in one
  # engine only; ordinary rows beside them keep their original arithmetic.
  X <- rbind(c(3, 4, 1e-200), c(4, 3, 1e-200), c(1e8, 1, 1e-200),
             c(3, 4, 0), c(4, 3, 0), c(0, 0, 0))
  for (squared in c(FALSE, TRUE)) {
    options(RdistNN.use_compiled = TRUE)
    native <- distance_l2(squared)$eval_batch(X, matrix(0, 2, 3))
    options(RdistNN.use_compiled = FALSE)
    fallback <- distance_l2(squared)$eval_batch(X, matrix(0, 2, 3))
    .expect_identical_fast(fallback, native)
    expect_equal(fallback[, 1], if (squared) c(25, 25, 1e16, 25, 25, 0) else c(5, 5, 1e8, 5, 5, 0), tolerance = 0)
    expect_identical(order(fallback[, 1]), c(6L, 1L, 2L, 4L, 5L, 3L))
  }
  # Independent literal ordinary reference: rounded products accumulated in
  # coordinate order. Include zeros so deferring underflow masks is exercised.
  set.seed(862)
  ordinary <- matrix(rnorm(3500), 500, 7)
  ordinary[seq(1, length(ordinary), by = 13)] <- 0
  expected <- numeric(nrow(ordinary))
  for (j in seq_len(ncol(ordinary))) {
    square <- ordinary[, j] * ordinary[, j]
    expected <- expected + square
  }
  for (compiled in c(TRUE, FALSE)) {
    options(RdistNN.use_compiled = compiled)
    .expect_identical_fast(dist_l2(ordinary, rep(0, 7), TRUE), expected)
    .expect_identical_fast(dist_l2(ordinary, rep(0, 7)), sqrt(expected))
  }
})

test_that("L2 uses the same exceptional range at the Stata finite boundary", {
  old <- options(RdistNN.use_compiled = TRUE)
  on.exit(options(old), add = TRUE)
  # Both squared coordinates are finite, but their sum crosses 2^1023.
  # The common scaled recurrence preserves the strict real-valued ordering;
  # this is not a claim of correctly rounded exact norms at the last bit.
  # Hexadecimal literals preserve the audited binary fixture on every parser.
  X <- rbind(c(0x1.1a07eed3827dp+511, 0x1.2ae696863ed85p+511),
             c(0x1.9af4345292aaep+511, 0))
  for (compiled in c(TRUE, FALSE)) {
    options(RdistNN.use_compiled = compiled)
    D <- distance_l2()$eval_batch(X, matrix(0, 2, 2))
    .expect_identical_fast(D[, 1], c(0x1.9af4345292aafp+511, 0x1.9af4345292aaep+511))
    expect_identical(order(D[, 1]), 2:1)
    expect_identical(order_by_l2(X, c(0, 0))$ord, 2:1)
    .expect_identical_fast(D[, 1], D[, 2])
    squared <- distance_l2(TRUE)$eval_batch(X, matrix(0, 2, 2))
    expect_true(all(is.finite(squared)))
    expect_equal(squared / 1e308, D^2 / 1e308, tolerance = 5e-16)
    XX <- matrix(c(0x1.c7b1f3cac7434p+1022, 0x1.c7b1f3cac7433p+1022))
    S <- distance_l2_scaled("custom", 10)
    ctx <- S$setup(XX)
    target <- -0x1.c1ffba583e5f2p+1022
    DD <- S$eval_batch(XX, matrix(rep(target, 2)), ctx)
    .expect_identical_fast(DD[, 1], c(0x1.6a47127468a76p+1020, 0x1.6a47127468a75p+1020))
    expect_identical(order(DD[, 1]), 2:1)
    .expect_identical_fast(S$eval(XX, target, ctx), DD[, 1])
  }
})
