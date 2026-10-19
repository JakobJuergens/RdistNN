test_that("scaled L1 recovers exceptional subtraction in every entry point", {
  old <- options(RdistNN.use_compiled = TRUE)
  on.exit(options(old), add = TRUE)
  expect_true(is.function(.rdistnn_compiled_fun("rustat_dist_l1_scaled_batch")))
  for (compiled in c(TRUE, FALSE)) {
    options(RdistNN.use_compiled = compiled)
    X <- matrix(c(8.000000000000001e307, 8e307))
    spec <- distance_l1_scaled("custom", 10)
    ctx <- spec$setup(X)
    D <- spec$eval_batch(X, matrix(rep(-7.9e307, 2)), ctx)
    expect_identical(D[, 1], c(1.5900000000000001e307, 1.59e307))
    expect_identical(order(D[, 1]), 2:1)
    expect_identical(D[, 1], spec$eval(X, -7.9e307, ctx))
    expect_identical(D[, 1], D[, 2])
    # Opposite finite inputs overflow IEEE subtraction; binary-halving first
    # gives the independent analytic result while retaining scaled rounding.
    expect_equal(dist_l1_scaled(matrix(1e308), -1e308, .25), 5e307, tolerance = 0)
    expect_identical(dist_l1_scaled(matrix(1e308), -1e308, 0), 0)
    expect_equal(dist_l1_scaled(matrix(1e308), -1e308, 1e-308), 2, tolerance = 4e-15)
    S <- distance_l1_scaled("custom", 4)
    expect_identical(S$eval_batch(matrix(1e308), matrix(c(-1e308, 1e308)), S$setup(matrix(1e308))),
      matrix(c(5e307, 0), 1))
    expect_identical(order_by_l1_scaled(matrix(1e308), -1e308, "custom", 4)$dx, 5e307)
    expect_identical(dist_l1(matrix(c(1.5e308, 1.5e308), 1), c(0, 0)), Inf)
    expect_identical(dist_l1_scaled(matrix(.Machine$double.xmax),
      -.Machine$double.xmax, 2^-1074), 1.7763568394002503e-15)
  }
})
