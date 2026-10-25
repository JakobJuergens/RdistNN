test_that("scaled L1 recovers exceptional subtraction in every entry point", {
  old <- options(RdistNN.use_compiled = TRUE)
  on.exit(options(old), add = TRUE)
  expect_true(is.function(.rdistnn_compiled_fun("rustat_dist_l1_scaled_batch")))
  # Fix the binary inputs: R's decimal parser uses different intermediate
  # precision on macOS ARM. These are 1e308 and the near-tie 8e307 fixture.
  huge <- 0x1.1ccf385ebc8ap+1023
  target <- -0x1.c1ffba583e5f2p+1022
  for (compiled in c(TRUE, FALSE)) {
    options(RdistNN.use_compiled = compiled)
    X <- matrix(c(0x1.c7b1f3cac7434p+1022, 0x1.c7b1f3cac7433p+1022))
    spec <- distance_l1_scaled("custom", 10)
    ctx <- spec$setup(X)
    D <- spec$eval_batch(X, matrix(rep(target, 2)), ctx)
    expect_identical(D[, 1], c(0x1.6a47127468a76p+1020, 0x1.6a47127468a75p+1020))
    expect_identical(order(D[, 1]), 2:1)
    expect_identical(D[, 1], spec$eval(X, target, ctx))
    expect_identical(D[, 1], D[, 2])
    # Opposite finite inputs overflow IEEE subtraction; binary-halving first
    # gives the independent analytic result while retaining scaled rounding.
    expect_equal(distance_l1_scaled()$eval(matrix(huge), -huge, list(inv_s = .25)), huge / 2, tolerance = 0)
    expect_identical(distance_l1_scaled()$eval(matrix(huge), -huge, list(inv_s = 0)), 0)
    expect_equal(distance_l1_scaled()$eval(matrix(huge), -huge, list(inv_s = 1e-308)), 2, tolerance = 4e-15)
    S <- distance_l1_scaled("custom", 4)
    expect_identical(S$eval_batch(matrix(huge), matrix(c(-huge, huge)), S$setup(matrix(huge))),
      matrix(c(huge / 2, 0), 1))
    expect_identical(.distance_ranking(matrix(huge), -huge, distance_l1_scaled(rule = "custom", scale_vec = 4))$dx, huge / 2)
    expect_identical(dist_l1(matrix(c(1.5e308, 1.5e308), 1), c(0, 0)), Inf)
    expect_identical(distance_l1_scaled()$eval(matrix(.Machine$double.xmax), -.Machine$double.xmax, list(inv_s = 2^-1074)), 0x1.fffffffffffffp-50)
  }
})
