test_that("fitted SD recovers finite scales after intermediate range failures", {
  .with_distance_tie_warning({
  for (h in c(1e154, 1e200, 1e-200)) {
    X <- cbind(c(-1, -1, 1, 1) * h, c(0, 1, 0, 1))
    for (builder in list(distance_l1_scaled, distance_l2_scaled)) {
      spec <- builder("sd", eps = 1e-250)
      ctx <- spec$setup(X)
      expect_equal((1 / ctx$inv_s)[1] / h, sqrt(4 / 3), tolerance = 5e-15)
      expect_true(all(is.finite(ctx$inv_s) & ctx$inv_s > 0))
      f <- np_dnn(X = X, Y = c(0, 100, 0, 0), x0 = c(0, 1),
        k = 1, s = 4, se = FALSE, distance = spec)
      # Reflection of the first coordinate gives two equal nearest
      # distances (rows 2 and 4). Size-n k=1 assigns half to each.
      expect_identical(as.numeric(f$W), c(0, .5, 0, .5))
      expect_equal(f$mu, mean(c(100, 0)), tolerance = 0)
      restored <- unserialize(serialize(f, NULL))
      expect_equal(predict(restored, X0 = rbind(c(0, 1), c(0, 1))), rep(mean(c(100, 0)), 2L), tolerance = 0)
    }
  }
  })
})

test_that("exceptional SD centering preserves small variation around huge offsets", {
  unit <- 2^650
  X <- matrix(2^700 + c(-1, -1, 1, 1) * unit)
  for (builder in list(distance_l1_scaled, distance_l2_scaled)) {
    ctx <- builder("sd")$setup(X)
    expect_equal((1 / ctx$inv_s) / unit, sqrt(4 / 3), tolerance = 5e-15)
  }
})

test_that("failed fitted scales never become zero inverse scales", {
  X <- matrix(c(-.Machine$double.xmax, .Machine$double.xmax))
  for (builder in list(distance_l1_scaled, distance_l2_scaled)) {
    expect_error(builder("range")$setup(X), "Fitted scales.*nonfinite")
    expect_error(builder("sd")$setup(X), "Fitted scales.*nonfinite")
    expect_error(builder("sd")$setup(matrix(1)), "at least two")
    for (eps in list(0, -1, NA_real_, Inf, c(1, 2))) {
      expect_error(builder("sd", eps = eps)$setup(matrix(1:3)), "eps.*positive finite")
    }
    expect_identical(builder("sd")$setup(matrix(rep(2, 4)))$inv_s, 1e12)
  }
})

test_that("SD recovery preserves the independently checked near-tie neighbor", {
  # Exact binary inputs isolate scale arithmetic from decimal-parser precision.
  X <- cbind(c(0, 0x1.7dddf6b095ff1p+511, -0x1.7dddf6b095ff1p+511, 0x1.ed596bd58630ap+511,
    -0x1.ed596bd58630ap+511, 0x1.46dc10de8429fp+512, -0x1.46dc10de8429fp+512, 0),
    c(0, 1, -3, 3, -3, 3, -3, 3))
  target <- c(0x1.ca3d8e6d80cbbp+509, 0x1.7adea375ab194p+0)
  spec <- distance_l2_scaled("sd")
  ctx <- spec$setup(X)
  .expect_identical_fast(ctx$inv_s[1], 1 / 0x1.e3040a1e6c41ap+511)
  old <- options(RdistNN.use_compiled = TRUE)
  on.exit(options(old), add = TRUE)
  for (compiled in c(TRUE, FALSE)) {
    options(RdistNN.use_compiled = compiled)
    dx <- spec$eval(X, target, ctx)
    expect_identical(order(dx)[1:2], c(2L, 1L))
    f <- np_dnn(X = X, Y = c(0, 100, rep(0, 6)), x0 = target,
      s = 8, se = FALSE, distance = spec)
    expect_equal(f$mu, 100, tolerance = 0)
  }
})
