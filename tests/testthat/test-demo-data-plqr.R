# PLQR demo-data oracle threshold density and its inference use.
test_that("plqr_demo_data ships a true oracle threshold density", {
  tau <- 0.3
  sigma <- 0.8
  demo <- plqr_demo_data(n = 60, p = 2, p_d = 2, tau = tau, sigma = sigma, seed = 5)

  expect_true(all(is.finite(demo$data$density_oracle)))
  expect_true(all(demo$data$density_oracle > 0))
  expect_equal(
    unique(demo$data$density_oracle),
    4 * tau * (1 - tau) / sigma,
    tolerance = 1e-12
  )
  expect_equal(demo$truth$density_fn(0), 4 * tau * (1 - tau) / sigma, tolerance = 1e-12)
  # The error density is continuous at the quantile threshold.
  expect_equal(
    demo$truth$density_fn(-1e-10),
    demo$truth$density_fn(1e-10),
    tolerance = 1e-6
  )

  fit <- plqr_dnn(
    data = demo$data,
    x_cols = c("X1", "X2"),
    D_cols = c("D1", "D2"),
    g_col = "g_oracle",
    mD_cols = c("mD1", "mD2"),
    y_col = "Y",
    density_col = "density_oracle",
    theta_eval = demo$truth$theta_true,
    tau = tau,
    x0 = c(0, 0),
    method = "exact",
    k = 1L, s = 8L,
    se = TRUE
  )
  expect_true(all(is.finite(fit$theta_se)))
  expect_identical(fit$meta$theta_se_method, "direct_sandwich")
})

test_that("one-coordinate PLQR bumps match scalar truth and restore seeded RNG", {
  .test_rng_scope({
    for (seeded in c(TRUE, FALSE)) {
      set.seed(276L)
      before <- .Random.seed
      if (!seeded) rm(".Random.seed", envir = .GlobalEnv)
      demo <- plqr_demo_data(n = 8L, p = 1L, design = "bump", x0 = .2,
        seed = 72L)
      expect_equal(demo$data$g_oracle,
        .2 + .8 * exp(-.5 * (demo$data$X1 - .5)^2), tolerance = 1e-15)
      expect_equal(demo$truth$g_fn(matrix(c(-.1, .2, .5))),
        .2 + .8 * exp(-.5 * (c(-.1, .2, .5) - .5)^2), tolerance = 1e-15)
      expect_equal(demo$g_x0, .2 + .8 * exp(-.5 * (.2 - .5)^2),
        tolerance = 1e-15)
      if (seeded) expect_identical(.Random.seed, before) else
        expect_false(exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
    }
  })
})
