test_that("new_gmm_model constructs and validates model objects", {
  mod <- new_gmm_model(
    name = "mean_zero",
    setup = function(sample, n, ...) list(g = sample$g),
    moment_obs = function(sample, theta, ctx) matrix(ctx$g - theta[1], ncol = 1),
    jacobian_obs = function(sample, theta, ctx) matrix(-1, length(ctx$g), 1),
    start = function(sample, ctx, w, x0) sum(w * ctx$g),
    coef_names = "theta",
    moment_names = "mean"
  )

  expect_true(is_gmm_model(mod))
  expect_equal(mod$p, 1L)
  expect_error(
    new_gmm_model(
      name = "",
      setup = function(sample, n, ...) list(),
      moment_obs = function(sample, theta, ctx) matrix(0, n, 1),
      coef_names = "theta"
    ),
    "name"
  )
})

test_that("gmm_model_dr_iv rejects missing model columns", {
  dgp <- sim_iv_data("dgp2_iv_homog_overid", n = 60L, kx = 2L, q = 2L, seed = 1L)
  model <- gmm_model_dr_iv(
    z_cols = c("Z1", "bad_Z"),
    pi_cols = c("piZ1_oracle", "piZ2_oracle"),
    muY_col = "muY_oracle",
    muD_col = "muD_oracle"
  )

  expect_error(
    gmm_dnn(
      data = dgp$data,
      x_cols = c("X1", "X2"),
      x0 = dgp$x0,
      model = model,
      k = 2L,
      s = 8L,
      se = FALSE
    ),
    "bad_Z"
  )
})
