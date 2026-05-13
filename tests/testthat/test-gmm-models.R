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

test_that("new_gmm_model validates function slots and name vectors", {
  base_args <- list(
    name = "ok",
    setup = function(sample, n, ...) list(),
    moment_obs = function(sample, theta, ctx) matrix(0, nrow(sample), 1),
    coef_names = "theta"
  )

  expect_error(do.call(new_gmm_model, modifyList(base_args, list(setup = NULL))), "setup")
  expect_error(do.call(new_gmm_model, modifyList(base_args, list(moment_obs = NULL))), "moment_obs")
  expect_error(
    do.call(new_gmm_model, modifyList(base_args, list(jacobian_obs = 1))),
    "jacobian_obs"
  )
  expect_error(do.call(new_gmm_model, modifyList(base_args, list(start = 1))), "start")
  expect_error(
    do.call(new_gmm_model, modifyList(base_args, list(coef_names = character()))),
    "coef_names"
  )
  expect_error(
    do.call(new_gmm_model, modifyList(base_args, list(moment_names = c("ok", NA)))),
    "moment_names"
  )
})

test_that("gmm_model_dr_iv moment, Jacobian, and start match manual formulas", {
  dgp <- sim_iv_data("dgp2_iv_homog_overid", n = 70L, kx = 2L, q = 2L, seed = 21L)
  model <- gmm_model_dr_iv(
    z_cols = c("Z1", "Z2"),
    pi_cols = c("piZ1_oracle", "piZ2_oracle"),
    muY_col = "muY_oracle",
    muD_col = "muD_oracle"
  )
  theta <- 0.7
  ctx <- model$setup(sample = dgp$data, n = nrow(dgp$data))

  z_resid <- as.matrix(dgp$data[, c("Z1", "Z2")]) -
    as.matrix(dgp$data[, c("piZ1_oracle", "piZ2_oracle")])
  resid <- dgp$data$Y - dgp$data$muY_oracle -
    (dgp$data$D - dgp$data$muD_oracle) * theta
  manual_g <- z_resid * resid
  manual_jac <- z_resid * (-(dgp$data$D - dgp$data$muD_oracle))

  expect_equal(model$moment_obs(dgp$data, theta, ctx), manual_g)
  expect_equal(model$jacobian_obs(dgp$data, theta, ctx), manual_jac)

  w <- rep(1 / nrow(dgp$data), nrow(dgp$data))
  start <- model$start(sample = dgp$data, ctx = ctx, w = w, x0 = dgp$x0)
  expect_length(start, 1L)
  expect_true(is.finite(start))
})

test_that("gmm_model_dr_iv supports numeric column references", {
  dgp <- sim_iv_data("dgp2_iv_homog_overid", n = 60L, kx = 2L, q = 2L, seed = 22L)
  nm <- names(dgp$data)
  model_names <- gmm_model_dr_iv(
    z_cols = c("Z1", "Z2"),
    pi_cols = c("piZ1_oracle", "piZ2_oracle"),
    muY_col = "muY_oracle",
    muD_col = "muD_oracle"
  )
  model_idx <- gmm_model_dr_iv(
    y_col = match("Y", nm),
    d_col = match("D", nm),
    z_cols = match(c("Z1", "Z2"), nm),
    pi_cols = match(c("piZ1_oracle", "piZ2_oracle"), nm),
    muY_col = match("muY_oracle", nm),
    muD_col = match("muD_oracle", nm)
  )

  ctx_names <- model_names$setup(sample = dgp$data, n = nrow(dgp$data))
  ctx_idx <- model_idx$setup(sample = dgp$data, n = nrow(dgp$data))

  expect_equal(ctx_idx$y, ctx_names$y)
  expect_equal(ctx_idx$d, ctx_names$d)
  expect_equal(ctx_idx$z_resid, ctx_names$z_resid)
  expect_equal(
    model_idx$moment_obs(dgp$data, 0.5, ctx_idx),
    model_names$moment_obs(dgp$data, 0.5, ctx_names)
  )
})

test_that("gmm_model_dr_iv validates setup data explicitly", {
  dgp <- sim_iv_data("dgp2_iv_homog_overid", n = 50L, kx = 2L, q = 2L, seed = 23L)
  model <- gmm_model_dr_iv(
    z_cols = c("Z1", "Z2"),
    pi_cols = c("piZ1_oracle", "piZ2_oracle"),
    muY_col = "muY_oracle",
    muD_col = "muD_oracle"
  )

  expect_error(gmm_model_dr_iv(z_cols = "Z1"), "must all be provided")
  expect_error(model$setup(sample = as.list(dgp$data), n = nrow(dgp$data)), "data.frame")
  expect_error(model$setup(sample = dgp$data[-1L, ], n = nrow(dgp$data)), "nrow")

  bad_len_model <- gmm_model_dr_iv(
    z_cols = c("Z1", "Z2"),
    pi_cols = "piZ1_oracle",
    muY_col = "muY_oracle",
    muD_col = "muD_oracle"
  )
  expect_error(bad_len_model$setup(sample = dgp$data, n = nrow(dgp$data)),
               "same number")

  bad_na <- dgp$data
  bad_na$Y[1L] <- NA_real_
  expect_error(model$setup(sample = bad_na, n = nrow(bad_na)), "Y.*NA")

  bad_inf <- dgp$data
  bad_inf$piZ1_oracle[1L] <- Inf
  expect_error(model$setup(sample = bad_inf, n = nrow(bad_inf)), "pi.*non-finite")
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
