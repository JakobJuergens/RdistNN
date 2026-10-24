# GMM estimator identities, Jacobian equivalence, and return shapes.
test_that("gmm_dnn handles multi-parameter models and 3D analytic Jacobians", {
  withr::local_preserve_seed()
  set.seed(25L)
  n <- 90L
  X <- matrix(rnorm(n * 2L), n, 2L)
  colnames(X) <- c("X1", "X2")
  df <- data.frame(
    X,
    Y1 = 1 + X[, 1L] + rnorm(n, sd = 0.1),
    Y2 = -0.5 + X[, 2L] + rnorm(n, sd = 0.1)
  )
  x0 <- c(0, 0)
  model <- .test_mean2_gmm_model()

  fit <- gmm_dnn(
    data = df,
    x_cols = c("X1", "X2"),
    x0 = x0,
    model = model,
    k = 2L,
    s = 10L,
    return = c("W", "g_obs", "omega", "jacobian", "ess")
  )
  expected <- as.numeric(crossprod(as.matrix(df[, c("Y1", "Y2")]), fit$W[, 1L]))

  expect_equal(dim(fit$theta), c(2L, 1L))
  expect_equal(unname(fit$theta[, 1L]), expected, tolerance = 1e-5)
  expect_equal(dim(fit$se), c(2L, 1L))
  expect_equal(dim(vcov(fit)), c(2L, 2L))
  expect_equal(unname(fit$jacobian[[1L]]), -diag(2), tolerance = 1e-10)
  expect_equal(dim(fit$g_obs[[1L]]), c(n, 2L))
  expect_true(is.na(fit$J[1L]))
  expect_true(is.na(fit$df[1L]))
})

test_that("gmm_dnn q=1 DR-IV matches weighted Wald ratio", {
  dgp <- sim_iv_data("dgp1_iv_justid_linear", n = 120L, kx = 2L, q = 1L, seed = 3L)
  model <- .test_dr_iv_model(dgp$q)
  x_cols <- c("X1", "X2")
  k <- 2L
  s <- 10L

  fit <- gmm_dnn(
    data = dgp$data,
    x_cols = x_cols,
    x0 = dgp$x0,
    model = model,
    k = k,
    s = s,
    return = "ess"
  )

  wt <- weights_dnn(
    X = as.matrix(dgp$data[, x_cols]),
    x0 = dgp$x0,
    k = k,
    s = s,
    return = "W"
  )
  w <- as.numeric(wt$W[, 1L])
  z_resid <- dgp$data$Z1 - dgp$data$piZ1_oracle
  a_bar <- sum(w * z_resid * (dgp$data$Y - dgp$data$muY_oracle))
  b_bar <- sum(w * z_resid * (dgp$data$D - dgp$data$muD_oracle))
  theta_wald <- a_bar / b_bar

  expect_equal(fit$theta[1L, 1L], theta_wald, tolerance = 1e-6)
  expect_true(is.na(fit$J[1L]))
  expect_true(is.na(fit$df[1L]))
})

test_that("gmm_dnn return extras have stable shapes", {
  dgp <- sim_iv_data("dgp2_iv_homog_overid", n = 95L, kx = 2L, q = 2L, seed = 26L)
  model <- .test_dr_iv_model(dgp$q)
  X0 <- rbind(dgp$x0, dgp$x0 + c(0.2, -0.1))

  fit <- gmm_dnn(
    data = dgp$data,
    x_cols = c("X1", "X2"),
    X0 = X0,
    model = model,
    k = 2L,
    s = 9L,
    return = c("W", "ranks", "ord", "dx", "g_obs", "omega", "jacobian",
               "ess", "neighbors")
  )

  expect_equal(dim(fit$W), c(95L, 2L))
  expect_equal(dim(fit$ranks), c(95L, 2L))
  expect_equal(dim(fit$ord), c(95L, 2L))
  expect_equal(dim(fit$dx), c(95L, 2L))
  expect_equal(dim(fit$neighbors), c(2L, 2L))
  expect_equal(length(fit$omega), 2L)
  expect_equal(length(fit$jacobian), 2L)
  expect_equal(length(fit$g_obs), 2L)
  expect_equal(unname(fit$ess), unname(fit$neff), tolerance = 1e-12)
})

test_that("gmm_dnn q>1 DR-IV matches oracle efficient benchmark", {
  dgp <- sim_iv_data("dgp2_iv_homog_overid", n = 150L, kx = 2L, q = 2L, seed = 4L)
  model <- .test_dr_iv_model(dgp$q)
  x_cols <- c("X1", "X2")
  k <- 2L
  s <- 12L
  alpha <- 0.05

  fit <- gmm_dnn(
    data = dgp$data,
    x_cols = x_cols,
    x0 = dgp$x0,
    model = model,
    k = k,
    s = s,
    return = "ess"
  )
  oracle <- getFromNamespace(".method_iv_oracle", "RdistNN")(
    dgp = dgp, x_cols = x_cols, k = k, s = s, alpha = alpha
  )

  expect_equal(fit$theta[1L, 1L], oracle[["theta_hat"]], tolerance = 1e-5)
  expect_equal(fit$se[1L, 1L], oracle[["se"]], tolerance = 1e-5)
  expect_equal(unname(fit$J[1L]), oracle[["J"]], tolerance = 1e-5)
  expect_equal(unname(fit$df[1L]), oracle[["df"]])
})

test_that("gmm_dnn finite-difference and analytic Jacobians agree", {
  dgp <- sim_iv_data("dgp2_iv_homog_overid", n = 120L, kx = 2L, q = 2L, seed = 5L)
  model <- .test_dr_iv_model(dgp$q)
  x_cols <- c("X1", "X2")

  fit_analytic <- gmm_dnn(
    data = dgp$data,
    x_cols = x_cols,
    x0 = dgp$x0,
    model = model,
    k = 2L,
    s = 10L,
    jacobian = "analytic"
  )
  fit_fd <- gmm_dnn(
    data = dgp$data,
    x_cols = x_cols,
    x0 = dgp$x0,
    model = model,
    k = 2L,
    s = 10L,
    jacobian = "fd"
  )

  expect_equal(fit_analytic$theta, fit_fd$theta, tolerance = 1e-8)
  expect_equal(fit_analytic$se, fit_fd$se, tolerance = 1e-5)
})

test_that("gmm_dnn multi-target output shape and return warning are stable", {
  dgp <- sim_iv_data("dgp2_iv_homog_overid", n = 100L, kx = 2L, q = 2L, seed = 6L)
  model <- .test_dr_iv_model(dgp$q)
  X0 <- rbind(dgp$x0, dgp$x0 + c(0.1, -0.1))

  fit <- NULL
  expect_warning(
    fit <- gmm_dnn(
      data = dgp$data,
      x_cols = c("X1", "X2"),
      X0 = X0,
      model = model,
      k = 2L,
      s = 9L,
      se = FALSE,
      return = c("ess", "bad_return")
    ),
    "Unknown `return`"
  )

  expect_equal(dim(fit$theta), c(1L, 2L))
  expect_null(fit$se)
  expect_equal(length(fit$ess), 2L)
  expect_equal(nrow(fit$convergence), 2L)
})
test_that("L-BFGS-B uses its supported objective tolerance control", {
  X <- matrix(1:8)
  y <- c(1, 4, 2, 7, 3, 8, 6, 5)
  model <- new_gmm_model("lbfgsb_tolerance", setup = function(sample, n, ...) list(y = y),
    moment_obs = function(sample, theta, ctx) matrix(ctx$y - theta[1L]),
    jacobian_obs = function(sample, theta, ctx) matrix(-1, length(ctx$y)),
    coef_names = "mean")
  expect_no_warning(fit <- gmm_dnn(X = X, X0 = matrix(c(.1, 4.1)), model = model,
    s = 3, optimizer = "L-BFGS-B"))
  expect_equal(as.numeric(fit$theta), as.numeric(crossprod(y, fit$W)), tolerance = 1e-7)
  expect_identical(fit$convergence$target, c("target_1", "target_2"))
  expect_length(fit$convergence$convergence_first, 2L)
  expect_length(fit$convergence$convergence_second, 2L)
  expect_identical(fit$meta$optimizer, "L-BFGS-B")
})
