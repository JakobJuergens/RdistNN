.test_dr_iv_model <- function(q) {
  gmm_model_dr_iv(
    z_cols = paste0("Z", seq_len(q)),
    pi_cols = paste0("piZ", seq_len(q), "_oracle"),
    muY_col = "muY_oracle",
    muD_col = "muD_oracle"
  )
}

test_that("gmm_dnn DR-IV data-frame and matrix interfaces agree", {
  dgp <- sim_iv_data("dgp2_iv_homog_overid", n = 120L, kx = 2L, q = 2L, seed = 2L)
  model <- .test_dr_iv_model(dgp$q)
  x_cols <- c("X1", "X2")

  fit_df <- gmm_dnn(
    data = dgp$data,
    x_cols = x_cols,
    x0 = dgp$x0,
    model = model,
    k = 2L,
    s = 10L,
    return = c("ess", "omega", "jacobian")
  )
  fit_mat <- gmm_dnn(
    X = as.matrix(dgp$data[, x_cols]),
    model_data = dgp$data,
    x0 = dgp$x0,
    model = model,
    k = 2L,
    s = 10L,
    return = "ess"
  )

  expect_s3_class(fit_df, "gmm_result")
  expect_equal(fit_df$theta, fit_mat$theta, tolerance = 1e-8)
  expect_equal(fit_df$se, fit_mat$se, tolerance = 1e-8)
  expect_true(is.finite(fit_df$neff[1L]))
  expect_true(is.matrix(fit_df$omega[[1L]]))
  expect_true(is.matrix(fit_df$jacobian[[1L]]))
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
