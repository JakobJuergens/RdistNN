.test_dr_iv_model <- function(q) {
  gmm_model_dr_iv(
    z_cols = paste0("Z", seq_len(q)),
    pi_cols = paste0("piZ", seq_len(q), "_oracle"),
    muY_col = "muY_oracle",
    muD_col = "muD_oracle"
  )
}

.test_mean2_gmm_model <- function() {
  new_gmm_model(
    name = "mean2",
    setup = function(sample, n, ...) {
      y <- as.matrix(sample[, c("Y1", "Y2"), drop = FALSE])
      list(y = y)
    },
    moment_obs = function(sample, theta, ctx) {
      sweep(ctx$y, 2L, as.numeric(theta), FUN = "-")
    },
    jacobian_obs = function(sample, theta, ctx) {
      n <- nrow(ctx$y)
      out <- array(0, dim = c(n, 2L, 2L))
      out[, 1L, 1L] <- -1
      out[, 2L, 2L] <- -1
      out
    },
    start = function(sample, ctx, w, x0) {
      as.numeric(crossprod(ctx$y, w))
    },
    coef_names = c("theta1", "theta2"),
    moment_names = c("mean1", "mean2")
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

test_that("gmm_dnn supports numeric columns and tibble data like data frames", {
  dgp <- sim_iv_data("dgp2_iv_homog_overid", n = 110L, kx = 2L, q = 2L, seed = 24L)
  nm <- names(dgp$data)
  model_names <- .test_dr_iv_model(dgp$q)
  model_idx <- gmm_model_dr_iv(
    y_col = match("Y", nm),
    d_col = match("D", nm),
    z_cols = match(c("Z1", "Z2"), nm),
    pi_cols = match(c("piZ1_oracle", "piZ2_oracle"), nm),
    muY_col = match("muY_oracle", nm),
    muD_col = match("muD_oracle", nm)
  )

  fit_names <- gmm_dnn(
    data = dgp$data,
    x_cols = c("X1", "X2"),
    x0 = dgp$x0,
    model = model_names,
    k = 2L,
    s = 9L
  )
  fit_idx <- gmm_dnn(
    data = dgp$data,
    x_cols = match(c("X1", "X2"), nm),
    x0 = dgp$x0,
    model = model_idx,
    k = 2L,
    s = 9L
  )

  expect_equal(fit_idx$theta, fit_names$theta, tolerance = 1e-8)
  expect_equal(fit_idx$se, fit_names$se, tolerance = 1e-8)

  testthat::skip_if_not_installed("tibble")
  fit_tbl <- gmm_dnn(
    data = tibble::as_tibble(dgp$data),
    x_cols = c("X1", "X2"),
    x0 = dgp$x0,
    model = model_names,
    k = 2L,
    s = 9L
  )
  expect_equal(fit_tbl$theta, fit_names$theta, tolerance = 1e-8)
  expect_equal(fit_tbl$se, fit_names$se, tolerance = 1e-8)
})

test_that("gmm_dnn handles multi-parameter models and 3D analytic Jacobians", {
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

test_that("gmm_result methods validate availability and selectors", {
  dgp <- sim_iv_data("dgp2_iv_homog_overid", n = 100L, kx = 2L, q = 2L, seed = 27L)
  model <- .test_dr_iv_model(dgp$q)
  fit <- gmm_dnn(
    data = dgp$data,
    x_cols = c("X1", "X2"),
    X0 = rbind(dgp$x0, dgp$x0 + c(0.1, 0.1)),
    model = model,
    k = 2L,
    s = 9L
  )

  expect_equal(coef(fit), fit$theta)
  expect_equal(length(vcov(fit)), 2L)
  expect_equal(vcov(fit, target = 1L), fit$vcov[[1L]])
  expect_error(vcov(fit, target = 3L), "target")

  ci <- confint(fit, parm = "theta", target = 1L)
  expect_equal(nrow(ci), 1L)
  expect_named(ci, c("target", "term", "conf.low", "conf.high"))
  expect_error(confint(fit, parm = 2L), "parm")
  expect_error(summary(fit, level = 1), "level")

  fit_no_se <- gmm_dnn(
    data = dgp$data,
    x_cols = c("X1", "X2"),
    x0 = dgp$x0,
    model = model,
    k = 2L,
    s = 9L,
    se = FALSE
  )
  expect_error(vcov(fit_no_se), "not available")
  expect_error(confint(fit_no_se), "not available")
})

test_that("gmm_dnn validates user-facing interface controls", {
  dgp <- sim_iv_data("dgp2_iv_homog_overid", n = 70L, kx = 2L, q = 2L, seed = 28L)
  model <- .test_dr_iv_model(dgp$q)
  X <- as.matrix(dgp$data[, c("X1", "X2")])

  expect_error(gmm_dnn(data = dgp$data, x_cols = c("X1", "X2"), x0 = dgp$x0,
                       model = list()), "model")
  expect_error(gmm_dnn(data = dgp$data, X = X, x_cols = c("X1", "X2"),
                       x0 = dgp$x0, model = model), "data-frame interface")
  expect_error(gmm_dnn(X = X, x0 = dgp$x0, model = model), "Y")
  expect_error(gmm_dnn(X = X, model_data = dgp$data[-1L, ], x0 = dgp$x0,
                       model = model), "nrow")
  expect_error(gmm_dnn(data = dgp$data, x_cols = c("X1", "X2"), x0 = dgp$x0,
                       model = model, omega_ridge = -1), "omega_ridge")
  expect_error(gmm_dnn(data = dgp$data, x_cols = c("X1", "X2"), x0 = dgp$x0,
                       model = model, fd_step = 0), "fd_step")
  expect_error(gmm_dnn(data = dgp$data, x_cols = c("X1", "X2"), x0 = dgp$x0,
                       model = model, start = c(0, 1)), "start")
  expect_error(gmm_dnn(data = dgp$data, x_cols = c("X1", "X2"), x0 = dgp$x0,
                       X0 = matrix(dgp$x0, nrow = 1L), model = model),
               "either `X0`")
})

test_that("gmm_dnn detects model moment-name mismatches", {
  bad_model <- new_gmm_model(
    name = "bad_names",
    setup = function(sample, n, ...) list(y = sample$Y),
    moment_obs = function(sample, theta, ctx) cbind(ctx$y - theta[1L],
                                                    2 * (ctx$y - theta[1L])),
    start = function(sample, ctx, w, x0) sum(w * ctx$y),
    coef_names = "theta",
    moment_names = "only_one"
  )
  df <- data.frame(X1 = rnorm(50L), X2 = rnorm(50L), Y = rnorm(50L))

  expect_error(
    gmm_dnn(data = df, x_cols = c("X1", "X2"), x0 = c(0, 0),
            model = bad_model, k = 1L, s = 5L),
    "moment_names"
  )
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
