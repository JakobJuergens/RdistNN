skip_if_not(
  all(vapply(
    c("plqr_dnn", "plqr_demo_data"),
    exists, logical(1),
    mode = "function"
  )),
  "plqr_dnn() or required helpers not found; skipping file"
)

.make_plqr_inference_data <- function(n = 90L, p = 2L, p_d = 2L, m = 2L, tau = 0.5) {
  demo <- plqr_demo_data(n = n, p = p, p_d = p_d, tau = tau, seed = 818)
  df <- demo$data
  X <- as.matrix(df[paste0("X", seq_len(p))])
  D <- as.matrix(df[paste0("D", seq_len(p_d))])
  mD <- as.matrix(df[paste0("mD", seq_len(p_d))])
  q <- as.numeric(df$q_oracle)
  Y <- as.numeric(df$Y)
  X0 <- matrix(rnorm(m * p), m, p)
  density <- 0.75 + 0.1 * seq_len(n) / n
  df$q_density <- density
  list(
    df = df, X = X, D = D, mD = mD, q = q, Y = Y, X0 = X0,
    density = density, beta = demo$truth$beta_true, tau = tau,
    x_cols = paste0("X", seq_len(p)),
    D_cols = paste0("D", seq_len(p_d)),
    mD_cols = paste0("mD", seq_len(p_d))
  )
}

.manual_plqr_sandwich <- function(d, w) {
  d_resid <- d$D - d$mD
  S <- d$tau * (1 - d$tau) * crossprod(d_resid, d_resid * w)
  J <- crossprod(d_resid, d$D * as.numeric(d$density * w))
  J_inv <- solve(J)
  omega <- J_inv %*% S %*% t(J_inv)
  vcov <- omega / (1 / sum(w^2))
  list(S = S, J = J, omega = omega, vcov = vcov)
}

test_that("plqr_dnn computes direct PLQR sandwich from supplied density", {
  d <- .make_plqr_inference_data()
  fit <- plqr_dnn(
    X = d$X, D = d$D, q = d$q, mD = d$mD, Y = d$Y,
    beta = d$beta, tau = d$tau, density = d$density,
    X0 = d$X0, method = "exact", k = 2L, s = 10L,
    se = TRUE, return = "W"
  )

  expect_s3_class(fit, "plqr_result")
  expect_equal(
    unname(fit$theta),
    matrix(d$beta, nrow = length(d$beta), ncol = nrow(d$X0))
  )
  expect_equal(length(fit$score_covariance), nrow(d$X0))
  expect_equal(length(fit$jacobian), nrow(d$X0))
  expect_equal(length(fit$omega), nrow(d$X0))
  expect_equal(length(fit$vcov), nrow(d$X0))
  expect_equal(dim(fit$se), c(length(d$beta), nrow(d$X0)))

  manual <- .manual_plqr_sandwich(d, fit$W[, 1L])
  expect_equal(unname(fit$score_covariance[[1L]]), unname(manual$S), tolerance = 1e-12)
  expect_equal(unname(fit$jacobian[[1L]]), unname(manual$J), tolerance = 1e-12)
  expect_equal(unname(fit$omega[[1L]]), unname(manual$omega), tolerance = 1e-10)
  expect_equal(unname(fit$vcov[[1L]]), unname(manual$vcov), tolerance = 1e-10)
  expect_equal(unname(fit$se[, 1L]), unname(sqrt(diag(manual$vcov))), tolerance = 1e-10)
  expect_equal(unname(fit$neff), unname(1 / colSums(fit$W^2)), tolerance = 1e-12)
  expect_identical(fit$meta$se_route, "direct_sandwich")

  txt <- capture.output(print(fit))
  expect_true(any(grepl("coefficient SEs available", txt, fixed = TRUE)))
  expect_true(any(grepl("coefficient SEs reported above", txt, fixed = TRUE)))
})

test_that("plqr_dnn direct sandwich agrees across matrix and data-frame interfaces", {
  d <- .make_plqr_inference_data(n = 80L, m = 3L)

  fit_mat <- plqr_dnn(
    X = d$X, D = d$D, q = d$q, mD = d$mD, Y = d$Y,
    beta = d$beta, tau = d$tau, density = d$density,
    X0 = d$X0, method = "exact", k = 2L, s = 9L,
    se = TRUE
  )
  fit_df <- plqr_dnn(
    data = d$df,
    x_cols = d$x_cols,
    D_cols = d$D_cols,
    q_col = "q_oracle",
    mD_cols = d$mD_cols,
    y_col = "Y",
    density_col = "q_density",
    beta = d$beta,
    tau = d$tau,
    X0 = d$X0,
    method = "exact", k = 2L, s = 9L,
    se = TRUE
  )

  expect_equal(fit_df$moment, fit_mat$moment, tolerance = 1e-12)
  expect_equal(fit_df$theta, fit_mat$theta, tolerance = 1e-12)
  expect_equal(fit_df$se, fit_mat$se, tolerance = 1e-12)
  expect_equal(fit_df$neff, fit_mat$neff, tolerance = 1e-12)
  expect_equal(fit_df$jacobian, fit_mat$jacobian, tolerance = 1e-12)
  expect_equal(fit_df$score_covariance, fit_mat$score_covariance, tolerance = 1e-12)
  expect_equal(fit_df$vcov, fit_mat$vcov, tolerance = 1e-12)
})

test_that("plqr_result methods report coefficient inference when direct sandwich is available", {
  d <- .make_plqr_inference_data(n = 85L, m = 2L)
  fit <- plqr_dnn(
    X = d$X, D = d$D, q = d$q, mD = d$mD, Y = d$Y,
    beta = d$beta, tau = d$tau, density = d$density,
    X0 = d$X0, method = "exact", k = 2L, s = 9L,
    se = TRUE
  )

  expect_equal(coef(fit), fit$theta)
  expect_equal(vcov(fit, target = 1L), fit$vcov[[1L]])

  s <- summary(fit)
  expect_s3_class(s, "summary.plqr_result")
  expect_true(is.data.frame(s$coefficients))
  expect_equal(
    colnames(s$coefficients),
    c("target", "term", "estimate", "std.error", "statistic", "p.value", "conf.low", "conf.high")
  )
  expect_equal(nrow(s$coefficients), length(d$beta) * nrow(d$X0))
  expect_true(all(is.finite(s$coefficients$std.error)))

  ci <- confint(fit, parm = 1L, target = 2L)
  expect_equal(nrow(ci), 1L)
  expect_named(ci, c("target", "term", "conf.low", "conf.high"))
  expect_error(vcov(fit, target = 3L), "target")
  expect_error(confint(fit, parm = "not_a_term"), "valid coefficient names", ignore.case = TRUE)
})

test_that("plqr_result methods fail clearly when direct sandwich is unavailable", {
  d <- .make_plqr_inference_data(n = 70L, m = 1L)
  fit <- plqr_dnn(
    X = d$X, D = d$D, q = d$q, mD = d$mD, Y = d$Y,
    beta = d$beta, tau = d$tau,
    X0 = d$X0, method = "exact", k = 2L, s = 8L,
    se = FALSE
  )

  expect_null(fit$theta_se)
  expect_error(vcov(fit), "not available", ignore.case = TRUE)
  expect_error(confint(fit), "not available", ignore.case = TRUE)

  s <- summary(fit)
  expect_false(isTRUE(s$has_inference))
  expect_true(all(is.na(s$coefficients$std.error)))
})

test_that("plqr_dnn validates direct sandwich density inputs", {
  d <- .make_plqr_inference_data(n = 50L, m = 1L)

  expect_error(
    plqr_dnn(
      X = d$X, D = d$D, q = d$q, mD = d$mD, Y = d$Y,
      beta = d$beta, tau = d$tau, density = d$density[-1L],
      X0 = d$X0, method = "exact", k = 1L, s = 6L, se = FALSE
    ),
    "density.*length",
    ignore.case = TRUE
  )

  density_zero <- d$density
  density_zero[1L] <- 0
  expect_error(
    plqr_dnn(
      X = d$X, D = d$D, q = d$q, mD = d$mD, Y = d$Y,
      beta = d$beta, tau = d$tau, density = density_zero,
      X0 = d$X0, method = "exact", k = 1L, s = 6L, se = FALSE
    ),
    "density.*positive",
    ignore.case = TRUE
  )

  expect_error(
    plqr_dnn(
      data = d$df,
      x_cols = d$x_cols,
      D_cols = d$D_cols,
      q_col = "q_oracle",
      mD_cols = d$mD_cols,
      y_col = "Y",
      density_col = "missing_density",
      beta = d$beta,
      tau = d$tau,
      X0 = d$X0,
      method = "exact", k = 1L, s = 6L, se = FALSE
    ),
    "density_col",
    ignore.case = TRUE
  )

  expect_error(
    plqr_dnn(
      X = d$X, D = d$D, q = d$q, mD = d$mD, Y = d$Y,
      beta = d$beta, tau = d$tau, density = d$density,
      X0 = d$X0, method = "exact", k = 1L,
      tuning_mode = "grid",
      se = TRUE
    ),
    "tuning_mode.*density",
    ignore.case = TRUE
  )
})

test_that("predict.plqr_result returns coefficient surfaces at new targets", {
  d <- .make_plqr_inference_data(n = 88L, m = 2L)
  new_X0 <- matrix(rnorm(3L * ncol(d$X)), 3L, ncol(d$X))
  fit <- plqr_dnn(
    X = d$X, D = d$D, q = d$q, mD = d$mD, Y = d$Y,
    beta = d$beta, tau = d$tau, density = d$density,
    X0 = d$X0, method = "exact", k = 2L, s = 9L,
    se = TRUE
  )
  refit <- plqr_dnn(
    X = d$X, D = d$D, q = d$q, mD = d$mD, Y = d$Y,
    beta = d$beta, tau = d$tau, density = d$density,
    X0 = new_X0, method = "exact", k = 2L, s = 9L,
    se = TRUE
  )

  pred <- predict(fit, X0 = new_X0)
  expect_equal(pred, refit$theta, tolerance = 1e-12)

  pred_se <- predict(fit, X0 = new_X0, se.fit = TRUE)
  expect_true(is.list(pred_se))
  expect_named(pred_se, c("fit", "se.fit"))
  expect_equal(pred_se$fit, refit$theta, tolerance = 1e-12)
  expect_equal(pred_se$se.fit, refit$se, tolerance = 1e-12)

  pred_ci <- predict(fit, X0 = new_X0, interval = "confidence")
  expect_true(is.list(pred_ci))
  expect_named(pred_ci, c("fit", "lwr", "upr"))
  z <- stats::qnorm(0.975)
  expect_equal(pred_ci$fit, refit$theta, tolerance = 1e-12)
  expect_equal(pred_ci$lwr, refit$theta - z * refit$se, tolerance = 1e-12)
  expect_equal(pred_ci$upr, refit$theta + z * refit$se, tolerance = 1e-12)

  pred_ci_se <- predict(fit, X0 = new_X0, interval = "confidence", se.fit = TRUE)
  expect_named(pred_ci_se, c("fit", "lwr", "upr", "se.fit"))
  expect_equal(pred_ci_se$se.fit, refit$se, tolerance = 1e-12)
})

test_that("predict.plqr_result accepts single x0 and data-frame newdata", {
  d <- .make_plqr_inference_data(n = 82L, m = 2L)
  fit <- plqr_dnn(
    data = d$df,
    x_cols = d$x_cols,
    D_cols = d$D_cols,
    q_col = "q_oracle",
    mD_cols = d$mD_cols,
    y_col = "Y",
    density_col = "q_density",
    beta = d$beta,
    tau = d$tau,
    X0 = d$X0,
    method = "exact", k = 2L, s = 9L,
    se = TRUE
  )

  pred_x0 <- predict(fit, x0 = d$X0[1L, ])
  expect_equal(dim(pred_x0), c(length(d$beta), 1L))
  expect_equal(unname(pred_x0[, 1L]), unname(d$beta), tolerance = 1e-12)

  nd <- data.frame(
    X1 = d$X0[, 1L],
    X2 = d$X0[, 2L],
    extra = seq_len(nrow(d$X0))
  )
  pred_df <- predict(fit, newdata = nd, se.fit = TRUE)
  pred_X0 <- predict(fit, X0 = as.matrix(nd[, d$x_cols, drop = FALSE]), se.fit = TRUE)

  expect_equal(pred_df$fit, pred_X0$fit, tolerance = 1e-12)
  expect_equal(pred_df$se.fit, pred_X0$se.fit, tolerance = 1e-12)
})

test_that("predict.plqr_result separates coefficient prediction from inference availability", {
  d <- .make_plqr_inference_data(n = 75L, m = 2L)
  new_X0 <- matrix(rnorm(2L * ncol(d$X)), 2L, ncol(d$X))
  fit <- plqr_dnn(
    X = d$X, D = d$D, q = d$q, mD = d$mD, Y = d$Y,
    beta = d$beta, tau = d$tau,
    X0 = d$X0, method = "exact", k = 2L, s = 8L,
    se = FALSE
  )

  pred <- predict(fit, X0 = new_X0)
  expect_equal(
    unname(pred),
    matrix(d$beta, nrow = length(d$beta), ncol = nrow(new_X0)),
    tolerance = 1e-12
  )
  expect_error(
    predict(fit, X0 = new_X0, se.fit = TRUE),
    "density",
    ignore.case = TRUE
  )
  expect_error(
    predict(fit, X0 = new_X0, interval = "confidence"),
    "density",
    ignore.case = TRUE
  )
})

test_that("predict.plqr_result validates inputs and fitted model state", {
  d <- .make_plqr_inference_data(n = 78L, m = 2L)
  fit <- plqr_dnn(
    X = d$X, D = d$D, q = d$q, mD = d$mD, Y = d$Y,
    beta = d$beta, tau = d$tau, density = d$density,
    X0 = d$X0, method = "exact", k = 2L, s = 8L,
    se = TRUE
  )

  expect_error(
    predict(fit),
    "newdata.*X0.*x0",
    ignore.case = TRUE
  )
  expect_error(
    predict(fit, X0 = d$X0, x0 = d$X0[1L, ]),
    "either `X0` or `x0`",
    fixed = TRUE
  )
  expect_error(
    predict(fit, newdata = d$X0, X0 = d$X0),
    "either `newdata` or `X0`/`x0`",
    fixed = TRUE
  )
  expect_error(
    predict(fit, X0 = d$X0, se.fit = "yes"),
    "se.fit.*TRUE or FALSE"
  )
  expect_error(
    predict(fit, X0 = d$X0, interval = "confidence", level = 1),
    "numeric scalar in \\(0,1\\)",
    ignore.case = TRUE
  )

  bad_state <- fit
  bad_state$.model <- NULL
  expect_error(
    predict(bad_state, X0 = d$X0),
    "fitted model state",
    ignore.case = TRUE
  )
})
