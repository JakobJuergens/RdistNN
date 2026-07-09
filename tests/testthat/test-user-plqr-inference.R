skip_if_not(
  all(vapply(
    c("plqr_dnn", "plqr_demo_data"),
    exists, logical(1),
    mode = "function"
  )),
  "plqr_dnn() or required helpers not found; skipping file"
)

test_that("plqr_dnn computes direct PLQR sandwich from supplied density", {
  d <- .make_plqr_inference_data()
  fit <- plqr_dnn(
    X = d$X, D = d$D, g = d$g, mD = d$mD, Y = d$Y,
    theta_eval = d$theta_eval, tau = d$tau, density = d$density,
    X0 = d$X0, method = "exact", k = 2L, s = 10L,
    se = TRUE, return = "W"
  )

  expect_s3_class(fit, "plqr_result")
  expect_equal(
    unname(fit$theta),
    matrix(d$theta_eval, nrow = length(d$theta_eval), ncol = nrow(d$X0))
  )
  expect_equal(length(fit$score_covariance), nrow(d$X0))
  expect_equal(length(fit$jacobian), nrow(d$X0))
  expect_equal(length(fit$theta_omega), nrow(d$X0))
  expect_equal(length(fit$theta_vcov), nrow(d$X0))
  expect_equal(dim(fit$theta_se), c(length(d$theta_eval), nrow(d$X0)))

  manual <- .manual_plqr_sandwich(d, fit$W[, 1L])
  expect_equal(unname(fit$score_covariance[[1L]]), unname(manual$S), tolerance = 1e-12)
  expect_equal(unname(fit$jacobian[[1L]]), unname(manual$J), tolerance = 1e-12)
  expect_equal(unname(fit$theta_omega[[1L]]), unname(manual$omega), tolerance = 1e-10)
  expect_equal(unname(fit$theta_vcov[[1L]]), unname(manual$vcov), tolerance = 1e-10)
  expect_equal(unname(fit$theta_se[, 1L]), unname(sqrt(diag(manual$vcov))), tolerance = 1e-10)
  expect_equal(unname(fit$neff), unname(1 / colSums(fit$W^2)), tolerance = 1e-12)
  expect_identical(fit$meta$theta_se_method, "direct_sandwich")
  expect_identical(fit$meta$moment_se_method, "none")
  expect_identical(fit$meta$theta_source, "supplied")

  txt <- capture.output(print(fit))
  expect_true(any(grepl("coefficient SEs available", txt, fixed = TRUE)))
  expect_true(any(grepl("coefficient SEs reported above", txt, fixed = TRUE)))
})

test_that("plqr_dnn direct sandwich agrees across matrix and data-frame interfaces", {
  d <- .make_plqr_inference_data(n = 80L, m = 3L)

  fit_mat <- plqr_dnn(
    X = d$X, D = d$D, g = d$g, mD = d$mD, Y = d$Y,
    theta_eval = d$theta_eval, tau = d$tau, density = d$density,
    X0 = d$X0, method = "exact", k = 2L, s = 9L,
    se = TRUE
  )
  fit_df <- plqr_dnn(
    data = d$df,
    x_cols = d$x_cols,
    D_cols = d$D_cols,
    g_col = "g_oracle",
    mD_cols = d$mD_cols,
    y_col = "Y",
    density_col = "density_oracle",
    theta_eval = d$theta_eval,
    tau = d$tau,
    X0 = d$X0,
    method = "exact", k = 2L, s = 9L,
    se = TRUE
  )

  expect_equal(fit_df$moment, fit_mat$moment, tolerance = 1e-12)
  expect_equal(fit_df$theta, fit_mat$theta, tolerance = 1e-12)
  expect_equal(fit_df$theta_se, fit_mat$theta_se, tolerance = 1e-12)
  expect_equal(fit_df$neff, fit_mat$neff, tolerance = 1e-12)
  expect_equal(fit_df$jacobian, fit_mat$jacobian, tolerance = 1e-12)
  expect_equal(fit_df$score_covariance, fit_mat$score_covariance, tolerance = 1e-12)
  expect_equal(fit_df$theta_vcov, fit_mat$theta_vcov, tolerance = 1e-12)
})

test_that("plqr_result methods report coefficient inference when direct sandwich is available", {
  d <- .make_plqr_inference_data(n = 85L, m = 2L)
  fit <- plqr_dnn(
    X = d$X, D = d$D, g = d$g, mD = d$mD, Y = d$Y,
    theta_eval = d$theta_eval, tau = d$tau, density = d$density,
    X0 = d$X0, method = "exact", k = 2L, s = 9L,
    se = TRUE
  )

  expect_equal(coef(fit), fit$theta)
  expect_equal(vcov(fit, target = 1L), fit$theta_vcov[[1L]])

  s <- summary(fit)
  expect_s3_class(s, "summary.plqr_result")
  expect_true(is.data.frame(s$coefficients))
  expect_equal(
    colnames(s$coefficients),
    c("target", "term", "estimate", "std.error", "statistic", "p.value", "conf.low", "conf.high")
  )
  expect_equal(nrow(s$coefficients), length(d$theta_eval) * nrow(d$X0))
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
    X = d$X, D = d$D, g = d$g, mD = d$mD, Y = d$Y,
    theta_eval = d$theta_eval, tau = d$tau,
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
      X = d$X, D = d$D, g = d$g, mD = d$mD, Y = d$Y,
      theta_eval = d$theta_eval, tau = d$tau, density = d$density[-1L],
      X0 = d$X0, method = "exact", k = 1L, s = 6L, se = FALSE
    ),
    "density.*length",
    ignore.case = TRUE
  )

  density_zero <- d$density
  density_zero[1L] <- 0
  expect_error(
    plqr_dnn(
      X = d$X, D = d$D, g = d$g, mD = d$mD, Y = d$Y,
      theta_eval = d$theta_eval, tau = d$tau, density = density_zero,
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
      g_col = "g_oracle",
      mD_cols = d$mD_cols,
      y_col = "Y",
      density_col = "missing_density",
      theta_eval = d$theta_eval,
      tau = d$tau,
      X0 = d$X0,
      method = "exact", k = 1L, s = 6L, se = FALSE
    ),
    "density_col",
    ignore.case = TRUE
  )

  expect_error(
    plqr_dnn(
      X = d$X, D = d$D, g = d$g, mD = d$mD, Y = d$Y,
      theta_eval = d$theta_eval, tau = d$tau, density = d$density,
      X0 = d$X0, method = "exact", k = 1L,
      tuning_mode = "grid",
      se = TRUE
    ),
    "tuning_mode.*density",
    ignore.case = TRUE
  )

  expect_error(
    plqr_dnn(
      X = d$X, D = d$D, g = d$g, mD = d$mD, Y = d$Y,
      theta_eval = d$theta_eval, tau = d$tau, density = d$density,
      X0 = d$X0, method = "asymp", k = 1L, s = 6L,
      normalize = FALSE, se = TRUE
    ),
    "normalize.*direct PLQR coefficient inference",
    ignore.case = TRUE
  )
})

test_that("predict.plqr_result returns coefficient surfaces at new targets", {
  d <- .make_plqr_inference_data(n = 88L, m = 2L)
  new_X0 <- matrix(rnorm(3L * ncol(d$X)), 3L, ncol(d$X))
  fit <- plqr_dnn(
    X = d$X, D = d$D, g = d$g, mD = d$mD, Y = d$Y,
    theta_eval = d$theta_eval, tau = d$tau, density = d$density,
    X0 = d$X0, method = "exact", k = 2L, s = 9L,
    se = TRUE
  )
  refit <- plqr_dnn(
    X = d$X, D = d$D, g = d$g, mD = d$mD, Y = d$Y,
    theta_eval = d$theta_eval, tau = d$tau, density = d$density,
    X0 = new_X0, method = "exact", k = 2L, s = 9L,
    se = TRUE
  )

  pred <- predict(fit, X0 = new_X0)
  expect_equal(pred, refit$theta, tolerance = 1e-12)

  pred_se <- predict(fit, X0 = new_X0, se.fit = TRUE)
  expect_true(is.list(pred_se))
  expect_named(pred_se, c("fit", "se.fit"))
  expect_equal(pred_se$fit, refit$theta, tolerance = 1e-12)
  expect_equal(pred_se$se.fit, refit$theta_se, tolerance = 1e-12)

  pred_ci <- predict(fit, X0 = new_X0, interval = "confidence")
  expect_true(is.list(pred_ci))
  expect_named(pred_ci, c("fit", "lwr", "upr"))
  z <- stats::qnorm(0.975)
  expect_equal(pred_ci$fit, refit$theta, tolerance = 1e-12)
  expect_equal(pred_ci$lwr, refit$theta - z * refit$theta_se, tolerance = 1e-12)
  expect_equal(pred_ci$upr, refit$theta + z * refit$theta_se, tolerance = 1e-12)

  pred_ci_se <- predict(fit, X0 = new_X0, interval = "confidence", se.fit = TRUE)
  expect_named(pred_ci_se, c("fit", "lwr", "upr", "se.fit"))
  expect_equal(pred_ci_se$se.fit, refit$theta_se, tolerance = 1e-12)
})

test_that("predict.plqr_result accepts single x0 and data-frame newdata", {
  d <- .make_plqr_inference_data(n = 82L, m = 2L)
  fit <- plqr_dnn(
    data = d$df,
    x_cols = d$x_cols,
    D_cols = d$D_cols,
    g_col = "g_oracle",
    mD_cols = d$mD_cols,
    y_col = "Y",
    density_col = "density_oracle",
    theta_eval = d$theta_eval,
    tau = d$tau,
    X0 = d$X0,
    method = "exact", k = 2L, s = 9L,
    se = TRUE
  )

  pred_x0 <- predict(fit, x0 = d$X0[1L, ])
  expect_equal(dim(pred_x0), c(length(d$theta_eval), 1L))
  expect_equal(unname(pred_x0[, 1L]), unname(d$theta_eval), tolerance = 1e-12)

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
    X = d$X, D = d$D, g = d$g, mD = d$mD, Y = d$Y,
    theta_eval = d$theta_eval, tau = d$tau,
    X0 = d$X0, method = "exact", k = 2L, s = 8L,
    se = FALSE
  )

  pred <- predict(fit, X0 = new_X0)
  expect_equal(
    unname(pred),
    matrix(d$theta_eval, nrow = length(d$theta_eval), ncol = nrow(new_X0)),
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
    X = d$X, D = d$D, g = d$g, mD = d$mD, Y = d$Y,
    theta_eval = d$theta_eval, tau = d$tau, density = d$density,
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

test_that("predict.plqr_result validates fitted model state before computing", {
  d <- .make_plqr_inference_data(n = 60L, m = 2L)
  fit <- plqr_dnn(
    X = d$X, D = d$D, g = d$g, mD = d$mD, Y = d$Y,
    theta_eval = d$theta_eval, tau = d$tau, X0 = d$X0,
    method = "exact", k = 2L, s = 6L, se = FALSE
  )

  bad_x <- fit
  bad_x$.model$x_mat[1L, 1L] <- Inf
  expect_error(
    predict(bad_x, X0 = d$X0),
    "object$.model$x_mat cannot contain non-finite values",
    fixed = TRUE
  )

  bad_d <- fit
  bad_d$.model$d_mat <- bad_d$.model$d_mat[-1L, , drop = FALSE]
  expect_error(
    predict(bad_d, X0 = d$X0),
    "object$.model$d_mat/object$.model$mD_mat",
    fixed = TRUE
  )

  bad_mD_cols <- fit
  bad_mD_cols$.model$mD_mat <- cbind(bad_mD_cols$.model$mD_mat, 1)
  expect_error(
    predict(bad_mD_cols, X0 = d$X0),
    "object$.model$mD_mat",
    fixed = TRUE
  )

  bad_theta <- fit
  bad_theta$.model$theta_eval <- bad_theta$.model$theta_eval[-1L]
  expect_error(
    predict(bad_theta, X0 = d$X0),
    "object$.model$theta_eval",
    fixed = TRUE
  )

  bad_k <- fit
  bad_k$.model$k <- bad_k$.model$s + 1L
  expect_error(
    predict(bad_k, X0 = d$X0),
    "object$.model$k must satisfy k <= s",
    fixed = TRUE
  )

  bad_s <- fit
  bad_s$.model$s <- nrow(bad_s$.model$x_mat) + 1L
  expect_error(
    predict(bad_s, X0 = d$X0),
    "object$.model$s must satisfy s <= n",
    fixed = TRUE
  )

  bad_method <- fit
  bad_method$.model$method <- "bogus"
  expect_error(
    predict(bad_method, X0 = d$X0),
    "object$.model$method must be one of",
    fixed = TRUE
  )

  bad_normalize <- fit
  bad_normalize$.model$normalize <- NA
  expect_error(
    predict(bad_normalize, X0 = d$X0),
    "object$.model$normalize must be TRUE or FALSE",
    fixed = TRUE
  )
})

test_that("diagnostics report target-level Jacobian conditioning for direct-sandwich fits", {
  d <- .make_plqr_inference_data(n = 70L, m = 2L)
  fit <- plqr_dnn(
    X = d$X, D = d$D, g = d$g, mD = d$mD, Y = d$Y,
    theta_eval = d$theta_eval, tau = d$tau, density = d$density,
    X0 = d$X0, method = "exact", k = 2L, s = 8L,
    se = TRUE, return = "W"
  )
  dd <- diagnostics(fit)

  cols <- c(
    "jacobian_available", "jacobian_min_singular_value",
    "jacobian_condition_number", "jacobian_rank",
    "density_min_weighted", "density_median_weighted",
    "density_max_weighted", "theta_vcov_available"
  )
  expect_true(all(cols %in% names(dd)))
  expect_true(all(dd$jacobian_available))
  expect_true(all(is.finite(dd$jacobian_min_singular_value)))
  expect_true(all(dd$jacobian_min_singular_value > 0))
  expect_true(all(is.finite(dd$jacobian_condition_number)))
  expect_true(all(dd$jacobian_rank == length(d$theta_eval)))
  expect_true(all(dd$theta_vcov_available))
  expect_true(all(dd$density_min_weighted >= min(d$density)))
  expect_true(all(dd$density_max_weighted <= max(d$density)))
  expect_true(all(is.finite(dd$density_median_weighted)))
})

test_that("diagnostics flag singular PLQR Jacobians and unavailable coefficient inference", {
  d <- .make_plqr_inference_data(n = 60L, m = 1L)
  # D == mD makes the residualized design identically zero, so J is singular.
  d_singular <- d
  d_singular$mD <- d_singular$D
  fit <- suppressWarnings(plqr_dnn(
    X = d_singular$X, D = d_singular$D, g = d_singular$g,
    mD = d_singular$mD, Y = d_singular$Y,
    theta_eval = d_singular$theta_eval, tau = d_singular$tau,
    density = d_singular$density,
    X0 = d_singular$X0, method = "exact", k = 2L, s = 8L,
    se = TRUE, return = "W"
  ))
  dd <- diagnostics(fit)

  expect_true(all(dd$jacobian_available))
  expect_true(all(dd$jacobian_condition_number == Inf))
  expect_true(all(dd$jacobian_rank == 0L))
  expect_false(any(dd$theta_vcov_available))

  # Fits without density report the Jacobian block as unavailable.
  fit_no_density <- plqr_dnn(
    X = d$X, D = d$D, g = d$g, mD = d$mD, Y = d$Y,
    theta_eval = d$theta_eval, tau = d$tau,
    X0 = d$X0, method = "exact", k = 2L, s = 8L,
    se = FALSE, return = "W"
  )
  dd_no <- diagnostics(fit_no_density)
  expect_false(any(dd_no$jacobian_available))
  expect_true(all(is.na(dd_no$jacobian_min_singular_value)))
  expect_false(any(dd_no$theta_vcov_available))
})

test_that("direct sandwich internals reject unnormalized weights", {
  d <- .make_plqr_inference_data(n = 40L, m = 1L)
  fit <- plqr_dnn(
    X = d$X, D = d$D, g = d$g, mD = d$mD, Y = d$Y,
    theta_eval = d$theta_eval, tau = d$tau, density = d$density,
    X0 = d$X0, method = "exact", k = 2L, s = 6L,
    se = TRUE, return = "W"
  )

  # Rescaled weight columns must fail clearly rather than silently producing
  # scale-dependent coefficient SEs.
  expect_error(
    RdistNN:::.plqr_direct_sandwich(
      d_mat = d$D, mD_mat = d$mD, density_vec = d$density,
      W = 2 * fit$W, tau = d$tau,
      coef_names = colnames(d$D), target_names = "target_1"
    ),
    "normalized weight columns",
    ignore.case = TRUE
  )

  # The jackknife-meat route shares the same internal guard.
  W_corrupt <- fit$W
  W_corrupt[1L, 1L] <- W_corrupt[1L, 1L] + 0.5
  expect_error(
    RdistNN:::.plqr_jackknife_meat_sandwich(
      x_mat = d$X, x0_mat = d$X0,
      score_mat = (d$tau - as.numeric(d$Y <= as.vector(d$D %*% d$theta_eval) + d$g)) * (d$D - d$mD),
      d_mat = d$D, mD_mat = d$mD, density_vec = d$density,
      W = W_corrupt,
      distance = "l2", method = "exact", k = 2L, s = 6L,
      normalize = TRUE, ctx = NULL, jackknife = "fast",
      precomputed_ranks = NULL,
      coef_names = colnames(d$D), target_names = "target_1"
    ),
    "normalized weight columns",
    ignore.case = TRUE
  )
})
