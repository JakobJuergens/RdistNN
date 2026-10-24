test_that("predict.rdistnn_plqr_result returns coefficient surfaces at new targets", {
  .test_rng_scope({
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
    expect_equal(pred, refit$theta, tolerance = 1e-12, ignore_attr = "theta_source")

    pred_se <- predict(fit, X0 = new_X0, se.fit = TRUE)
    expect_true(is.list(pred_se))
    expect_named(pred_se, c("fit", "se.fit"))
    expect_equal(pred_se$fit, refit$theta, tolerance = 1e-12, ignore_attr = "theta_source")
    expect_equal(pred_se$se.fit, refit$theta_se, tolerance = 1e-12)

    pred_ci <- predict(fit, X0 = new_X0, interval = "confidence")
    expect_true(is.list(pred_ci))
    expect_named(pred_ci, c("fit", "lwr", "upr"))
    z <- stats::qnorm(0.975)
    expect_equal(pred_ci$fit, refit$theta, tolerance = 1e-12, ignore_attr = "theta_source")
    expect_equal(pred_ci$lwr, refit$theta - z * refit$theta_se, tolerance = 1e-12, ignore_attr = "theta_source")
    expect_equal(pred_ci$upr, refit$theta + z * refit$theta_se, tolerance = 1e-12, ignore_attr = "theta_source")

    pred_ci_se <- predict(fit, X0 = new_X0, interval = "confidence", se.fit = TRUE)
    expect_named(pred_ci_se, c("fit", "lwr", "upr", "se.fit"))
    expect_equal(pred_ci_se$se.fit, refit$theta_se, tolerance = 1e-12)
  }, seed = 818L)
})

test_that("predict.rdistnn_plqr_result accepts single x0 and data-frame newdata", {
  .test_rng_scope({
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
  }, seed = 818L)
})

test_that("predict.rdistnn_plqr_result separates coefficient prediction from inference availability", {
  .test_rng_scope({
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
      tolerance = 1e-12,
      ignore_attr = "theta_source"
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
  }, seed = 818L)
})

test_that("predict.rdistnn_plqr_result validates inputs and fitted model state", {
  .test_rng_scope({
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
    error <- expect_error(predict(bad_state, X0 = d$X0), class = "rdistnn_state_error")
    expect_identical(error$code, "state_invalid")
    expect_match(error$argument, "object$.model", fixed = TRUE)
    expect_match(conditionMessage(error), "[Rr]efit")
  }, seed = 818L)
})

test_that("predict.rdistnn_plqr_result validates fitted model state before computing", {
  .test_rng_scope({
    d <- .make_plqr_inference_data(n = 60L, m = 2L)
    fit <- plqr_dnn(
      X = d$X, D = d$D, g = d$g, mD = d$mD, Y = d$Y,
      theta_eval = d$theta_eval, tau = d$tau, X0 = d$X0,
      method = "exact", k = 2L, s = 6L, se = FALSE
    )

    bad_x <- fit
    bad_x$.model$x_mat[1L, 1L] <- Inf
    error <- expect_error(predict(bad_x, X0 = d$X0), class = "rdistnn_state_error")
    expect_identical(error$code, "state_invalid")
    expect_match(conditionMessage(error), "[Rr]efit")

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
    error <- expect_error(
      predict(bad_k, X0 = d$X0),
      "meta\\$k/\\.model\\$k.*inconsistent.*[Rr]efit",
      class = "rdistnn_state_error"
    )
    expect_identical(error$code, "state_invalid")

    bad_s <- fit
    bad_s$.model$s <- nrow(bad_s$.model$x_mat) + 1L
    expect_error(
      predict(bad_s, X0 = d$X0),
      "meta\\$s/\\.model\\$s.*inconsistent.*[Rr]efit",
      class = "rdistnn_state_error"
    )

    bad_method <- fit
    bad_method$.model$method <- "bogus"
    expect_error(
      predict(bad_method, X0 = d$X0),
      "method.*exact.*asymp.*[Rr]efit",
      class = "rdistnn_state_error"
    )

    bad_normalize <- fit
    bad_normalize$.model$normalize <- NA
    expect_error(
      predict(bad_normalize, X0 = d$X0),
      "normalize.*TRUE or FALSE.*[Rr]efit",
      class = "rdistnn_state_error"
    )
    bad_names <- fit
    bad_names$.model$coef_names <- NULL
    expect_error(predict(bad_names, X0 = d$X0), "object\\$\\.model\\$coef_names")
  }, seed = 818L)
})

test_that("PLQR numeric selectors retain named prediction columns after reordering", {
  d <- .make_plqr_inference_data(n = 40L, m = 2L, p_d = 1L)
  args <- list(data = d$df, x_cols = match(d$x_cols, names(d$df)),
    D_cols = d$D_cols, g_col = "g_oracle", mD_cols = d$mD_cols,
    y_col = "Y", theta_eval = d$theta_eval, X0 = d$X0, s = 5L)
  direct <- do.call(plqr_dnn, c(args, list(density_col = "density_oracle")))
  nd <- data.frame(extra = c(20, -20), X2 = d$X0[, 2L], X1 = d$X0[, 1L])
  expected <- predict(direct, X0 = d$X0, se.fit = TRUE)
  for (fit in list(direct, unserialize(serialize(direct, NULL)))) {
    expect_equal(predict(fit, newdata = nd, se.fit = TRUE), expected,
      tolerance = 1e-12)
    expect_equal(predict(fit, newdata = nd[c("X1", "extra", "X2")],
      se.fit = TRUE), expected, tolerance = 1e-12)
  }

  args$s <- NULL
  args$tuning <- new_tuning("plqr_prediction_selectors",
    eval = function(...) list(s = 5L),
    eval_grid = function(...) data.frame(s = c(3L, 5L, 7L)))
  args$tuning_mode <- "grid"
  grid <- do.call(plqr_dnn, args)
  for (fit in candidate_fits(unserialize(serialize(grid, NULL)))) {
    # Grid fits cannot carry density; this prediction route must still select
    # the reporting coordinates by their original names.
    expect_equal(predict(fit, newdata = nd), predict(fit, X0 = d$X0))
    expect_error(predict(fit, newdata = nd["extra"]), "Column.*X1.*not found")
  }
})

test_that("PLQR saved prediction freezes registered geometry and reuses setup", {
  d <- .make_plqr_inference_data(n = 36L, m = 2L, p_d = 1L)
  setup_calls <- 0L
  setup_sizes <- integer()
  metric <- new_distance("plqr_saved_geometry",
    setup = function(X) {
      setup_calls <<- setup_calls + 1L
      setup_sizes <<- c(setup_sizes, nrow(X))
      list(scale = c(1, 6))
    },
    eval = function(X, x0, ctx) {
      rowSums(sweep(sweep(X, 2L, x0, "-"), 2L, ctx$scale, "/")^2)
    })
  old <- get_distance("l2")
  on.exit(register_distance("l2", old), add = TRUE)
  register_distance("l2", metric)
  args <- list(X = d$X, D = d$D, g = d$g, mD = d$mD, Y = d$Y,
    density = d$density, theta_eval = d$theta_eval, X0 = d$X0,
    s = 5L, return = character())
  fit <- do.call(plqr_dnn, args)
  expect_identical(setup_calls, 1L)
  expected <- predict(fit, X0 = d$X0, se.fit = TRUE)
  expect_identical(setup_calls, 1L)
  expect_equal(expected$se.fit, fit$theta_se, tolerance = 1e-12)
  register_distance("l2", distance_l1())
  restored <- unserialize(serialize(fit, NULL))
  expect_equal(predict(restored, X0 = d$X0, se.fit = TRUE), expected,
    tolerance = 1e-12)

  # Moment-only grid candidates reuse the same fitted geometry, too.
  args$density <- args$s <- NULL
  args$distance <- metric
  args$tuning <- new_tuning("plqr_prediction_geometry",
    eval = function(...) list(s = 5L),
    eval_grid = function(...) data.frame(s = c(3L, 5L, 7L)))
  args$tuning_mode <- "grid"
  setup_calls <- 0L
  grid <- do.call(plqr_dnn, args)
  expect_identical(setup_calls, 1L)
  for (candidate in candidate_fits(grid)) predict(candidate, X0 = d$X0)
  expect_identical(setup_calls, 1L)

  # The requested LOO route keeps its deletion-specific setup behavior.
  args$tuning <- args$tuning_mode <- NULL
  args$s <- 5L
  args$jackknife <- "loo"
  setup_calls <- 0L
  setup_sizes <- integer()
  invisible(do.call(plqr_dnn, args))
  expect_gt(setup_calls, 1L)
  expect_identical(sum(setup_sizes == nrow(d$X)), 1L)
  expect_true(all(setup_sizes[-1L] == nrow(d$X) - 1L))
})
