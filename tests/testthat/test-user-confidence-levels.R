test_that("valid extreme confidence levels retain finite normal bounds across result families", {
  X <- matrix(1:8)
  y <- c(1, 4, 2, 7, 3, 8, 6, 5)
  X0 <- matrix(c(.1, 4.1))
  model <- new_gmm_model("confidence_level_mean", setup = function(sample, n, ...) list(y = y),
    moment_obs = function(sample, theta, ctx) matrix(ctx$y - theta[1L]),
    jacobian_obs = function(sample, theta, ctx) matrix(-1, length(ctx$y)),
    coef_names = "mean")
  fits <- list(
    np = np_dnn(X = X, Y = y, X0 = X0, s = 3),
    tdnn = np_tdnn(X = X, Y = y, X0 = X0, s1 = 2, s2 = 4),
    cate = cate_dnn(X = X, pseudo_outcome = y, X0 = X0, s = 3),
    ratio = affine_ratio_dnn(X = X, numerator = y, denominator = rep(1, 8), X0 = X0, s = 3),
    plqr = plqr_dnn(X = X, Y = y, D = matrix(1, 8), mD = matrix(0, 8),
      g = rep(0, 8), theta_eval = 0, density = rep(1, 8), X0 = X0, s = 3),
    gmm = gmm_dnn(X = X, X0 = X0, model = model, s = 3)
  )
  for (level in c(1 - .Machine$double.eps / 2, 1e-12, .95)) {
    for (fit in fits) {
      ci <- confint(fit, level = level)
      if (is.data.frame(ci)) ci <- ci[, c("conf.low", "conf.high")]
      expect_true(all(is.finite(as.matrix(ci))))
      tab <- summary(fit, level = level)$coefficients
      expect_true(all(is.finite(tab$conf.low) & is.finite(tab$conf.high)))
      quantity <- if (inherits(fit, "rdistnn_plqr_result")) list(quantity = "theta") else list()
      frame <- do.call(plot_data_estimates, c(list(fit, level = level), quantity))
      expect_equal(frame$ci_lower, as.numeric(ci[, 1L]), tolerance = 1e-13)
      expect_equal(frame$ci_upper, as.numeric(ci[, 2L]), tolerance = 1e-13)
      if (!inherits(fit, "rdistnn_gmm_result")) {
        pred <- predict(fit, X0 = X0, interval = "confidence", level = level)
        if (inherits(fit, "rdistnn_plqr_result")) {
          expect_true(all(is.finite(pred$lwr) & is.finite(pred$upr)))
        } else {
          expect_true(all(is.finite(pred[, "lwr"]) & is.finite(pred[, "upr"])))
        }
      }
    }
    grid <- np_dnn(X = X, Y = y, X0 = X0,
      tuning = new_tuning("confidence_grid", eval = function(...) list(s = 3L),
        eval_grid = function(...) data.frame(s = 2:4)), tuning_mode = "grid")
    frame <- plot_data_estimates(grid, candidates = "all", level = level)
    expect_true(all(is.finite(frame$ci_lower) & is.finite(frame$ci_upper)))
    if (requireNamespace("ggplot2", quietly = TRUE)) {
      calibration <- plot_calibration(fits$plqr, level = level)$data
      expect_true(all(is.finite(calibration$band_lower) & is.finite(calibration$band_upper)))
    }
  }
})

test_that("normal critical values preserve tiny central probabilities", {
  for (level in c(.Machine$double.xmin * .Machine$double.eps, 1e-300,
                 1e-20, 1e-12, 1e-5, .5, 1 - .Machine$double.eps / 2)) {
    critical <- RdistNN:::.dnn_normal_critical(level)
    expect_true(is.finite(critical) && critical > 0)
    if (level < 1e-4) {
      central <- stats::integrate(stats::dnorm, -critical, critical, rel.tol = 1e-12,
        abs.tol = 0)$value
      expect_equal(central / level, 1, tolerance = 1e-14)
    } else {
      expect_identical(critical, stats::qnorm((1 - level) / 2, lower.tail = FALSE))
    }
  }
})

test_that("simulation metric intervals retain represented extreme tail probabilities", {
  for (alpha in c(.Machine$double.xmin * .Machine$double.eps, 1e-300, 1e-20, .05,
                  1 - .Machine$double.eps / 2)) {
    cate <- RdistNN:::.sim_metric(0, 1, 0, alpha)
    iv <- RdistNN:::.sim_metric_iv(0, 1, 0, alpha)
    expect_true(is.finite(cate$ci_len) && cate$ci_len > 0)
    expect_identical(cate$ci_len, unname(iv["ci_len"]))
    if (alpha < .9) {
      expect_equal(cate$ci_len / 2,
        stats::qnorm(log(alpha) - log(2), lower.tail = FALSE, log.p = TRUE), tolerance = 0)
    }
  }
})
