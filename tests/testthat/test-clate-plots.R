.clate_plot_fit <- function(method = "exact", se = TRUE, perfect = FALSE) {
  a <- .clate_fit_arguments()
  a$X0 <- cbind(seq(.1, .9, length.out = 6L), 0)
  a$return <- c("W", "dx", "ranks")
  if (perfect) {
    a$D <- a$A
    a$m_d0 <- rep(0, 8L)
    a$m_d1 <- rep(1, 8L)
  }
  do.call(clate_dnn, c(a, list(method = method, se = se)))
}

test_that("CLATE estimate plot styles and generic dispatch mirror CATE", {
  skip_if_not_installed("ggplot2")
  for (method in c("exact", "asymp")) {
    fit <- .clate_plot_fit(method)
    profile <- plot_estimates(fit)
    expect_s3_class(profile, "ggplot")
    expect_true(all(c("GeomHline", "GeomRibbon", "GeomLine") %in% .plot_layer_geoms(profile)))
    expect_equal(profile$data$estimate, fit$theta)
    expect_equal(profile$data$ci_lower, as.numeric(confint(fit)[, 1L]))
    expect_equal(profile$data$ci_upper, as.numeric(confint(fit)[, 2L]))
    expect_s3_class(ggplot2::autoplot(fit), "ggplot")
    expect_true("GeomVline" %in% .plot_layer_geoms(plot_estimates(fit, style = "forest")))
    expect_true("GeomPointrange" %in% .plot_layer_geoms(plot_estimates(fit, style = "pointrange")))
    for (type in c("profile", "cumulative", "kernel")) {
      p <- plot_weights(fit, type = type)
      expect_s3_class(p, "ggplot")
      expect_s3_class(ggplot2::ggplot_build(p), "ggplot_built")
    }
  }
  point <- plot_estimates(.clate_plot_fit(se = FALSE))
  expect_false("GeomRibbon" %in% .plot_layer_geoms(point))
  expect_match(point$labels$caption, "no standard errors")
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)
  ret <- withVisible(plot(.clate_plot_fit(perfect = TRUE)))
  expect_false(ret$visible)
  expect_s3_class(ret$value, "ggplot")
})

test_that("CLATE overlap plots use the instrument and its retained propensity", {
  skip_if_not_installed("ggplot2")
  fit <- .clate_plot_fit()
  iv <- fit$.model$iv
  p <- plot_overlap(fit, type = "mass")
  expect_identical(unique(p$data$arm), c("instrument 1", "instrument 0"))
  expect_equal(p$data$mass[p$data$arm == "instrument 1"], as.numeric(crossprod(iv$A, fit$W)))
  expect_equal(p$data$mass[p$data$arm == "instrument 0"], as.numeric(crossprod(1 - iv$A, fit$W)))
  expect_gt(max(abs(as.numeric(crossprod(iv$A - iv$D, fit$W)))), .01)
  density <- plot_overlap(fit, type = "propensity", targets = 2L)
  expect_equal(density$data$e, iv$propensity$used)
  expect_equal(sum(density$data$w), 1)
  expect_match(density$labels$x, "Instrument")
  expect_error(plot_overlap(fit, treatment = iv$D), "treatment.*instrument")
  expect_error(plot_overlap(fit, share_limits = c(.9, .1)), "share_limits")
})

test_that("CLATE plot data survive serialization and candidate extraction", {
  skip_if_not_installed("ggplot2")
  a <- .clate_fit_arguments()
  a$s <- NULL
  a$tuning <- .clate_alignment_tuning()
  a$return <- c("W", "dx")
  fit <- do.call(clate_dnn, c(a, list(se = FALSE)))
  file <- tempfile(fileext = ".rds")
  on.exit(unlink(file))
  saveRDS(fit, file)
  saved <- readRDS(file)
  sensitivity <- plot_tuning(saved, se = TRUE)
  expect_true(all(is.finite(sensitivity$data$se)))
  expect_equal(sort(unique(sensitivity$data$s)), c(2L, 3L, 6L))
  expect_warning(grid <- do.call(clate_dnn, c(a, list(tuning_mode = "grid"))), "infeasible")
  view <- candidate_fit(grid, 1L)
  expect_equal(plot_estimates(view, style = "pointrange")$data$estimate, view$theta)
  for (p in list(plot_estimates(saved), plot_overlap(saved), plot_influence(saved),
                  sensitivity, plot_tuning(view))) {
    expect_s3_class(ggplot2::ggplot_build(p), "ggplot_built")
  }
  bare <- fit
  bare$W <- NULL
  expect_error(plot_overlap(bare), "return.*W")
  expect_error(plot_influence(bare), "return.*W")
  bad <- .clate_plot_fit()
  bad$theta[1L] <- NA_real_
  expect_error(plot_influence(bad, targets = 1L), "undefined|finite")
})

test_that("ratio plot bounds retain per-target Wald availability", {
  skip_if_not_installed("ggplot2")
  expect_warning(fit <- affine_ratio_dnn(X = matrix(1:8), X0 = matrix(c(0, 9)),
    numerator = rep(2, 8), denominator = rep(1, 8),
    tuning = .clate_alignment_tuning(), tuning_mode = "grid"), "infeasible")
  data <- plot_data_estimates(fit, candidates = "all")
  expect_true(all(is.na(data$ci_lower)))
  expect_true(all(is.na(data$ci_upper)))
  expect_true(all(data$estimate == 2))
  expect_s3_class(ggplot2::ggplot_build(plot_estimates(fit, candidates = "all")), "ggplot_built")
})

test_that("ratio plots explain unavailable uncertainty and draw finite points cleanly", {
  skip_if_not_installed("ggplot2")
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)
  zero <- .ratio_fixture(denominator = rep(0, 4L))
  p <- plot_estimates(zero, style = "pointrange")
  expect_match(p$labels$caption, "1 estimate.*unavailable")
  expect_false(grepl("se = FALSE", p$labels$caption, fixed = TRUE))
  expect_warning(ggplot2::ggplotGrob(p), NA)
  constant <- .ratio_fixture(numerator = rep(2, 4L), denominator = rep(1, 4L))
  p <- plot_estimates(constant, style = "forest")
  expect_match(p$labels$caption, "1 Wald interval.*unavailable")
  expect_warning(ggplot2::ggplotGrob(p), NA)
  expect_equal(p$data$estimate, 2)
  expect_true(is.na(p$data$ci_lower))
})

test_that("CLATE influence plots retain paired score contributions", {
  skip_if_not_installed("ggplot2")
  fit <- .clate_plot_fit()
  p <- plot_influence(fit, targets = 2L)
  signal <- fit$.model$signal
  reference <- fit$W[, 2L] * (signal[, 1L] - fit$theta[2L] * signal[, 2L])
  expect_equal(p$data$contribution, reference[p$data$obs_id], tolerance = 1e-12)
  expect_equal(sum(p$data$contribution), 0, tolerance = 1e-12)
  expect_match(p$labels$y, "CLATE")
  expect_error(plot_influence(fit, pseudo_outcome = rep(0, 8L)), "paired|pseudo_outcome")
  perfect <- .clate_plot_fit(perfect = TRUE)
  cate <- cate_dnn(X = perfect$.model$x_mat, X0 = perfect$x0,
    pseudo_outcome = perfect$.model$signal[, 1L], s = 3L, return = c("W", "dx"))
  expect_equal(plot_influence(perfect, targets = 2L)$data$contribution,
                 plot_influence(cate, targets = 2L)$data$contribution, tolerance = 1e-12)
})

test_that("CLATE grid plots and dashboards use stored candidate data", {
  skip_if_not_installed("ggplot2")
  a <- .clate_fit_arguments()
  a$s <- NULL
  a$tuning <- .clate_alignment_tuning()
  a$return <- c("W", "dx")
  expect_warning(fit <- do.call(clate_dnn, c(a, list(tuning_mode = "grid", method = "asymp"))), "infeasible")
  p <- plot_estimates(fit, candidates = "all", style = "pointrange")
  expect_identical(unique(p$data$coef), "theta")
  expect_true("GeomHline" %in% .plot_layer_geoms(p))
  stored <- fit$tuning_grid$estimates
  lookup <- match(paste(p$data$candidate_id, p$data$target_id), paste(stored$candidate_id, stored$target_id))
  expect_equal(p$data$estimate, stored$estimate[lookup])
  expect_equal(p$data$se, stored$se[lookup])
  expect_s3_class(plot_weights(fit, candidates = "all"), "ggplot")
  tuning <- plot_tuning(fit)
  expect_equal(sort(unique(tuning$data$s)), c(2L, 3L, 6L))
  expect_identical(unique(tuning$data$coef), "theta")
  expect_s3_class(ggplot2::autoplot(diagnostics(fit, candidates = "all")), "ggplot")
  testthat::local_mocked_bindings(.plot_pkg_installed = function(pkg) pkg != "patchwork", .package = "RdistNN")
  expect_message(panels <- plot_dashboard(fit), "patchwork")
  expect_true(all(c("estimates", "overlap", "propensity", "weights", "tuning") %in% names(panels)))
  expect_true(all(vapply(panels, inherits, logical(1L), what = "ggplot")))
})

test_that("CLATE plots preserve optional dependency and missing-weight guards", {
  fit <- .clate_plot_fit()
  testthat::local_mocked_bindings(.plot_pkg_installed = function(pkg) FALSE, .package = "RdistNN")
  expect_error(plot_estimates(fit), "ggplot2")
  expect_error(plot_overlap(fit), "ggplot2")
  expect_error(plot_influence(fit), "ggplot2")
  expect_error(plot_dashboard(fit), "ggplot2")
})
