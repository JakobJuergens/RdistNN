test_that("public estimate and dashboard plots issue one duplicate-policy warning", {
  skip_if_not_installed("ggplot2")
  X <- matrix(c(0, 0, 0, 1, 2, 2, 3, 4))
  Y <- c(1, 7, -3, 4, 10, -2, 5, 11)
  targets <- matrix(c(0, 1, 2, 3))
  fits <- list(
    np = suppressWarnings(np_dnn(X = X, X0 = targets, Y = Y, s = 5,
                                return = c("W", "dx"))),
    cate = suppressWarnings(cate_dnn(X = X, X0 = targets,
      pseudo_outcome = Y, s = 5, return = c("W", "dx")))
  )
  capture <- function(expr) {
    warnings <- character()
    value <- withCallingHandlers(expr, warning = function(w) {
      warnings <<- c(warnings, conditionMessage(w))
      invokeRestart("muffleWarning")
    })
    list(value = value, count = sum(grepl(
      "Exact distance ties detected", warnings, fixed = TRUE)))
  }
  for (fit in fits) {
    estimate <- capture(plot_estimates(fit))
    expect_s3_class(estimate$value, "ggplot")
    expect_identical(estimate$count, 1L)
    dashboard <- capture(suppressMessages(plot_dashboard(fit)))
    if (requireNamespace("patchwork", quietly = TRUE)) {
      expect_s3_class(dashboard$value, "patchwork")
    } else {
      expect_true(is.list(dashboard$value))
      expect_true(all(vapply(dashboard$value, inherits, logical(1), "ggplot")))
    }
    expect_identical(dashboard$count, 1L)
  }
  # A finished public call must not suppress the next independent call.
  expect_identical(capture(plot_estimates(fits$np))$count, 1L)

  panels <- list(
    overlap = function() plot_overlap(fits$cate, treatment = rep(0:1, 4)),
    influence = function() plot_influence(fits$cate),
    autoplot = function() ggplot2::autoplot(fits$cate)
  )
  tdnn <- suppressWarnings(np_tdnn(X = X, X0 = targets, Y = Y,
                                  s1 = 3L, s2 = 5L, se = FALSE))
  plqr <- suppressWarnings(plqr_dnn(X = X, X0 = targets,
    D = matrix(seq_len(8)), Y = Y, g = rep(0, 8), mD = matrix(0, 8, 1),
    theta_eval = 0, tau = 0.5, s = 5L, se = FALSE))
  panels$tdnn <- function() plot_tdnn_scales(tdnn)
  panels$calibration <- function() plot_calibration(plqr)
  panels$conditioning <- function() plot_conditioning(plqr)
  for (panel in panels) {
    result <- capture(panel())
    expect_s3_class(result$value, "ggplot")
    expect_identical(result$count, 1L)
  }
})

test_that("repeated candidate orders retain distinct public plotting identities", {
  .with_distance_tie_warning({
  skip_if_not_installed("ggplot2")
  tuning <- new_tuning("plot_audit_repeated_orders",
    eval = function(...) list(s = 5L),
    eval_grid = function(...) data.frame(s = c(7L, 3L, 5L, 3L)))
  fit <- np_dnn(X = matrix(seq_len(12)), Y = sin(seq_len(12)),
    X0 = matrix(c(2, 6, 10)), tuning = tuning, tuning_mode = "grid",
    return = c("W", "ranks", "dx"))
  expected <- fit$tuning_grid$candidates
  frame <- plot_data_estimates(fit, candidates = "all")
  identities <- unique(frame[, c("candidate_id", "candidate_label", "s")])
  expect_setequal(identities$candidate_id, expected$candidate_id)
  expect_equal(nrow(identities), nrow(expected))
  expect_equal(length(unique(identities$candidate_label)), nrow(expected))
  for (builder in list(plot_estimates, plot_weights)) {
    plot <- builder(fit, candidates = "all")
    expect_s3_class(plot, "ggplot")
    expect_setequal(plot$data$candidate_id, expected$candidate_id)
    expect_equal(length(unique(plot$data$candidate_label)), nrow(expected))
    # Distinct rows with identical resolved orders still have identical values.
    repeated <- plot$data[plot$data$s == 3L, ]
    expect_equal(length(unique(repeated$candidate_id)), 2L)
  }
  expect_s3_class(plot_tuning(fit), "ggplot")
  })
})

test_that("plot target selectors reject values that would change during coercion", {
  .with_distance_tie_warning({
  skip_if_not_installed("ggplot2")
  fit <- np_dnn(X = matrix(seq_len(12)), Y = sin(seq_len(12)),
    X0 = matrix(c(2, 6, 10)), s = 4L, return = c("W", "dx"))
  for (invalid in list(1.5, "1", TRUE, 1 + 1i, NA_real_, Inf, numeric(), 0, 4)) {
    for (builder in list(plot_weights, plot_tuning)) {
      expect_error(builder(fit, targets = invalid), "targets")
    }
  }
  for (builder in list(plot_weights, plot_tuning)) {
    plot <- builder(fit, targets = c(3, 1, 3))
    expect_setequal(plot$data$target_id, c(3L, 1L))
  }
  })
})

test_that("plot scope checks preserve actionable class validation", {
  skip_if_not_installed("ggplot2")
  for (builder in list(plot_overlap, plot_influence, plot_calibration,
                       plot_conditioning, plot_tdnn_scales)) {
    expect_error(builder(1), "supports .*fits; got class")
  }
})
