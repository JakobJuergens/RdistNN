.plot_recovery_fits <- function() {
  X <- matrix(seq_len(24L), ncol = 1L)
  y <- as.numeric(seq_len(nrow(X)))
  list(
    np = np_dnn(X = X, x0 = 0, Y = y, s = 4L, se = FALSE),
    cate = cate_dnn(X = X, x0 = 0, pseudo_outcome = y, s = 4L, se = FALSE),
    ratio = affine_ratio_dnn(X = X, x0 = 0, numerator = y,
      denominator = rep(1, length(y)), s = 4L, se = FALSE),
    plqr = plqr_dnn(X = X, x0 = 0, Y = y, D = matrix(1, length(y), 1L),
      g = rep(0, length(y)), mD = matrix(0, length(y), 1L),
      theta_eval = 0, s = 4L, se = FALSE)
  )
}

test_that("all tuning plot adapters preserve unrelated warnings and caller controls", {
  skip_if_not_installed("ggplot2")
  fits <- .plot_recovery_fits()
  tuning <- new_tuning("plot_recovery", eval = function(...) list(s = 4L),
    eval_grid = function(...) data.frame(lambda = c(.5, 1, 2), s = c(3L, 4L, 5L)))
  for (kind in names(fits)) {
    name <- switch(kind, np = ".np_dnn_tuning_grid", cate = ".cate_dnn_tuning_grid",
      ratio = ".affine_ratio_tuning_grid", plqr = ".plqr_dnn_tuning_grid")
    wrapper <- local({
      original <- get(name, envir = asNamespace("RdistNN"))
      function(...) {
        warning("unrelated tuning warning")
        original(...)
      }
    })
    switch(kind,
      np = local_mocked_bindings(.np_dnn_tuning_grid = wrapper, .package = "RdistNN"),
      cate = local_mocked_bindings(.cate_dnn_tuning_grid = wrapper, .package = "RdistNN"),
      ratio = local_mocked_bindings(.affine_ratio_tuning_grid = wrapper, .package = "RdistNN"),
      plqr = local_mocked_bindings(.plqr_dnn_tuning_grid = wrapper, .package = "RdistNN"))
    expect_warning(plot_tuning(fits[[kind]], tuning = tuning, se = FALSE),
      "unrelated tuning warning")
    withr::with_options(list(warn = 2L), {
      expect_error(plot_tuning(fits[[kind]], tuning = tuning, se = FALSE),
        "unrelated tuning warning")
      seen <- 0L
      plot <- withCallingHandlers(plot_tuning(fits[[kind]], tuning = tuning, se = FALSE),
        warning = function(w) {
          if (identical(conditionMessage(w), "unrelated tuning warning")) {
            seen <<- seen + 1L
            invokeRestart("muffleWarning")
          }
        })
      expect_identical(seen, 1L)
      expect_true(inherits(plot, "ggplot"))
    })
  }
})

test_that("captions reuse retained localization facts without diagnostic evaluation", {
  skip_if_not_installed("ggplot2")
  fits <- .plot_recovery_fits()
  local_mocked_bindings(plot_data_diagnostics = function(...) {
    stop("caption must not evaluate diagnostics")
  }, .package = "RdistNN")
  for (fit in fits) {
    concentration <- fit$meta$target_notices$weight_concentration
    expected <- sprintf("median n_eff = %.1f (%.1f-%.1f)",
      stats::median(concentration), min(concentration), max(concentration))
    ordinary <- plot_estimates(fit)
    expect_match(ordinary$labels$caption, expected, fixed = TRUE)
    fit$W <- NULL
    missing_extra <- plot_estimates(fit)
    expect_identical(missing_extra$labels$caption, ordinary$labels$caption)
    fit$meta$target_notices$weight_concentration[1L] <- NA_real_
    expect_error(plot_estimates(fit), class = "rdistnn_state_error")
  }
})
