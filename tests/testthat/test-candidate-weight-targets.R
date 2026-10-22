test_that("candidate weight profiles honor target selection and reject unknown targets", {
  skip_if_not_installed("ggplot2")
  X <- matrix(c(0, 0, 0, 1, 2, 2, 3, 4), ncol = 1L)
  tuning <- new_tuning("candidate_target_profiles", eval = function(...) list(s = 3L),
    eval_grid = function(...) data.frame(s = c(2L, 3L, 5L)))
  make <- function(X) suppressWarnings(np_dnn(X = X, Y = c(1, 7, -3, 4, 10, -2, 5, 11),
    X0 = matrix(c(0, 2), ncol = 1L), tuning = tuning, tuning_mode = "grid",
    se = FALSE, return = "ranks"))
  for (training in list(X, matrix(seq_len(nrow(X)), ncol = 1L))) {
    fit <- make(training)
    expect_null(fit$W)
    plot <- suppressWarnings(plot_weights(fit, targets = 2L, candidates = "all"))
    expect_identical(unique(plot$data$target_id), 2L)
    expect_equal(length(unique(plot$data$candidate_id)), 3L)
    panels <- ggplot2::ggplot_build(plot)$layout$layout
    expect_identical(unique(panels$target_id), 2L)
    expect_error(suppressWarnings(plot_weights(fit, targets = 999L, candidates = "all")),
                 "targets.*1\\.\\.2")
  }
})
