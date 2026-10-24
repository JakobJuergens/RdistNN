# TDNN scale-decomposition panel (plan panel NP-5): the two-scale estimate
# overlaid on its single-scale DNN components. The single-scale curves are
# obtained by internal refits (plan decision 10.2); nothing new is stored
# on the fitted object.

#' Plot the TDNN two-scale decomposition
#'
#' Overlays the two-scale TDNN estimate (bold) with the single-scale DNN
#' estimates at `s1` and `s2` (muted) obtained by internal refits on the
#' stored training data. The comparison shows what the two-scale
#' bias-cancelling combination does relative to its components.
#'
#' @param object A `rdistnn_tdnn_np_result` from [np_tdnn()].
#' @param ... Ignored; reserved for future options.
#'
#' @return A ggplot object.
#'
#' @examples
#' if (requireNamespace("ggplot2", quietly = TRUE)) {
#'   set.seed(1)
#'   n <- 100
#'   x_mat <- matrix(rnorm(n * 2), n, 2,
#'                   dimnames = list(NULL, c("X1", "X2")))
#'   y_vec <- x_mat[, 1] + rnorm(n, sd = 0.3)
#'   fit <- np_tdnn(
#'     X = x_mat, X0 = cbind(seq(-1, 1, length.out = 6), 0), Y = y_vec,
#'     s1 = 5L, s2 = 10L, se = FALSE
#'   )
#'   plot_tdnn_scales(fit)
#' }
#' @export
plot_tdnn_scales <- function(object, ...) {
  .dnn_reuse_scope(object)
  .plot_require_ggplot2("plot_tdnn_scales")
  if (!inherits(object, "rdistnn_tdnn_np_result")) {
    .abort_bad_arg(
      "object",
      sprintf(
        "plot_tdnn_scales() supports rdistnn_tdnn_np_result fits; got class '%s'",
        class(object)[1L]
      ),
      quoted = FALSE
    )
  }
  model <- object$.model
  x_mat <- model$x_mat
  y_vec <- model$y_vec
  distance <- model$distance
  method <- model$method
  normalize <- model$normalize
  ctx <- model$ctx

  s1 <- as.integer(object$meta$s1)
  s2 <- as.integer(object$meta$s2)
  refit <- function(s) {
    np_dnn(
      X = x_mat, X0 = object$x0, Y = y_vec,
      distance = distance, method = method,
      k = 1L, s = s, normalize = normalize, ctx = ctx,
      se = FALSE, return = character(0)
    )
  }
  fit_s1 <- refit(s1)
  fit_s2 <- refit(s2)

  df_est <- function(fit, label) {
    d <- plot_data_estimates(fit)
    d$estimator <- label
    d
  }
  levels_est <- c(
    sprintf("DNN(s = %d)", s1),
    sprintf("DNN(s = %d)", s2),
    sprintf("TDNN(%d, %d)", s1, s2)
  )
  df <- rbind(
    df_est(fit_s1, levels_est[1L]),
    df_est(fit_s2, levels_est[2L]),
    df_est(object, levels_est[3L])
  )
  df$estimator <- factor(df$estimator, levels = levels_est)

  axis <- .plot_profile_axis(plot_data_estimates(object))
  xcol <- if (!is.null(axis$var)) axis$var else "target_id"

  p <- ggplot2::ggplot(
    df,
    ggplot2::aes(
      x = .data[[xcol]], y = .data$estimate,
      colour = .data$estimator, linewidth = .data$estimator
    )
  ) +
    ggplot2::geom_line() +
    ggplot2::geom_point(size = 1) +
    ggplot2::scale_colour_manual(
      values = stats::setNames(
        c("grey55", "grey35", .plot_col_treated), levels_est
      )
    ) +
    ggplot2::scale_linewidth_manual(
      values = stats::setNames(c(0.4, 0.4, 1), levels_est),
      guide = "none"
    ) +
    ggplot2::labs(
      x = axis$label, y = "Estimate", colour = NULL,
      caption = .plot_honesty_caption(object)
    )
  p + .rdistnn_theme()
}
