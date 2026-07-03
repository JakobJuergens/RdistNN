# Diagnostics panel (plan panel NP-4): n_eff, maximum weight, and top-5
# weight share across targets, drawn from a dnn_diagnostics frame. The
# n_eff facet carries the degenerate-weights guide at 0.9 * n so the panel
# and the runtime warning tell one story.

#' Plot per-target weight diagnostics
#'
#' `autoplot()` on a [diagnostics()] frame draws three facets across
#' targets: effective sample size (`n_eff`, with a dotted guide at
#' `0.9 * n`, the degenerate-weights warning threshold), the maximum
#' weight, and the top-5 weight share. `plot()` prints the same panel.
#'
#' @param object A `dnn_diagnostics` object from [diagnostics()].
#' @param x A `dnn_diagnostics` object (for the `plot()` method).
#' @param ... Ignored; reserved for future options.
#'
#' @return A ggplot object.
#'
#' @examples
#' if (requireNamespace("ggplot2", quietly = TRUE)) {
#'   set.seed(1)
#'   n <- 80
#'   x_mat <- matrix(rnorm(n * 2), n, 2)
#'   y_vec <- x_mat[, 1] + rnorm(n, sd = 0.3)
#'   fit <- np_dnn(
#'     X = x_mat, X0 = cbind(seq(-1, 1, length.out = 5), 0), Y = y_vec,
#'     k = 1L, s = 12L, se = FALSE
#'   )
#'   ggplot2::autoplot(diagnostics(fit))
#' }
#' @name autoplot.dnn_diagnostics
NULL

#' @rdname autoplot.dnn_diagnostics
#' @exportS3Method ggplot2::autoplot dnn_diagnostics
autoplot.dnn_diagnostics <- function(object, ...) {
  .plot_require_ggplot2("autoplot.dnn_diagnostics")

  dd <- as.data.frame(object, stringsAsFactors = FALSE)
  class(dd) <- "data.frame"
  m <- nrow(dd)
  metrics <- c(
    n_eff = "n_eff",
    max_weight = "max weight",
    top5_weight_share = "top-5 weight share"
  )
  long <- do.call(rbind, lapply(names(metrics), function(col) {
    data.frame(
      target_id = seq_len(m),
      metric = metrics[[col]],
      value = as.numeric(dd[[col]]),
      stringsAsFactors = FALSE
    )
  }))
  long$metric <- factor(long$metric, levels = unname(metrics))

  guide <- data.frame(
    metric = factor(metrics[["n_eff"]], levels = unname(metrics)),
    yintercept = 0.9 * dd$n[1L]
  )

  p <- ggplot2::ggplot(
    long,
    ggplot2::aes(x = .data$target_id, y = .data$value)
  ) +
    ggplot2::geom_hline(
      data = guide,
      ggplot2::aes(yintercept = .data$yintercept),
      linetype = "dotted"
    ) +
    ggplot2::geom_line(colour = .plot_col_single) +
    ggplot2::geom_point(colour = .plot_col_single) +
    ggplot2::facet_wrap(ggplot2::vars(.data$metric), scales = "free_y") +
    ggplot2::labs(
      x = "target", y = NULL,
      caption = sprintf(
        "n = %d; dotted guide: n_eff = 0.9 n (degenerate-weights warning threshold)",
        as.integer(dd$n[1L])
      )
    )
  p + .rdistnn_theme()
}

#' @rdname autoplot.dnn_diagnostics
#' @export
plot.dnn_diagnostics <- function(x, ...) {
  p <- autoplot.dnn_diagnostics(x, ...)
  print(p)
  invisible(p)
}
