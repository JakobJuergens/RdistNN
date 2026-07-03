# Weight-geometry panels (plan panels NP-2, NP-3): weight decay by rank,
# cumulative concentration, and the weight-vs-distance "effective kernel"
# view. The data layer is plot_data_weights().

#' Plot DNN weight profiles, concentration, and effective kernels
#'
#' Three views of the localization weights stored on a fitted object.
#' `type = "profile"` draws weight against rank with a guide at the exact
#' support edge `n - s + 1`; `type = "cumulative"` draws the cumulative
#' weight share with guides at 0.5 and 0.9 and markers at the top-1/5/10
#' shares; `type = "kernel"` draws weight against stored distance — the
#' effective data-driven kernel — with a guide at the rank-`s` radius.
#'
#' When more than eight targets are stored, panels default to a subset
#' (first, median effective sample size, last); override with `targets =`.
#' Two-scale (TDNN) fits have signed weights: a dashed zero line is added
#' automatically whenever negative weights are shown.
#'
#' @param object A fitted result object holding the weight matrix `W`
#'   (stored by default; `type = "kernel"` additionally requires stored
#'   distances via `return = c("W", "dx")`).
#' @param type One of `"profile"`, `"cumulative"`, or `"kernel"`.
#' @param targets Optional integer vector of target ids to display.
#' @param log_weights Logical (default `FALSE`); for rank-domain views,
#'   drop non-positive weights and use a log10 y scale (the number of
#'   dropped weights is reported in the caption).
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
#'     X = x_mat, x0 = c(0, 0), Y = y_vec,
#'     k = 1L, s = 12L, se = FALSE, return = c("W", "dx")
#'   )
#'   plot_weights(fit, type = "kernel")
#' }
#' @export
plot_weights <- function(object,
                         type = c("profile", "cumulative", "kernel"),
                         targets = NULL,
                         log_weights = FALSE,
                         ...) {
  .plot_require_ggplot2("plot_weights")
  type <- match.arg(type)

  dw <- plot_data_weights(object)
  n <- max(dw$rank)
  all_targets <- sort(unique(dw$target_id))
  m <- length(all_targets)

  if (type == "kernel" && all(is.na(dw$distance))) {
    .abort_bad_arg(
      "object",
      "must store distances for type = \"kernel\"; refit with `return = c(\"W\", \"dx\")`",
      quoted = FALSE
    )
  }

  neff_by_target <- vapply(
    all_targets,
    function(j) 1 / sum(dw$weight[dw$target_id == j]^2),
    numeric(1L)
  )
  targets <- .plot_select_targets(all_targets, neff_by_target, targets)
  dw <- dw[dw$target_id %in% targets, , drop = FALSE]
  dw$target <- factor(dw$target_id)
  multi <- length(targets) > 1L

  captions <- character(0)
  signed <- any(dw$weight < 0)
  if (signed) {
    captions <- c(captions, "two-scale weights are signed")
  }

  base_aes <- function(xcol, ycol) {
    if (multi) {
      ggplot2::aes(
        x = .data[[xcol]], y = .data[[ycol]], colour = .data$target
      )
    } else {
      ggplot2::aes(x = .data[[xcol]], y = .data[[ycol]])
    }
  }
  col_args <- if (multi) list() else list(colour = .plot_col_single)

  s_meta <- object$meta$s

  if (type == "profile") {
    if (isTRUE(log_weights)) {
      dropped <- sum(dw$weight <= 0)
      dw <- dw[dw$weight > 0, , drop = FALSE]
      captions <- c(
        captions,
        sprintf("log scale: %d non-positive weights dropped", dropped)
      )
    }
    p <- ggplot2::ggplot(dw, base_aes("rank", "weight")) +
      do.call(ggplot2::geom_step, col_args)
    if (signed && !isTRUE(log_weights)) {
      p <- p + ggplot2::geom_hline(yintercept = 0, linetype = "dashed")
    }
    if (!is.null(s_meta)) {
      p <- p + ggplot2::geom_vline(
        xintercept = n - as.integer(s_meta) + 1L,
        linetype = "dotted"
      )
    }
    if (isTRUE(log_weights)) {
      p <- p + ggplot2::scale_y_log10()
    }
    p <- p + ggplot2::labs(x = "Rank", y = "DNN weight")
  } else if (type == "cumulative") {
    p <- ggplot2::ggplot(dw, base_aes("rank", "cum_share")) +
      do.call(ggplot2::geom_line, col_args) +
      ggplot2::geom_hline(yintercept = c(0.5, 0.9), linetype = "dotted") +
      ggplot2::geom_point(
        data = dw[dw$rank %in% c(1L, 5L, 10L), , drop = FALSE]
      ) +
      ggplot2::labs(x = "Rank", y = "Cumulative weight share")
  } else {
    p <- ggplot2::ggplot(dw, base_aes("distance", "weight")) +
      do.call(ggplot2::geom_step, col_args)
    if (signed) {
      p <- p + ggplot2::geom_hline(yintercept = 0, linetype = "dashed")
    }
    if (!is.null(s_meta)) {
      radius <- vapply(
        targets,
        function(j) dw$distance[dw$target_id == j][as.integer(s_meta)],
        numeric(1L)
      )
      p <- p + ggplot2::geom_vline(
        xintercept = radius, linetype = "dotted"
      )
    }
    p <- p + ggplot2::labs(x = "Distance", y = "DNN weight")
  }

  if (multi) {
    p <- p + ggplot2::labs(colour = "target")
  }
  honesty <- .plot_honesty_caption(object)
  if (nzchar(honesty)) captions <- c(captions, honesty)
  if (length(captions) > 0L) {
    p <- p + ggplot2::labs(caption = paste(captions, collapse = "\n"))
  }
  p + .rdistnn_theme()
}
