# Estimate panels (plan panels NP-1, C-1, C-2) and the autoplot()/plot()
# methods for np_result and cate_result. ggplot2 is guarded as a suggested
# dependency; the data layer is plot_data_estimates().

#' Plot estimates with confidence bounds across target points
#'
#' The estimate panel of the plotting layer: point estimates with normal
#' confidence bounds from the stored jackknife standard errors, drawn as a
#' profile band when the targets vary along a single stored coordinate, as
#' a forest layout for small non-profile target sets, and as pointranges
#' otherwise. `autoplot()` and `plot()` on `np_result` and `cate_result`
#' objects delegate here.
#'
#' Every panel carries an honesty caption with the median effective sample
#' size across targets (when weight diagnostics are available) and the
#' kernel order `s` with its provenance (`paper_default` or `manual`).
#' CATE panels draw a dashed zero reference line.
#'
#' @param object A fitted result object accepted by
#'   [plot_data_estimates()].
#' @param level Confidence level for the bounds (default `0.95`).
#' @param style One of `"auto"`, `"profile"`, `"pointrange"`, or
#'   `"forest"`. `"auto"` picks a profile band when a profile axis exists
#'   and there are at least four targets, a forest layout for `2..12`
#'   non-profile targets, and pointranges otherwise.
#' @param x_var Optional profile-axis override: a stored coordinate name
#'   (see `attr(plot_data_estimates(object), "x0_names")`), an `x0_*`
#'   column name, or a coordinate index.
#' @param shade_significant Logical (default `FALSE`); for profile-style
#'   CATE panels, add a rug along the axis where the confidence interval
#'   excludes zero. Intervals are marked, never point estimates.
#' @param ... Passed to methods (for `autoplot()`/`plot()`) or ignored.
#'
#' @return A ggplot object.
#'
#' @examples
#' if (requireNamespace("ggplot2", quietly = TRUE)) {
#'   set.seed(1)
#'   n <- 80
#'   x_mat <- matrix(rnorm(n * 2), n, 2,
#'                   dimnames = list(NULL, c("X1", "X2")))
#'   y_vec <- x_mat[, 1] + rnorm(n, sd = 0.3)
#'   X0 <- cbind(seq(-1, 1, length.out = 6), 0)
#'   fit <- np_dnn(X = x_mat, X0 = X0, Y = y_vec, k = 1L, s = 12L)
#'   plot_estimates(fit)
#' }
#' @export
plot_estimates <- function(object, level = 0.95,
                           style = c("auto", "profile", "pointrange", "forest"),
                           x_var = NULL,
                           shade_significant = FALSE,
                           ...) {
  .plot_require_ggplot2("plot_estimates")
  style <- match.arg(style)

  df <- plot_data_estimates(object, level = level)
  m <- length(unique(df$target_id))
  n_coef <- length(unique(df$coef))
  axis <- .plot_profile_axis(df, x_var = x_var)
  mode <- .plot_estimates_mode(style, has_profile = !is.null(axis$var), m = m)
  if (mode == "profile" && is.null(axis$var)) {
    .abort_bad_arg(
      "style",
      "= \"profile\" requires a profile axis: targets varying in exactly one stored coordinate, or an explicit `x_var`",
      quoted = FALSE
    )
  }

  has_se <- !all(is.na(df$se))
  captions <- character(0)
  if (!has_se) {
    captions <- c(captions, "no standard errors stored (se = FALSE)")
  }
  honesty <- .plot_honesty_caption(object)
  if (nzchar(honesty)) captions <- c(captions, honesty)

  # Zero reference lines are mandatory where zero is meaningful (plan
  # 3.0): treatment effects and PLQR moment coordinates.
  zero_line <- inherits(object, "cate_result") || inherits(object, "plqr_result")

  if (mode == "forest") {
    labs_map <- .plot_target_labels(df)
    df$target_label <- factor(
      labs_map[as.character(df$target_id)],
      levels = rev(unname(labs_map))
    )
    p <- ggplot2::ggplot(
      df,
      ggplot2::aes(x = .data$estimate, y = .data$target_label)
    )
    if (zero_line) {
      p <- p + ggplot2::geom_vline(xintercept = 0, linetype = "dashed")
    }
    if (has_se) {
      p <- p + ggplot2::geom_pointrange(
        ggplot2::aes(xmin = .data$ci_lower, xmax = .data$ci_upper),
        colour = .plot_col_single
      )
    } else {
      p <- p + ggplot2::geom_point(colour = .plot_col_single)
    }
    p <- p + ggplot2::labs(x = "Estimate", y = NULL)
  } else {
    xcol <- if (!is.null(axis$var)) axis$var else "target_id"
    p <- ggplot2::ggplot(df, ggplot2::aes(x = .data[[xcol]], y = .data$estimate))
    if (zero_line) {
      p <- p + ggplot2::geom_hline(yintercept = 0, linetype = "dashed")
    }
    if (mode == "profile") {
      if (has_se) {
        p <- p + ggplot2::geom_ribbon(
          ggplot2::aes(ymin = .data$ci_lower, ymax = .data$ci_upper),
          alpha = 0.25, fill = .plot_col_single
        )
      }
      p <- p + ggplot2::geom_line(colour = .plot_col_single)
      if (isTRUE(shade_significant) && has_se) {
        sig <- df[!is.na(df$ci_lower) &
                    sign(df$ci_lower) == sign(df$ci_upper) &
                    df$ci_lower != 0, , drop = FALSE]
        if (nrow(sig) > 0L) {
          p <- p + ggplot2::geom_rug(
            data = sig, sides = "b", colour = .plot_col_single,
            length = ggplot2::unit(0.03, "npc")
          )
        }
      }
    } else {
      if (has_se) {
        p <- p + ggplot2::geom_pointrange(
          ggplot2::aes(ymin = .data$ci_lower, ymax = .data$ci_upper),
          colour = .plot_col_single
        )
      } else {
        p <- p + ggplot2::geom_point(colour = .plot_col_single)
      }
    }
    p <- p + ggplot2::labs(x = axis$label, y = "Estimate")
  }

  if (n_coef > 1L) {
    p <- p + ggplot2::facet_wrap(ggplot2::vars(.data$coef), scales = "free_y")
  }
  if (length(captions) > 0L) {
    p <- p + ggplot2::labs(caption = paste(captions, collapse = "\n"))
  }
  p + .rdistnn_theme()
}

#' @rdname plot_estimates
#' @exportS3Method ggplot2::autoplot np_result
autoplot.np_result <- function(object, ...) {
  plot_estimates(object, ...)
}

#' @rdname plot_estimates
#' @exportS3Method ggplot2::autoplot cate_result
autoplot.cate_result <- function(object, ...) {
  plot_estimates(object, ...)
}

#' @rdname plot_estimates
#' @exportS3Method ggplot2::autoplot plqr_result
autoplot.plqr_result <- function(object, ...) {
  p <- plot_estimates(object, ...)
  beta_txt <- paste(format(as.numeric(object$beta), digits = 3),
                    collapse = ", ")
  p + ggplot2::labs(
    subtitle = sprintf("PLQR moments at beta = (%s), tau = %.2f",
                       beta_txt, object$tau)
  )
}

#' @rdname plot_estimates
#' @exportS3Method ggplot2::autoplot gmm_result
autoplot.gmm_result <- function(object, ...) {
  p <- plot_estimates(object, ...)
  if (!is.null(object$meta$model)) {
    p <- p + ggplot2::labs(
      subtitle = sprintf("GMM model: %s", object$meta$model)
    )
  }
  p
}

#' @rdname plot_estimates
#' @param x A fitted result object (for the `plot()` methods).
#' @export
plot.np_result <- function(x, ...) {
  p <- plot_estimates(x, ...)
  print(p)
  invisible(p)
}

#' @rdname plot_estimates
#' @export
plot.cate_result <- function(x, ...) {
  p <- plot_estimates(x, ...)
  print(p)
  invisible(p)
}

#' @rdname plot_estimates
#' @export
plot.plqr_result <- function(x, ...) {
  p <- autoplot.plqr_result(x, ...)
  print(p)
  invisible(p)
}

#' @rdname plot_estimates
#' @export
plot.gmm_result <- function(x, ...) {
  p <- autoplot.gmm_result(x, ...)
  print(p)
  invisible(p)
}
