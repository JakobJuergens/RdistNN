# PLQR Gram-conditioning panel (plan panel P-2): smallest eigenvalue and
# condition number of the weighted residualized-design Gram matrix across
# targets, from the stored PLQR diagnostics.

#' Plot local design conditioning for PLQR fits (experimental)
#'
#' Two facets across targets: the smallest eigenvalue of the weighted
#' residualized-design Gram matrix (linear scale, guide at zero) and its
#' condition number (log10 scale, dashed guides at `cond_limits`, labelled
#' caution and ill-conditioned). This is the software face of the local
#' design-conditioning diagnostics: a well-conditioned local Gram is what
#' makes the moment evaluation (and later inversion) trustworthy.
#'
#' @section Lifecycle:
#' **Experimental.** PLQR functionality is outside the CATE and regression
#' workflow focus of the initial public release. Its API and inference
#' interfaces may change in later releases.
#' @param object A `plqr_result` from [plqr_dnn()].
#' @param cond_limits Numeric length-2 vector (default `c(1e3, 1e6)`);
#'   dashed guides drawn at these condition numbers.
#' @param ... Ignored; reserved for future options.
#'
#' @return A ggplot object.
#'
#' @examples
#' if (requireNamespace("ggplot2", quietly = TRUE)) {
#'   demo <- plqr_demo_data(n = 90, p = 2, p_d = 2, tau = 0.5, seed = 1)
#'   df <- demo$data
#'   fit <- plqr_dnn(
#'     X = as.matrix(df[paste0("X", 1:2)]),
#'     D = as.matrix(df[paste0("D", 1:2)]),
#'     g = as.numeric(df$g_oracle),
#'     mD = as.matrix(df[paste0("mD", 1:2)]),
#'     Y = as.numeric(df$Y),
#'     theta_eval = demo$truth$theta_true, tau = 0.5,
#'     X0 = rbind(c(0, 0), c(0.5, 0.5)), k = 1L, s = 8L, se = FALSE
#'   )
#'   plot_conditioning(fit)
#' }
#' @export
plot_conditioning <- function(object, cond_limits = c(1e3, 1e6), ...) {
  .plot_require_ggplot2("plot_conditioning")
  if (!inherits(object, "plqr_result")) {
    .abort_bad_arg(
      "object",
      sprintf(
        "plot_conditioning() supports plqr_result fits; got class '%s'",
        class(object)[1L]
      ),
      quoted = FALSE
    )
  }
  if (!is.numeric(cond_limits) || length(cond_limits) != 2L ||
      any(!is.finite(cond_limits)) || any(cond_limits <= 0) ||
      cond_limits[1L] >= cond_limits[2L]) {
    .abort_bad_arg("cond_limits", "must be two increasing positive numbers")
  }

  dd <- .plot_diag_frame(object)
  if (!"gram_available" %in% names(dd) || !any(dd$gram_available)) {
    .abort_bad_state(
      paste0(
        "no Gram conditioning diagnostics are available for this fit; ",
        "the weighted residualized design could not be formed."
      )
    )
  }
  axis <- .plot_profile_axis(dd)
  xcol <- if (!is.null(axis$var)) axis$var else "target_id"

  avail <- dd[dd$gram_available, , drop = FALSE]
  # Free facet scales cannot mix axis transforms, so the condition-number
  # facet carries log10 values with an explicit facet label.
  metric_levels <- c("smallest eigenvalue", "log10 condition number")
  long <- rbind(
    data.frame(
      x = avail[[xcol]], metric = metric_levels[1L],
      value = avail$gram_min_eig, stringsAsFactors = FALSE
    ),
    data.frame(
      x = avail[[xcol]], metric = metric_levels[2L],
      value = log10(avail$gram_condition_number), stringsAsFactors = FALSE
    )
  )
  long$metric <- factor(long$metric, levels = metric_levels)

  guides <- data.frame(
    metric = factor(
      c(metric_levels[1L], metric_levels[2L], metric_levels[2L]),
      levels = metric_levels
    ),
    yintercept = c(0, log10(cond_limits)),
    stringsAsFactors = FALSE
  )

  p <- ggplot2::ggplot(long, ggplot2::aes(x = .data$x, y = .data$value)) +
    ggplot2::geom_hline(
      data = guides,
      ggplot2::aes(yintercept = .data$yintercept),
      linetype = "dashed", colour = "grey40"
    ) +
    ggplot2::geom_line(colour = .plot_col_single) +
    ggplot2::geom_point(colour = .plot_col_single) +
    ggplot2::facet_wrap(ggplot2::vars(.data$metric), scales = "free_y") +
    ggplot2::labs(
      x = axis$label, y = NULL,
      caption = paste(
        c(
          sprintf(
            "condition-number guides at %g (caution) and %g (ill-conditioned)",
            cond_limits[1L], cond_limits[2L]
          ),
          .plot_honesty_caption(object)
        ),
        collapse = "\n"
      )
    )

  p + .rdistnn_theme()
}
