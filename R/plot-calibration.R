# PLQR local quantile calibration panel (plan panel P-3): the weighted
# share of below-threshold observations per target against the quantile
# level tau, with a normal-approximation reference band.

#' Plot local quantile calibration for PLQR fits (experimental)
#'
#' For each target, the weighted share of observations below the quantile
#' threshold, `sum_i w_ij 1{Y_i <= D_i' theta + g_i}`, is plotted against
#' the nominal level `tau` with a reference band
#' `tau +/- z * sqrt(tau (1 - tau) / n_eff)`. This is a descriptive calibration
#' diagnostic for the supplied threshold; the band does not test the full
#' coefficient moment condition or account for fitted-nuisance uncertainty.
#'
#' The below-threshold indicator is taken from the fit's retained state.
#' Refit if this indicator is missing or invalid.
#'
#' @section Lifecycle:
#' **Experimental.** The PLQR API and inference interfaces may change.
#' @param object A `rdistnn_plqr_result` with the weight matrix `W` stored (the
#'   default).
#' @param level Confidence level of the reference band (default `0.95`).
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
#'   plot_calibration(fit)
#' }
#' @export
plot_calibration <- function(object, level = 0.95) {
  .dnn_reuse_scope(object)
  .plot_require_ggplot2("plot_calibration")
  if (!inherits(object, "rdistnn_plqr_result")) {
    .abort_bad_arg(
      "object",
      sprintf(
        "plot_calibration() supports rdistnn_plqr_result fits; got class '%s'",
        class(object)[1L]
      ),
      quoted = FALSE
    )
  }
  level <- .validate_scalar_prob(level, arg = "level")

  W <- object$W
  if (is.null(W)) {
    .abort_bad_arg(
      "object",
      "must store the weight matrix `W` for calibration; refit with `return = c(\"W\")`",
      quoted = FALSE
    )
  }
  W <- as.matrix(W)
  n <- nrow(W)
  model <- object$.model

  below <- model$below_threshold
  if (!is.numeric(below) || is.complex(below) || !is.null(dim(below)) ||
      length(below) != n || any(!is.finite(below)) || any(!below %in% c(0, 1))) {
    .dnn_abort("state_invalid", paste0(
      "Retained `object$.model$below_threshold` must be a 0/1 numeric vector aligned with the observations. ",
      "Refit with plqr_dnn()."
    ))
  }

  tau <- object$tau
  share <- as.numeric(crossprod(below, W))
  neff <- 1 / colSums(W^2)
  z <- .dnn_normal_critical(level)
  half <- z * sqrt(tau * (1 - tau) / neff)

  dd <- .plot_diag_frame(object)
  axis <- .plot_profile_axis(dd)
  xcol <- if (!is.null(axis$var)) axis$var else "target_id"
  df <- data.frame(
    x = if (!is.null(axis$var)) dd[[axis$var]] else seq_along(share),
    share = share,
    band_lower = tau - half,
    band_upper = tau + half,
    stringsAsFactors = FALSE
  )
  n_out <- sum(share < df$band_lower | share > df$band_upper)

  p <- ggplot2::ggplot(df, ggplot2::aes(x = .data$x)) +
    ggplot2::geom_ribbon(
      ggplot2::aes(ymin = .data$band_lower, ymax = .data$band_upper),
      alpha = 0.15, fill = .plot_col_single
    ) +
    ggplot2::geom_hline(yintercept = tau, linetype = "dashed") +
    ggplot2::geom_line(ggplot2::aes(y = .data$share),
                       colour = .plot_col_treated) +
    ggplot2::geom_point(ggplot2::aes(y = .data$share),
                        colour = .plot_col_treated) +
    ggplot2::labs(
      x = axis$label, y = "Weighted below-threshold share",
      caption = paste(
        c(
          sprintf(
            "reference band: tau +/- %.2f-level normal bounds; %d target(s) outside",
            level, n_out
          ),
          .plot_honesty_caption(object)
        ),
        collapse = "\n"
      )
    )
  p + .rdistnn_theme()
}
