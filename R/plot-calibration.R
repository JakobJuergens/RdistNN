# PLQR local quantile calibration panel (plan panel P-3): the weighted
# share of below-threshold observations per target against the quantile
# level tau, with a normal-approximation reference band.

#' Plot local quantile calibration for PLQR fits
#'
#' For each target, the weighted share of observations below the quantile
#' threshold, `sum_i w_ij 1{Y_i <= D_i' theta + g_i}`, is plotted against
#' the nominal level `tau` with a reference band
#' `tau +/- z * sqrt(tau (1 - tau) / n_eff)`. Shares outside the band flag
#' local miscalibration of the supplied `theta_eval` or first stages — the
#' applied reading of the quantile moment condition.
#'
#' The below-threshold indicator is stored on the fit; alternatively,
#' supply `residuals` (`Y - D' theta - g`, length `n`), whose non-positive
#' entries define the indicator.
#'
#' @param object A `plqr_result` with the weight matrix `W` stored (the
#'   default).
#' @param residuals Optional numeric vector of length `n` of quantile
#'   residuals; used when the fit does not store the indicator.
#' @param level Confidence level of the reference band (default `0.95`).
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
#'   plot_calibration(fit)
#' }
#' @export
plot_calibration <- function(object, residuals = NULL, level = 0.95, ...) {
  .plot_require_ggplot2("plot_calibration")
  if (!inherits(object, "plqr_result")) {
    .abort_bad_arg(
      "object",
      sprintf(
        "plot_calibration() supports plqr_result fits; got class '%s'",
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

  below <- if (is.list(model) && !is.null(model$below_threshold)) {
    model$below_threshold
  } else if (!is.null(residuals)) {
    r <- as.numeric(residuals)
    if (length(r) != n || any(!is.finite(r))) {
      .abort_bad_arg(
        "residuals",
        sprintf("must be a finite numeric vector of length n = %d", n),
        quoted = FALSE
      )
    }
    as.numeric(r <= 0)
  } else {
    .abort_bad_arg(
      "residuals",
      paste0(
        "is required: this fit predates stored calibration indicators; ",
        "refit, or supply `residuals` (Y - D'theta - g, length n)"
      ),
      quoted = FALSE
    )
  }

  tau <- object$tau
  share <- as.numeric(crossprod(below, W))
  neff <- 1 / colSums(W^2)
  z <- stats::qnorm(1 - (1 - level) / 2)
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
