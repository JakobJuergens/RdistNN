# Data extractors for the plotting layer (Phase 0).
#
# Design contract: every extractor returns a plain data.frame that the
# ggplot2-based builders (later phases) consume, and that users can consume
# directly to build their own figures. No plotting dependency is used here.
# See markdown_notes/plotting_layer_implementation_plan.md.

# ── Internal helpers ──────────────────────────────────────────────────────────

.plot_data_se_qxm <- function(se, q, m) {
  if (is.null(se)) {
    return(matrix(NA_real_, nrow = q, ncol = m))
  }
  se_mat <- if (is.matrix(se)) se else matrix(as.numeric(se), nrow = q, ncol = m)
  if (nrow(se_mat) != q || ncol(se_mat) != m) {
    .abort_bad_state(sprintf(
      "stored `se` has dimensions %d x %d; expected %d x %d.",
      nrow(se_mat), ncol(se_mat), q, m
    ))
  }
  se_mat
}

.plot_data_estimates_frame <- function(est, se, x0, coef_names, level) {
  level <- .validate_scalar_prob(level, arg = "level")
  est <- as.matrix(est)
  q <- nrow(est)
  m <- ncol(est)
  se_mat <- .plot_data_se_qxm(se, q = q, m = m)
  z <- stats::qnorm(1 - (1 - level) / 2)

  out <- data.frame(
    target_id = rep(seq_len(m), each = q),
    coef = rep(as.character(coef_names), times = m),
    estimate = as.numeric(est),
    se = as.numeric(se_mat),
    stringsAsFactors = FALSE
  )
  out$ci_lower <- out$estimate - z * out$se
  out$ci_upper <- out$estimate + z * out$se

  if (!is.null(x0)) {
    x0 <- as.matrix(x0)
    p <- ncol(x0)
    x0_df <- as.data.frame(matrix(
      as.numeric(x0[rep(seq_len(m), each = q), , drop = FALSE]),
      nrow = q * m, ncol = p
    ))
    names(x0_df) <- paste0("x0_", seq_len(p))
    out <- cbind(out, x0_df)
    attr(out, "x0_names") <- colnames(x0)
  }
  attr(out, "level") <- level
  out
}

# ── Estimates ─────────────────────────────────────────────────────────────────

#' Extract a tidy estimates frame from a fitted RdistNN object
#'
#' Returns one row per target point and coefficient with the point estimate,
#' stored jackknife standard error, and symmetric normal confidence bounds.
#' This is the data layer of the plotting interface: plot builders consume
#' these frames, and users can consume them directly for custom figures.
#'
#' @param object A fitted result object (`np_result`, `cate_result`,
#'   `plqr_result`, or `gmm_result`; `tdnn_np_result` dispatches through
#'   `np_result`).
#' @param level Confidence level for the normal bounds (default `0.95`).
#' @param ... Passed to methods; not used.
#'
#' @return A `data.frame` with columns `target_id` (integer, `1..m`),
#'   `coef` (character; `"mu"`/`"tau"` for scalar fits, moment or
#'   coefficient names for vector fits), `estimate`, `se`, `ci_lower`,
#'   `ci_upper`, and — when the fit stores evaluation points — target
#'   coordinates `x0_1..x0_p`, with the original covariate names recorded in
#'   `attr(., "x0_names")`. Fits without stored standard errors return `NA`
#'   in `se` and the bound columns. The confidence level is recorded in
#'   `attr(., "level")`.
#'
#' @examples
#' set.seed(1)
#' n <- 60
#' x_mat <- matrix(rnorm(n * 2), n, 2)
#' y_vec <- x_mat[, 1] + rnorm(n, sd = 0.3)
#' fit <- np_dnn(
#'   X = x_mat, X0 = rbind(c(0, 0), c(1, 0)), Y = y_vec,
#'   k = 1L, s = 10L, se = FALSE
#' )
#' plot_data_estimates(fit)
#' @export
plot_data_estimates <- function(object, level = 0.95, ...) {
  UseMethod("plot_data_estimates")
}

#' @rdname plot_data_estimates
#' @export
plot_data_estimates.default <- function(object, level = 0.95, ...) {
  .abort_bad_arg(
    "object",
    sprintf(
      "plot_data_estimates() is not implemented for class '%s'",
      class(object)[1L]
    ),
    quoted = FALSE
  )
}

#' @rdname plot_data_estimates
#' @export
plot_data_estimates.np_result <- function(object, level = 0.95, ...) {
  mu <- matrix(as.numeric(object$mu), nrow = 1L)
  .plot_data_estimates_frame(
    est = mu, se = object$se, x0 = object$x0,
    coef_names = "mu", level = level
  )
}

#' @rdname plot_data_estimates
#' @export
plot_data_estimates.cate_result <- function(object, level = 0.95, ...) {
  tau <- matrix(as.numeric(object$tau), nrow = 1L)
  .plot_data_estimates_frame(
    est = tau, se = object$se, x0 = object$x0,
    coef_names = "tau", level = level
  )
}

#' @rdname plot_data_estimates
#' @export
plot_data_estimates.plqr_result <- function(object, level = 0.95, ...) {
  est <- as.matrix(object$moment)
  coef_names <- rownames(est)
  if (is.null(coef_names)) {
    coef_names <- paste0("moment_", seq_len(nrow(est)))
  }
  .plot_data_estimates_frame(
    est = est, se = object$se, x0 = object$x0,
    coef_names = coef_names, level = level
  )
}

#' @rdname plot_data_estimates
#' @export
plot_data_estimates.gmm_result <- function(object, level = 0.95, ...) {
  est <- as.matrix(object$theta)
  coef_names <- object$meta$coef_names
  if (is.null(coef_names)) coef_names <- rownames(est)
  if (is.null(coef_names)) coef_names <- paste0("theta_", seq_len(nrow(est)))
  .plot_data_estimates_frame(
    est = est, se = object$se, x0 = object$x0,
    coef_names = coef_names, level = level
  )
}

# ── Weights ───────────────────────────────────────────────────────────────────

#' Extract a tidy weight-profile frame from a fitted RdistNN object
#'
#' Returns one row per (target point, rank) with the DNN weight, cumulative
#' weight share, the originating observation index, and — when distances were
#' stored — the distance to the target. Ranks are assigned by ascending
#' stored distance when `dx` is available and by descending weight otherwise;
#' the two orderings coincide for monotone single-scale DNN weights, while
#' two-scale (TDNN) fits should request `return = c("W", "dx")` for
#' distance-based ranks. For exact DNN fits the finite-population support edge
#' is rank `n - s + k`, not rank `s`; `s` is the subsample/localization order.
#'
#' @param object A fitted result object holding the weight matrix `W`
#'   (stored by default via `return = c("W")`).
#'
#' @return A `data.frame` with columns `target_id` (integer), `rank`
#'   (integer, `1..n`), `obs_id` (integer row index into the training data),
#'   `weight`, `cum_share` (cumulative weight share in rank order), and
#'   `distance` (`NA` when distances were not stored).
#'
#' @examples
#' set.seed(1)
#' n <- 60
#' x_mat <- matrix(rnorm(n * 2), n, 2)
#' y_vec <- x_mat[, 1] + rnorm(n, sd = 0.3)
#' fit <- np_dnn(
#'   X = x_mat, x0 = c(0, 0), Y = y_vec,
#'   k = 1L, s = 10L, se = FALSE, return = c("W", "dx")
#' )
#' head(plot_data_weights(fit))
#' @export
plot_data_weights <- function(object) {
  W <- object$W
  if (is.null(W)) {
    .abort_bad_arg(
      "object",
      "must store the weight matrix `W`; refit with `return = c(\"W\")` (add \"dx\" for distance-based ranks)",
      quoted = FALSE
    )
  }
  W <- as.matrix(W)
  n <- nrow(W)
  m <- ncol(W)
  dx <- object$dx
  if (!is.null(dx)) dx <- as.matrix(dx)

  frames <- vector("list", m)
  for (j in seq_len(m)) {
    ord <- if (!is.null(dx)) {
      order(dx[, j], -W[, j])
    } else {
      order(W[, j], decreasing = TRUE)
    }
    w_j <- W[ord, j]
    total <- sum(w_j)
    cum_share <- if (is.finite(total) && abs(total) > 0) {
      cumsum(w_j) / total
    } else {
      rep(NA_real_, n)
    }
    frames[[j]] <- data.frame(
      target_id = rep(j, n),
      rank = seq_len(n),
      obs_id = as.integer(ord),
      weight = as.numeric(w_j),
      cum_share = as.numeric(cum_share),
      distance = if (is.null(dx)) rep(NA_real_, n) else as.numeric(dx[ord, j]),
      stringsAsFactors = FALSE
    )
  }
  do.call(rbind, frames)
}

# ── Diagnostics ───────────────────────────────────────────────────────────────

#' Extract the per-target diagnostics frame from a fitted RdistNN object
#'
#' A tidy wrapper around [diagnostics()]: the same per-target frame
#' (`n_eff`, weight concentration, distance radii, tuning provenance, and
#' estimator-specific columns), returned as a plain `data.frame` with an
#' integer `target_id` consistent with the other `plot_data_*()` extractors.
#'
#' @param object A fitted result object with a [diagnostics()] method.
#'
#' @return A `data.frame` with one row per target point; the original
#'   diagnostics label is kept in `target_label`.
#'
#' @examples
#' set.seed(1)
#' n <- 60
#' x_mat <- matrix(rnorm(n * 2), n, 2)
#' y_vec <- x_mat[, 1] + rnorm(n, sd = 0.3)
#' fit <- np_dnn(
#'   X = x_mat, x0 = c(0, 0), Y = y_vec,
#'   k = 1L, s = 10L, se = FALSE
#' )
#' plot_data_diagnostics(fit)
#' @export
plot_data_diagnostics <- function(object) {
  d <- diagnostics(object)
  out <- as.data.frame(d, stringsAsFactors = FALSE)
  class(out) <- "data.frame"
  if ("target_id" %in% names(out)) {
    out$target_label <- as.character(out$target_id)
  }
  out$target_id <- seq_len(nrow(out))
  out
}
