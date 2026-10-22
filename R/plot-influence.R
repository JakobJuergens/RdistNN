# CATE influence panel (plan panel C-5): per-observation weighted score
# contributions, flagging single-observation dominance. This is the
# plot-level counterpart of the delete-one jackknife's sensitivity logic.

#' Plot per-observation contributions to a CATE or CLATE estimating equation
#'
#' For a target `j`, each observation contributes
#' `w_ij (psi_i - tau_hat_j)` to the localized estimating equation; large
#' contributions mean the estimate leans on few observations. Drawn
#' against rank (or stored distance when available), with the largest
#' three absolute contributions labelled by observation index and the
#' weight-concentration diagnostics in the caption.
#'
#' The pseudo-outcome is taken from the fit; supply `pseudo_outcome` for
#' fits that do not store it.
#' For CLATE the contribution is w_ij (phi_Yi - theta_j phi_Di), using both
#' retained signals at the target's fitted ratio. These are weighted score
#' residuals, not delete-one pseudo-values or estimated individual effects.
#'
#' @param object A `cate_result` or `clate_result` with the weight matrix `W` stored (the
#'   default).
#' @param targets Optional integer vector of target ids (defaults to the
#'   first/median-`n_eff`/last subset above eight targets).
#' @param pseudo_outcome Optional numeric vector of length `n`; used when
#'   the fit does not store the pseudo-outcome.
#' @param n_label Integer; number of largest absolute contributions to
#'   label per target (default `3`).
#' @param ... Ignored; reserved for future options.
#'
#' @return A ggplot object.
#'
#' @examples
#' if (requireNamespace("ggplot2", quietly = TRUE)) {
#'   set.seed(1)
#'   n <- 100
#'   X <- matrix(rnorm(n * 2), n, 2)
#'   psi <- rnorm(n, mean = X[, 1])
#'   fit <- cate_dnn(
#'     X = X, x0 = c(0, 0), pseudo_outcome = psi,
#'     k = 1L, s = 15L, se = FALSE, return = c("W", "dx")
#'   )
#'   plot_influence(fit)
#' }
#' @export
plot_influence <- function(object,
                           targets = NULL,
                           pseudo_outcome = NULL,
                           n_label = 3L,
                           ...) {
  .dnn_reuse_scope(object)
  .plot_require_ggplot2("plot_influence")
  clate <- inherits(object, "clate_result")
  if (!inherits(object, "cate_result") && !clate) {
    .abort_bad_arg(
      "object",
      sprintf(
        "plot_influence() supports cate_result or clate_result fits; got class '%s'",
        class(object)[1L]
      ),
      quoted = FALSE
    )
  }
  n_label <- .validate_scalar_int(n_label, "n_label")

  model <- object$.model
  if (clate && !is.null(pseudo_outcome)) {
    .abort_bad_arg("pseudo_outcome", "is a CATE input; CLATE uses its paired stored signals")
  }
  psi <- if (clate) {
    .affine_prediction_state(object)$signal
  } else if (is.list(model) && !is.null(model$psi)) {
    model$psi
  } else {
    pseudo_outcome
  }
  if (is.null(psi)) {
    .abort_bad_arg(
      "pseudo_outcome",
      "is required: the fit does not store the pseudo-outcome",
      quoted = FALSE
    )
  }

  dw <- plot_data_weights(object)
  n <- max(dw$rank)
  if (!clate) psi <- as.numeric(psi)
  if ((!clate && length(psi) != n) || (clate && nrow(psi) != n) || any(!is.finite(psi))) {
    .abort_bad_arg(
      "pseudo_outcome",
      sprintf("must be a finite numeric vector of length n = %d", n),
      quoted = FALSE
    )
  }

  all_targets <- sort(unique(dw$target_id))
  neff_by_target <- vapply(
    all_targets,
    function(j) 1 / sum(dw$weight[dw$target_id == j]^2),
    numeric(1L)
  )
  targets <- .plot_select_targets(all_targets, neff_by_target, targets)
  dw <- dw[dw$target_id %in% targets, , drop = FALSE]

  if (clate) {
    if (any(!is.finite(object$theta[targets]))) {
      .abort_bad_arg("targets", paste0(
        "must select targets with a finite CLATE ratio; inspect ratio_status in diagnostics(object)"
      ))
    }
  }
  dw$contribution <- 0
  for (j in targets) {
    rows <- which(dw$target_id == j)
    obs <- dw$obs_id[rows]
    w <- dw$weight[rows]
    left <- if (clate) psi[obs, 1L] else psi[obs]
    right <- if (clate) numeric(length(rows)) else rep(object$tau[j], length(rows))
    if (clate) {
      active <- which(w > 0)
      right[active] <- object$theta[j] * psi[obs[active], 2L]
    }
    dw$contribution[rows] <- .dnn_score_contributions(w, left, right)$contribution
  }
  if (any(!is.finite(dw$contribution))) {
    .abort_bad_arg("object", paste0(
      "has unrepresentable weighted score contributions; rescale the retained signals"
    ))
  }
  has_dx <- !all(is.na(dw$distance))
  dw$xval <- if (has_dx) dw$distance else dw$rank
  dw$target <- factor(dw$target_id)
  multi <- length(targets) > 1L

  lab_rows <- do.call(rbind, lapply(split(dw, dw$target_id), function(dj) {
    dj[order(-abs(dj$contribution)), , drop = FALSE][
      seq_len(min(n_label, nrow(dj))), , drop = FALSE
    ]
  }))

  mapping <- if (multi) {
    ggplot2::aes(x = .data$xval, y = .data$contribution,
                 colour = .data$target)
  } else {
    ggplot2::aes(x = .data$xval, y = .data$contribution)
  }
  point_args <- if (multi) {
    list(alpha = 0.7)
  } else {
    list(alpha = 0.7, colour = .plot_col_single)
  }

  dd <- plot_data_diagnostics(object)
  cap <- paste(
    c(
      sprintf(
        "max weight %.3f; top-1 weight share %.3f (across shown targets)",
        max(dd$max_weight[targets]), max(dd$top1_weight_share[targets])
      ),
      .plot_honesty_caption(object)
    ),
    collapse = "\n"
  )

  p <- ggplot2::ggplot(dw, mapping) +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed") +
    do.call(ggplot2::geom_point, point_args) +
    ggplot2::geom_text(
      data = lab_rows,
      ggplot2::aes(label = .data$obs_id),
      vjust = -0.6, size = 3, show.legend = FALSE
    ) +
    ggplot2::labs(
      x = if (has_dx) "Distance" else "Rank",
      y = if (clate) "Weighted CLATE score contribution" else "Weighted score contribution",
      colour = "target",
      caption = cap
    )
  p + .rdistnn_theme()
}
