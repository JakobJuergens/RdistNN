#' Paired component and first-stage diagnostics
#'
#' @param object An affine_ratio_result or clate_result.
#' @param candidates "default", "all", or candidate ids from a grid-mode fit.
#'   Candidate diagnostics do not provide candidate-selection inference.
#' @param ... Not used.
#' @details Reports component uncertainty, denominator diagnostics, geometry and
#'   centered score pseudo-value concentration. Instrument and propensity fields
#'   are available only from CLATE observations/predictions. Training declarations
#'   are descriptive. Effective information is not a complier count.
#' @return A dnn_diagnostics data frame. Component covariance is estimator-scale.
#'   component_covariance reports the paired numerator-denominator covariance;
#'   geometry_rule and tie_rule identify the frozen geometry and tie convention.
#'   Zero denominator variance gives NA standardized magnitude/correlation.
#'   Disabled covariance gives unavailable uncertainty/concentration fields.
#'   physical_radius_90 uses original-coordinate Euclidean distance; radius_90
#'   uses the fitted localization metric. Both enclose 90 percent of weight.
#' @examples
#' f <- affine_ratio_dnn(
#'   X = matrix(1:4, ncol = 1), x0 = 0,
#'   numerator = c(2, 1, 4, 3), denominator = c(1, -1, 2, 1), s = 2
#' )
#' diagnostics(f)
#' @name clate_diagnostics
NULL

#' @rdname clate_diagnostics
#' @export
diagnostics.affine_ratio_result <- function(object, candidates = "default", ...) {
  if (!.plot_candidates_is_default(candidates)) {
    return(.dnn_candidate_diagnostics(object, candidates))
  }
  model <- .affine_prediction_state(object)
  geometry <- model$geometry
  if (is.null(geometry$ranks)) {
    .abort_bad_arg("object$.model", "is missing the frozen geometry; refit the estimator")
  }
  ranks <- .validate_permutation_matrix(geometry$ranks, "object$.model$geometry$ranks")
  rw <- .affine_rank_weight(model)
  W <- rank_weight_eval_batch(rw, ranks)
  m <- ncol(W)
  out <- .dnn_weight_diagnostics(W, object$meta, geometry$dx, object$meta$tuning)
  ns <- .affine_component_se(object, 1L)
  ds <- .affine_component_se(object, 2L)
  denom_stat <- correlation <- cross_covariance <- rep(NA_real_, m)
  variance_status <- rep("not_computed", m)
  if (!is.null(object$component_vcov)) {
    for (j in seq_len(m)) {
      cross_covariance[j] <- object$component_vcov[[j]][1L, 2L]
      variance_status[j] <- if (!is.finite(ds[j])) "numerical_failure" else if (ds[j] == 0) "zero" else "positive"
      if (is.finite(ds[j]) && ds[j] > 0) denom_stat[j] <- object$components[j, 2L] / ds[j]
      if (is.finite(ns[j]) && ns[j] > 0 && is.finite(ds[j]) && ds[j] > 0) {
        correlation[j] <- object$component_vcov[[j]][1L, 2L] / ns[j] / ds[j]
      }
    }
  }
  top1 <- top5 <- rep(NA_real_, m)
  if (!is.null(object$component_vcov) && model$s < nrow(model$x_mat)) {
    loo <- .dnn_loo_centered(
      model$signal, ranks,
      .rank_weight_table(rw, nrow(model$x_mat) - 1L),
      tab_full = .rank_weight_table(rw, nrow(model$x_mat))
    )$loo
    for (j in seq_len(m)) {
      if (!is.finite(object$theta[j])) next
      centered <- sweep(loo[, , j], 2L, colMeans(loo[, , j]), "-")
      q <- object$theta[j] * centered[, 2L] - centered[, 1L]
      scale <- max(abs(q))
      if (!is.finite(scale) || scale == 0) next
      squares <- (q / scale)^2
      shares <- sort(squares / sum(squares), decreasing = TRUE)
      top1[j] <- shares[1L]
      top5[j] <- sum(utils::head(shares, 5L))
    }
  }
  iv <- model$iv
  iv_available <- !is.null(iv)
  mass1 <- clipped_mass <- rep(NA_real_, m)
  if (iv_available) {
    mass1 <- as.numeric(crossprod(iv$A, W))
    clipped_mass <- as.numeric(crossprod(as.numeric(iv$propensity$clipped), W))
  }
  roles <- object$meta$component_roles
  extra <- list(
    numerator = object$components[, 1L], denominator = object$components[, 2L],
    numerator_se = ns, denominator_se = ds,
    component_covariance = cross_covariance,
    component_correlation = correlation, denominator_stat = denom_stat,
    denominator_variance_status = variance_status,
    theta = object$theta, se = if (is.null(object$se)) rep(NA_real_, m) else object$se,
    neff = object$neff,
    physical_radius_90 = .dnn_physical_radius(W, model$x_mat, as.matrix(object$x0)),
    radius_90 = if (is.null(geometry$dx)) {
      rep(NA_real_, m)
    } else {
      vapply(seq_len(m), function(j) .dnn_weighted_quantile(geometry$dx[, j], W[, j], .9), numeric(1L))
    },
    jk_top1_sq_share = top1, jk_top5_sq_share = top5,
    score_concentration_available = is.finite(top1),
    instrument_available = rep(iv_available, m),
    instrument_1_mass = mass1, instrument_0_mass = if (iv_available) 1 - mass1 else rep(NA_real_, m),
    propensity_clipped_mass = clipped_mass,
    training = rep(if (iv_available) iv$training else "supplied", m),
    numerator_role = rep(roles[[1L]], m), denominator_role = rep(roles[[2L]], m),
    geometry_rule = rep("frozen_full_sample", m), tie_rule = rep(object$meta$tie_rule, m)
  )
  for (nm in setdiff(names(object$status), "target_id")) extra[[nm]] <- object$status[[nm]]
  .dnn_diag_append(out, extra)
}

#' @rdname clate_diagnostics
#' @export
diagnostics.clate_result <- function(object, candidates = "default", ...) {
  diagnostics.affine_ratio_result(object, candidates = candidates, ...)
}
