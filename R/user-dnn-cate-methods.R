# Methods for rdistnn_cate_result objects.

#' Joint jackknife covariance across CATE reporting profiles
#'
#' @param object A fitted `rdistnn_cate_result` from [cate_dnn()].
#' @param ... Not used.
#' @return An `m` by `m` covariance matrix, with axes in `coef(object)` order.
#' @details Uses common deleted observations with frozen supplied signals.
#'   The diagonal agrees with the fit's squared standard errors. The `fast`
#'   route freezes fitted geometry; `loo` retains deletion-specific geometry
#'   refitting. Covariance is calculated on request from retained model state,
#'   independently of optional `return` extras, including after serialization.
#'   Fitted uncertainty and an active jackknife are required. Singular matrices
#'   and exact zero covariance are retained. Joint covariance supplies the
#'   arithmetic for prespecified contrasts; valid contrast inference requires
#'   the corresponding joint assumptions and does not supply simultaneous bands.
#' @examples
#' fit <- cate_dnn(X = matrix(1:7, ncol = 1),
#'   X0 = matrix(c(2, 5), ncol = 1),
#'   pseudo_outcome = c(2, 1, 4, 3, 6, 5, 7), s = 3,
#'   return = character())
#' vcov(fit)
#' @export
vcov.rdistnn_cate_result <- function(object, ...) {
  .dnn_reuse_scope(object)
  jackknife <- object$meta$jackknife
  if (is.null(object$se) || is.null(jackknife) ||
      !(jackknife %in% c("fast", "loo"))) {
    .dnn_abort("extra_unavailable", "Covariance is not available. Refit with `se = TRUE` and an active jackknife.")
  }
  model <- .cate_validate_prediction_model(object$.model)
  required <- if (identical(jackknife, "loo")) {
    c("inference_distance", "inference_ctx")
  } else {
    c("distance", "ctx")
  }
  if (!all(required %in% names(model))) {
    .dnn_abort("state_invalid", "This object lacks fitted geometry for covariance. Refit with the current cate_dnn().")
  }
  x0_mat <- .dnn_resolve_targets(model$x_mat, X0 = object$x0)
  loo <- .dnn_jackknife_loo_array(
    x_mat = model$x_mat, x0_mat = x0_mat, signal = model$psi,
    distance = model[[required[[1L]]]], ctx = model[[required[[2L]]]],
    method = model$method, k = model$k, s = model$s,
    normalize = model$normalize, signal_arg = "pseudo_outcome",
    jackknife = jackknife,
    precomputed_ranks = if (identical(jackknife, "fast")) .dnn_retained_ranking(object)$ranks else NULL,
    precomputed_groups = model$groups
  )
  n <- nrow(model$x_mat)
  V <- .covariance_matrix(matrix(loo, nrow = n), (n - 1) / n)
  dimnames(V) <- rep(list(names(stats::coef(object))), 2L)
  V
}

#' Summarize a `rdistnn_cate_result` object
#'
#' @param object A `rdistnn_cate_result` object from [cate_dnn()].
#' @param level Confidence level used for interval columns.
#' @param ... Not used.
#'
#' @return An object of class `rdistnn_summary_cate_result` with a coefficient table
#'   and `meta`, including the original fit's signal provenance and creation
#'   version.
#' @examples
#' set.seed(1)
#' n <- 30
#' x_mat <- matrix(rnorm(n * 2), n, 2)
#' out <- cate_dnn(X = x_mat, x0 = c(0, 0), pseudo_outcome = rnorm(n))
#' summary(out)
#' @export
summary.rdistnn_cate_result <- function(object, level = 0.95, ...) {
  .dnn_check_report_state(object)
  if (!is.numeric(level) || length(level) != 1L || !is.finite(level) || level <= 0 || level >= 1) {
    .abort_bad_arg("level", "must be a numeric scalar in (0,1)")
  }

  est <- as.numeric(object$tau)
  m <- length(est)
  se <- object$se
  has_se <- is.numeric(se) && length(se) == m && all(is.finite(se))

  z <- pval <- lwr <- upr <- rep(NA_real_, m)
  se_col <- rep(NA_real_, m)
  if (has_se) {
    se_col <- as.numeric(se)
    z <- est / se_col
    pval <- 2 * stats::pnorm(-abs(z))
    crit <- .dnn_normal_critical(level)
    lwr <- est - crit * se_col
    upr <- est + crit * se_col
  }

  coef_tab <- data.frame(
    estimate = est,
    std.error = se_col,
    statistic = z,
    p.value = pval,
    conf.low = lwr,
    conf.high = upr,
    row.names = paste0("target_", seq_len(m))
  )

  out <- list(
    coefficients = coef_tab,
    level = as.numeric(level),
    has_se = has_se,
    meta = .cate_result_meta(object$meta)
  )
  class(out) <- "rdistnn_summary_cate_result"
  out
}

#' Print a `rdistnn_summary_cate_result` object
#'
#' @param x A `rdistnn_summary_cate_result` object.
#' @param digits Number of significant digits for printing.
#' @param ... Not used.
#'
#' @return `x`, invisibly.
#' @export
print.rdistnn_summary_cate_result <- function(x, digits = max(3L, getOption("digits") - 3L), ...) {
  .dnn_check_tie_policy(x$meta$tie_rule)
  meta <- .cate_result_meta(x$meta)
  cat("<rdistnn_summary_cate_result>\n", sep = "")
  cat("  signal: ", meta$signal_source, "\n", sep = "")
  .dnn_print_localization(x$meta)
  cat("  targets: ", nrow(x$coefficients), "\n", sep = "")
  if (!isTRUE(x$has_se)) {
    cat("  note   : standard errors are not available in this object.\n", sep = "")
  }
  print(signif(x$coefficients, digits = digits))
  invisible(x)
}

#' Confidence intervals for `rdistnn_cate_result`
#'
#' @param object A `rdistnn_cate_result` object from [cate_dnn()].
#' @param parm Optional integer vector selecting target indices.
#' @param level Confidence level.
#' @param ... Not used.
#'
#' @return A matrix with lower and upper confidence bounds.
#' @examples
#' set.seed(2)
#' n <- 25
#' x_mat <- matrix(rnorm(n * 2), n, 2)
#' out <- cate_dnn(X = x_mat, x0 = c(0, 0), pseudo_outcome = rnorm(n), se = TRUE)
#' confint(out)
#' @export
confint.rdistnn_cate_result <- function(object, parm = NULL, level = 0.95, ...) {
  .dnn_check_report_state(object)
  if (!is.numeric(level) || length(level) != 1L || !is.finite(level) || level <= 0 || level >= 1) {
    .abort_bad_arg("level", "must be a numeric scalar in (0,1)")
  }
  if (is.null(object$se)) {
    .dnn_abort("extra_unavailable", "Standard errors are not available. Refit with `se = TRUE`.")
  }

  m <- length(object$tau)
  idx <- if (is.null(parm)) {
    seq_len(m)
  } else {
    .validate_index_vector(
      parm,
      upper = m,
      arg = "parm",
      label = "target indices"
    )
  }

  est <- as.numeric(object$tau[idx])
  se <- as.numeric(object$se[idx])
  crit <- .dnn_normal_critical(level)
  ci <- cbind(
    conf.low = est - crit * se,
    conf.high = est + crit * se
  )
  rownames(ci) <- paste0("target_", idx)
  ci
}

#' Predict from a `rdistnn_cate_result` object at new targets
#'
#' Uses retained observations, supplied signals and fitted geometry, including
#' after serialization and without optional `return` extras. Prediction
#' uncertainty uses frozen-geometry rank deletions, even if the fit used
#' `jackknife = "loo"`. It can be requested from a fit with `se = FALSE`
#' when `s < n`.
#' Prediction localizes the retained signal without establishing new
#' inferential guarantees.
#'
#' @param object A `rdistnn_cate_result` object from [cate_dnn()].
#' @param newdata Optional matrix/data frame of target covariates.
#' @param X0 Optional numeric matrix of targets (rows).
#' @param x0 Optional numeric vector for a single target.
#' @param se.fit Logical; if `TRUE`, also return jackknife standard errors.
#' @param interval One of `"none"` or `"confidence"`.
#' @param level Confidence level for interval output.
#' @param ... Not used.
#'
#' @return A numeric vector when no uncertainty is requested; a list with
#'   `fit` and `se.fit` when `se.fit = TRUE`. With `interval = "confidence"`,
#'   `fit` is a matrix with columns `fit`, `lwr`, `upr`, returned directly
#'   unless `se.fit = TRUE`.
#' @examples
#' set.seed(3)
#' n <- 40
#' x_mat <- matrix(rnorm(n * 2), n, 2)
#' fit <- cate_dnn(X = x_mat, x0 = c(0, 0), pseudo_outcome = rnorm(n), se = TRUE)
#' predict(fit, x0 = c(0.2, -0.1))
#' @export
predict.rdistnn_cate_result <- function(
    object,
    newdata = NULL,
    X0 = NULL,
    x0 = NULL,
    se.fit = FALSE,
    interval = c("none", "confidence"),
    level = 0.95,
    ...) {
  .dnn_reuse_scope(object, groups = NULL)
  interval <- .dnn_match_arg(interval)
  if (!is.numeric(level) || length(level) != 1L || !is.finite(level) || level <= 0 || level >= 1) {
    .abort_bad_arg("level", "must be a numeric scalar in (0,1)")
  }
  if (!is.logical(se.fit) || length(se.fit) != 1L || is.na(se.fit)) {
    .abort_bad_arg("se.fit", "must be TRUE or FALSE")
  }

  model <- .cate_validate_prediction_model(object$.model)

  x0_mat <- .dnn_resolve_prediction_targets(model, newdata = newdata, X0 = X0, x0 = x0)

  need_se <- se.fit || interval == "confidence"
  wt_return <- if (need_se) c("W", "ranks") else "W"
  wt <- weights_dnn(
    X = model$x_mat,
    X0 = x0_mat,
    distance = model$distance,
    method = model$method,
    k = model$k,
    s = model$s,
    normalize = model$normalize,
    ctx = model$ctx,
    return = wt_return
  )
  fit <- as.numeric(.dnn_estimate_from_signal(model$psi, wt$W,
    n = nrow(model$x_mat), arg = "pseudo_outcome"))

  if (!need_se) {
    return(fit)
  }

  if (as.integer(model$s) >= nrow(model$x_mat)) {
    .abort_bad_arg(
      "model$s",
      "must satisfy s < n when jackknife-based prediction intervals are requested",
      expected = sprintf("s < n=%d", nrow(model$x_mat)),
      actual = sprintf("s=%d", as.integer(model$s)),
      quoted = FALSE
    )
  }
  se_hat <- .cate_dnn_jackknife_se_loo(
    x_mat = model$x_mat,
    x0_mat = x0_mat,
    psi = model$psi,
    distance = model$distance,
    method = model$method,
    k = as.integer(model$k),
    s = as.integer(model$s),
    normalize = isTRUE(model$normalize),
    ctx = model$ctx,
    jackknife = "fast",
    precomputed_ranks = wt$ranks, precomputed_groups = wt$groups
  )

  if (interval == "none") {
    if (!se.fit) {
      return(fit)
    }
    return(list(fit = fit, se.fit = se_hat))
  }

  crit <- .dnn_normal_critical(level)
  pred_mat <- cbind(
    fit = fit,
    lwr = fit - crit * se_hat,
    upr = fit + crit * se_hat
  )
  if (!se.fit) {
    return(pred_mat)
  }
  list(fit = pred_mat, se.fit = se_hat)
}

#' Extract coefficients from a `rdistnn_cate_result`
#'
#' @param object A `rdistnn_cate_result` object.
#' @param ... Not used.
#' @return Named numeric vector of point estimates (`object$tau`), one entry
#'   per target point.
#' @export
coef.rdistnn_cate_result <- function(object, ...) {
  .dnn_check_report_state(object)
  stats::setNames(as.numeric(object$tau), paste0("target_", seq_along(object$tau)))
}

#' Print a `rdistnn_cate_result` object
#'
#' @param x An object returned by [cate_dnn()].
#' @param ... Not used.
#' @return `x`, invisibly.
#' @examples
#' set.seed(1)
#' n <- 20
#' x_mat <- matrix(rnorm(n * 2), n, 2)
#' out <- cate_dnn(
#'   X = x_mat,
#'   x0 = c(0, 0),
#'   pseudo_outcome = rnorm(n),
#'   se = FALSE
#' )
#' print(out)
#' @export
print.rdistnn_cate_result <- function(x, ...) {
  .dnn_check_report_state(x)
  meta <- .cate_result_meta(x$meta)
  cat("<rdistnn_cate_result>\n", sep = "")
  cat("  signal: ", meta$signal_source, "\n", sep = "")
  .dnn_print_localization(x$meta)
  cat("  data  : n=", x$meta$n, ", p=", x$meta$p, ", m=", x$meta$m, "\n", sep = "")
  cat("  DNN   : method=", x$meta$method, ", k=", x$meta$k,
    ", distance=", x$meta$distance, "\n",
    sep = ""
  )
  if (!is.null(x$se)) {
    cat("  est   : tau[1:min(5,m)] = ",
      paste0(round(utils::head(x$tau, 5), 4), collapse = ", "),
      "  (SEs available)\n",
      sep = ""
    )
  } else {
    cat("  est   : tau[1:min(5,m)] = ",
      paste0(round(utils::head(x$tau, 5), 4), collapse = ", "),
      "  (SEs not requested)\n",
      sep = ""
    )
  }
  invisible(x)
}
