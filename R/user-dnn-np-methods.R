# Methods for np_result objects.

#' Print an `np_result` object
#'
#' @param x An object returned by [np_dnn()].
#' @param ... Not used.
#' @return `x`, invisibly.
#' @export
print.np_result <- function(x, ...) {
  .dnn_check_report_state(x)
  cat("<np_result>\n", sep = "")
  .dnn_print_localization(x$meta)
  cat("  data  : n=", x$meta$n, ", p=", x$meta$p, ", m=", x$meta$m, "\n", sep = "")
  cat("  DNN   : method=", x$meta$method, ", k=", x$meta$k,
    ", distance=", x$meta$distance, "\n",
    sep = ""
  )
  if (!is.null(x$se)) {
    cat("  est   : mu[1:min(5,m)] = ",
      paste0(round(utils::head(x$mu, 5), 4), collapse = ", "),
      "  (SEs available)\n",
      sep = ""
    )
  } else {
    cat("  est   : mu[1:min(5,m)] = ",
      paste0(round(utils::head(x$mu, 5), 4), collapse = ", "),
      "  (SEs not requested)\n",
      sep = ""
    )
  }
  invisible(x)
}

.np_wald_bounds <- function(mu, se, level) {
  ci <- matrix(NA_real_, length(mu), 2L,
    dimnames = list(paste0("target_", seq_along(mu)), c("conf.low", "conf.high")))
  if (!is.null(se)) {
    ok <- is.finite(se) & se >= 0 & is.finite(mu)
    delta <- .dnn_normal_critical(level) * se[ok]
    ci[ok, ] <- cbind(mu[ok] - delta, mu[ok] + delta)
  }
  ci
}

#' Reporting and prediction for nonparametric regression results
#'
#' Methods for `np_result`, including the `tdnn_np_result` subclass returned
#' by [np_tdnn()]. Estimates are conditional means, not CATE coefficients.
#'
#' @param object An object from [np_dnn()] or [np_tdnn()].
#' @param x A summary object for printing.
#' @param level Confidence level in (0,1).
#' @param parm Optional integer vector selecting target indices.
#' @param newdata Optional matrix/data frame of target covariates. Data-frame
#'   fits extract the retained covariate column names, allowing extra columns.
#' @param X0,x0 Alternative target matrix or single target vector.
#' @param se.fit Logical; return jackknife standard errors with predictions.
#' @param interval One of `"none"` or `"confidence"`.
#' @param digits Printing precision.
#' @param ... Not used.
#'
#' @details Fits retain canonical outcomes, training covariates, weight controls
#'   and the resolved distance specification and fitted context independently
#'   of optional `return` extras. Prediction uses this full training state and
#'   fitted geometry, including after serialization; it never fits on newdata.
#'   Custom distance functions must respect their retained context to freeze
#'   fitted preprocessing. A named distance is resolved at fitting time.
#'
#'   `vcov()` returns the full joint covariance across targets from common
#'   delete-one observations. Its diagonal agrees with the fit's squared SEs.
#'   The `jackknife = "loo"` fit and covariance recompute each deletion's
#'   geometry using the original context. Prediction uncertainty always uses
#'   frozen training geometry and rank deletions, as does `jackknife = "fast"`.
#'   Shared scaled covariance arithmetic recovers representable results after
#'   intermediate overflow or underflow. Unrepresentable nonzero variances
#'   raise a rescaling error instead of giving infinite or false-zero SEs.
#'
#'   Intervals are pointwise normal-reference intervals. A zero SE gives
#'   collapsed Wald bounds and, for a nonzero mean, an infinite statistic and
#'   zero p-value; zero divided by zero gives undefined statistics/p-values.
#'   Zero estimated variance does not establish negligible bias or coverage.
#'   Missing fitted uncertainty gives `NA` summary columns and an
#'   error from `confint()`/`vcov()`. Prediction can request uncertainty even if
#'   the original fit used `se = FALSE`, provided the orders permit deletion.
#'
#' @return `coef()` returns a named vector of local means; `vcov()` an `m x m`
#'   covariance matrix in target order; `confint()` a two-column matrix.
#'   `summary()` returns class `summary.np_result` with `coefficients`, `level`,
#'   `has_se` and `meta`. Prediction returns a vector, a `fit`/`se.fit` list,
#'   or a `fit`/`lwr`/`upr` matrix (inside the list when `se.fit = TRUE`).
#'   Printing returns its input invisibly.
#' @name np_result
#' @examples
#' fit <- np_dnn(X = matrix(1:7, ncol = 1), X0 = matrix(c(2, 5), ncol = 1),
#'               Y = c(2, 1, 4, 3, 6, 5, 7), s = 3, return = character())
#' coef(fit)
#' summary(fit)
#' vcov(fit)
#' confint(fit)
#' predict(fit, x0 = 4, se.fit = TRUE)
NULL

#' @rdname np_result
#' @export
coef.np_result <- function(object, ...) {
  .dnn_check_report_state(object)
  stats::setNames(as.numeric(object$mu), paste0("target_", seq_along(object$mu)))
}

#' @rdname np_result
#' @export
vcov.np_result <- function(object, ...) {
  .dnn_reuse_scope(object)
  if (is.null(object$se)) {
    .abort_bad_state("Covariance is not available. Refit with `se = TRUE` and an active jackknife.")
  }
  model <- .np_prediction_state(object)
  ranking <- if (identical(object$meta$jackknife, "fast")) .dnn_retained_ranking(object) else NULL
  .np_joint_covariance(model, object$x0, ranking$ranks, object$meta$jackknife, model$groups)
}

#' @rdname np_result
#' @export
summary.np_result <- function(object, level = .95, ...) {
  .dnn_check_report_state(object)
  level <- .validate_scalar_prob(level, "level")
  mu <- as.numeric(object$mu)
  se <- object$se
  m <- length(mu)
  z <- p <- rep(NA_real_, m)
  if (!is.null(se)) {
    ok <- is.finite(se) & se >= 0 & is.finite(mu)
    z[ok] <- mu[ok] / se[ok]
    p[ok] <- 2 * stats::pnorm(-abs(z[ok]))
  }
  ci <- .np_wald_bounds(mu, se, level)
  tab <- data.frame(estimate = mu, std.error = if (is.null(se)) rep(NA_real_, m) else se,
    statistic = z, p.value = p, conf.low = ci[, 1L], conf.high = ci[, 2L],
    row.names = rownames(ci))
  structure(list(coefficients = tab, level = level, has_se = !is.null(se), meta = object$meta),
    class = "summary.np_result")
}

#' @rdname np_result
#' @export
print.summary.np_result <- function(x, digits = max(3L, getOption("digits") - 3L), ...) {
  .dnn_check_tie_policy(x$meta$tie_rule)
  cat("<summary.np_result>\n")
  .dnn_print_localization(x$meta)
  if (!isTRUE(x$has_se)) cat("  Standard errors are not available.\n")
  print(signif(x$coefficients, digits))
  invisible(x)
}

#' @rdname np_result
#' @export
confint.np_result <- function(object, parm = NULL, level = .95, ...) {
  .dnn_check_report_state(object)
  level <- .validate_scalar_prob(level, "level")
  if (is.null(object$se)) {
    .abort_bad_state("Standard errors are not available. Refit with `se = TRUE` and an active jackknife.")
  }
  idx <- if (is.null(parm)) {
    seq_along(object$mu)
  } else {
    .validate_index_vector(parm, length(object$mu), "parm", label = "target indices")
  }
  .np_wald_bounds(object$mu, object$se, level)[idx, , drop = FALSE]
}

#' @rdname np_result
#' @export
predict.np_result <- function(object, newdata = NULL, X0 = NULL, x0 = NULL,
    se.fit = FALSE, interval = c("none", "confidence"), level = .95, ...) {
  .dnn_reuse_scope(object, groups = NULL)
  interval <- match.arg(interval)
  se.fit <- .validate_scalar_logical(se.fit, "se.fit")
  level <- .validate_scalar_prob(level, "level")
  model <- .np_prediction_state(object)
  x0_mat <- .dnn_resolve_prediction_targets(model, newdata, X0, x0)
  need_se <- se.fit || interval == "confidence"
  wt <- .np_prediction_weights(model, x0_mat, if (need_se) c("W", "ranks") else "W")
  fit <- as.numeric(.dnn_estimate_from_signal(model$y_vec, wt$W, nrow(model$x_mat), "Y"))
  if (!need_se) return(fit)
  se <- sqrt(pmax(.np_marginal_variance(model, x0_mat, wt$ranks, groups = wt$groups), 0))
  if (interval == "confidence") {
    ci <- .np_wald_bounds(fit, se, level)
    fit <- cbind(fit = fit, lwr = ci[, 1L], upr = ci[, 2L])
  }
  if (se.fit) return(list(fit = fit, se.fit = unname(se)))
  fit
}
