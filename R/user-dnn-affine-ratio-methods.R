#' Reporting and prediction for paired affine-ratio results
#'
#' These methods apply to affine_ratio_result and its clate_result subclass.
#' Intervals are pointwise normal-reference intervals with frozen-score standard
#' errors. A zero or numerically unavailable ratio variance gives NA interval
#' bounds, not a zero-width confidence interval.
#'
#' @param object,x A fitted affine_ratio_result or clate_result.
#' @param target Optional single target index for vcov.
#' @param parm Optional numeric target indices for confint.
#' @param level Confidence level in (0,1).
#' @param newdata Optional reporting-target matrix/data frame. For a data-frame
#'   fit, only the retained reporting columns are extracted.
#' @param X0,x0 Target matrix or single target vector, alternatives to newdata.
#' @param se.fit Whether prediction returns standard errors.
#' @param interval Either "none" or "confidence".
#' @param digits Printing precision.
#' @param ... Not used.
#'
#' @return coef returns a named vector; vcov returns one marginal matrix or a
#'   list of marginal matrices; confint returns a two-column matrix.
#'   summary returns a summary.affine_ratio_result data frame.
#'   predict returns a vector, a fit/se.fit list or a fit/lwr/upr matrix.
#'   Printing returns its input invisibly.
#' @details Covariance matrices are on the estimator scale. Prediction changes
#'   reporting targets only, holding the observations, signals and order fixed.
#'   It never trains nuisance models. A supplied-training label does not verify
#'   the statistical or causal assumptions underlying inference.
#' @name affine_ratio_result
#' @examples
#' f <- affine_ratio_dnn(
#'   X = matrix(1:4, ncol = 1), x0 = 0,
#'   numerator = c(2, 1, 4, 3), denominator = c(1, -1, 2, 1), s = 2
#' )
#' coef(f)
#' vcov(f)
#' confint(f)
#' summary(f)
#' predict(f, x0 = 5, se.fit = TRUE)
NULL

#' @rdname affine_ratio_result
#' @export
coef.affine_ratio_result <- function(object, ...) {
  stats::setNames(object$theta, paste0("target_", seq_along(object$theta)))
}

#' @rdname affine_ratio_result
#' @export
vcov.affine_ratio_result <- function(object, target = NULL, ...) {
  if (is.null(object$theta_vcov)) {
    .abort_bad_arg("object", "has no covariance; refit with se = TRUE and an active jackknife")
  }
  if (is.null(target)) {
    if (length(object$theta_vcov) == 1L) return(object$theta_vcov[[1L]])
    return(object$theta_vcov)
  }
  target <- .affine_integer(target, "target")
  if (target > length(object$theta_vcov)) .abort_bad_arg("target", "must select an existing target")
  object$theta_vcov[[target]]
}

.affine_wald_bounds <- function(object, level) {
  level <- .validate_scalar_prob(level, "level")
  m <- length(object$theta)
  ci <- matrix(NA_real_, m, 2L, dimnames = list(
    paste0("target_", seq_len(m)), c("conf.low", "conf.high")
  ))
  if (!is.null(object$se)) {
    ok <- object$status$wald_available
    delta <- stats::qnorm((1 - level) / 2, lower.tail = FALSE) * object$se[ok]
    ci[ok, 1L] <- object$theta[ok] - delta
    ci[ok, 2L] <- object$theta[ok] + delta
  }
  ci
}

#' @rdname affine_ratio_result
#' @export
confint.affine_ratio_result <- function(object, parm = NULL, level = .95, ...) {
  level <- .validate_scalar_prob(level, "level")
  if (is.null(object$se)) {
    .abort_bad_arg("object", "has no covariance; refit with se = TRUE and an active jackknife")
  }
  idx <- if (is.null(parm)) seq_along(object$theta) else {
    .validate_index_vector(parm, length(object$theta), "parm", label = "target indices")
  }
  .affine_wald_bounds(object, level)[idx, , drop = FALSE]
}

.affine_component_se <- function(object, component) {
  if (is.null(object$component_vcov)) return(rep(NA_real_, length(object$theta)))
  vapply(object$component_vcov, function(v) {
    value <- v[component, component]
    if (!is.finite(value) || value < 0) NA_real_ else sqrt(value)
  }, numeric(1L))
}

#' @rdname affine_ratio_result
#' @export
summary.affine_ratio_result <- function(object, level = .95, ...) {
  ci <- .affine_wald_bounds(object, level)
  m <- length(object$theta)
  out <- data.frame(
    target_id = seq_len(m), numerator = object$components[, 1L],
    denominator = object$components[, 2L],
    numerator_se = .affine_component_se(object, 1L),
    denominator_se = .affine_component_se(object, 2L),
    theta = object$theta, se = if (is.null(object$se)) rep(NA_real_, m) else object$se,
    conf.low = ci[, 1L], conf.high = ci[, 2L], neff = object$neff,
    ratio_status = object$status$ratio_status,
    negative_denominator = object$status$negative_denominator,
    near_zero_denominator = object$status$near_zero_denominator,
    covariance_status = object$status$covariance_status,
    wald_available = object$status$wald_available, row.names = NULL
  )
  attr(out, "component_roles") <- object$meta$component_roles
  attr(out, "level") <- level
  class(out) <- c("summary.affine_ratio_result", "data.frame")
  out
}

#' @rdname affine_ratio_result
#' @export
print.affine_ratio_result <- function(x, digits = max(3L, getOption("digits") - 3L), ...) {
  kind <- if (inherits(x, "clate_result")) "clate_result" else "affine_ratio_result"
  cat("<", kind, "> paired exact-DNN ratio; ", length(x$theta), " target(s)\n", sep = "")
  print(summary(x), digits = digits)
  invisible(x)
}

#' @rdname affine_ratio_result
#' @export
print.summary.affine_ratio_result <- function(x, digits = max(3L, getOption("digits") - 3L), ...) {
  roles <- attr(x, "component_roles")
  if (!is.null(roles)) {
    cat("Components: numerator = ", roles[[1L]], "; denominator = ", roles[[2L]], "\n", sep = "")
  }
  counts <- table(x$ratio_status)
  cat("Ratio statuses: ", paste(paste(names(counts), counts, sep = "="), collapse = ", "), "\n", sep = "")
  print.data.frame(x, digits = digits, row.names = FALSE)
  invisible(x)
}

.affine_prediction_state <- function(object) {
  model <- object$.model
  if (!is.list(model) || is.null(model$x_mat) || is.null(model$signal) ||
      is.null(model$s) || is.null(model$row_ids)) {
    .abort_bad_arg("object$.model", "is missing frozen model state; refit with affine_ratio_dnn or clate_dnn")
  }
  model$x_mat <- .validate_matrix(model$x_mat, "object$.model$x_mat", numeric = TRUE, finite = TRUE)
  model$signal <- .validate_matrix(model$signal, "object$.model$signal", numeric = TRUE, finite = TRUE)
  if (nrow(model$signal) != nrow(model$x_mat) || ncol(model$signal) != 2L) {
    .abort_bad_arg("object$.model$signal", "must be an n by 2 matrix aligned with the reporting design")
  }
  model
}

#' @rdname affine_ratio_result
#' @export
predict.affine_ratio_result <- function(
  object, newdata = NULL, X0 = NULL, x0 = NULL, se.fit = FALSE,
  interval = c("none", "confidence"), level = .95, ...
) {
  interval <- match.arg(interval)
  se.fit <- .validate_scalar_logical(se.fit, "se.fit")
  level <- .validate_scalar_prob(level, "level")
  model <- .affine_prediction_state(object)
  if (!is.null(newdata)) {
    if (!is.null(X0) || !is.null(x0)) {
      .abort_bad_arg("newdata/X0/x0", "must supply only one target interface", quoted = FALSE)
    }
    X0 <- if (is.data.frame(newdata) && !is.null(model$x_cols)) {
      .dnn_col_mat(model$x_cols, newdata, "x_cols", out_arg = "newdata")
    } else {
      .validate_matrix(newdata, "newdata", numeric = TRUE, finite = TRUE)
    }
  }
  need_se <- se.fit || interval == "confidence"
  engine <- object$meta$jackknife
  if (is.null(engine) || !engine %in% c("fast", "loo")) engine <- "fast"
  fit <- affine_ratio_dnn(
    X = model$x_mat, X0 = X0, x0 = x0,
    numerator = model$signal[, 1L], denominator = model$signal[, 2L],
    distance = model$distance, ctx = model$ctx, s = model$s,
    se = need_se, jackknife = engine, denominator_tol = model$denominator_tol,
    row_ids = model$row_ids, return = character()
  )
  result <- if (interval == "confidence") {
    ci <- .affine_wald_bounds(fit, level)
    matrix(c(fit$theta, ci[, 1L], ci[, 2L]), ncol = 3L,
           dimnames = list(NULL, c("fit", "lwr", "upr")))
  } else fit$theta
  if (se.fit) list(fit = result, se.fit = fit$se) else result
}
