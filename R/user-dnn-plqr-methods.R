# Methods for plqr_result objects.

.plqr_result_coef_table <- function(object, level = 0.95) {
  if (!is.numeric(level) || length(level) != 1L || !is.finite(level) ||
      level <= 0 || level >= 1) {
    .abort_bad_arg("level", "must be a numeric scalar in (0,1)")
  }

  theta <- as.matrix(object$theta)
  p <- nrow(theta)
  m <- ncol(theta)
  has_inference <- identical(object$meta$theta_se_method, "direct_sandwich") &&
    is.matrix(object$theta_se) && all(dim(object$theta_se) == dim(theta))
  crit <- stats::qnorm(0.5 + level / 2)

  coef_names <- rownames(theta)
  if (is.null(coef_names)) {
    coef_names <- paste0("theta_", seq_len(p))
  }
  target_names <- colnames(theta)
  if (is.null(target_names)) {
    target_names <- paste0("target_", seq_len(m))
  }

  # A z-statistic and p-value for a supplied evaluation value would read as
  # a test of an estimated coefficient; suppress them unless a future
  # root-solving mode sets theta_source = "local_root".
  supplied <- identical(object$meta$theta_source, "supplied")

  rows <- vector("list", p * m)
  idx <- 1L
  for (target in seq_len(m)) {
    for (coef_idx in seq_len(p)) {
      est <- theta[coef_idx, target]
      se_val <- if (has_inference) object$theta_se[coef_idx, target] else NA_real_
      stat <- if (!supplied && is.finite(se_val) && se_val > 0) est / se_val else NA_real_
      pval <- if (is.finite(stat)) 2 * stats::pnorm(-abs(stat)) else NA_real_
      rows[[idx]] <- data.frame(
        target = target_names[[target]],
        term = coef_names[[coef_idx]],
        estimate = est,
        std.error = se_val,
        statistic = stat,
        p.value = pval,
        conf.low = if (is.finite(se_val)) est - crit * se_val else NA_real_,
        conf.high = if (is.finite(se_val)) est + crit * se_val else NA_real_,
        stringsAsFactors = FALSE
      )
      idx <- idx + 1L
    }
  }
  do.call(rbind, rows)
}

#' Print a `plqr_result` object (experimental)
#'
#' @section Lifecycle:
#' **Experimental.** The PLQR API and inference interfaces may change.
#' @param x An object returned by [plqr_dnn()].
#' @param ... Not used.
#' @return `x`, invisibly.
#' @seealso [plqr_dnn()].
#' @export
print.plqr_result <- function(x, ...) {
  cat("<plqr_result>\n", sep = "")
  cat("  data  : n=", x$meta$n, ", p=", x$meta$p, ", p_d=", x$meta$p_d,
    ", m=", x$meta$m, "\n",
    sep = ""
  )
  cat("  model : tau=", x$tau, ", theta_eval=[",
    paste0(round(x$theta_eval, 4), collapse = ", "),
    "]\n",
    sep = ""
  )
  cat("  DNN   : method=", x$meta$method, ", k=", x$meta$k, ", s=", x$meta$s,
    ", distance=", x$meta$distance, "\n",
    sep = ""
  )
  if (!is.null(x$theta)) {
    cat("  theta : supplied evaluation value; theta[,1] = ",
      paste0(round(x$theta[, 1], 4), collapse = ", "),
      if (identical(x$meta$theta_se_method, "direct_sandwich")) {
        "  (coefficient SEs available)\n"
      } else {
        "  (coefficient SEs not available)\n"
      },
      sep = ""
    )
  }
  moment_note <- if (!is.null(x$moment_se)) {
    "  (moment SEs available)\n"
  } else if (identical(x$meta$theta_se_method, "direct_sandwich")) {
    "  (coefficient SEs reported above)\n"
  } else {
    "  (moment SEs not requested)\n"
  }
  cat("  moment: moment[,1] = ",
    paste0(round(x$moment[, 1], 4), collapse = ", "),
    moment_note,
    sep = ""
  )
  invisible(x)
}

#' Extract coefficients from a `plqr_result` (experimental)
#'
#' Returns the coefficient value represented by the object. For current
#' [plqr_dnn()] fits this is the supplied `theta_eval` repeated across
#' targets (`meta$theta_source = "supplied"`); no local PLQR root is
#' solved. The source is attached as `attr(., "theta_source")`.
#'
#' @section Lifecycle:
#' **Experimental.** The PLQR API and inference interfaces may change.
#' @param object A `plqr_result` object.
#' @param ... Not used.
#' @return Matrix of PLQR coefficient values with a `theta_source`
#'   attribute.
#' @seealso [plqr_dnn()].
#' @export
coef.plqr_result <- function(object, ...) {
  out <- object$theta
  attr(out, "theta_source") <- object$meta$theta_source
  out
}

#' Variance-covariance matrices for `plqr_result` (experimental)
#'
#' @section Lifecycle:
#' **Experimental.** The PLQR API and inference interfaces may change.
#' @param object A `plqr_result` object.
#' @param target Optional target index. When omitted, returns one matrix for a
#'   single-target fit and a named list for a multi-target fit.
#' @param ... Not used.
#' @return A covariance matrix or list of covariance matrices.
#' @seealso [plqr_dnn()].
#' @export
vcov.plqr_result <- function(object, target = NULL, ...) {
  if (is.null(object$theta_vcov)) {
    .abort_bad_state(
      paste0(
        "Direct PLQR variance-covariance matrices are not available. ",
        "Refit with `density` or `density_col` and `se = TRUE`."
      )
    )
  }
  if (is.null(target)) {
    if (length(object$theta_vcov) == 1L) {
      return(object$theta_vcov[[1L]])
    }
    return(object$theta_vcov)
  }
  target <- .validate_scalar_int(target, arg = "target", min = 1L)
  if (target > length(object$theta_vcov)) {
    .abort_bad_arg(
      "target",
      "must select an existing target",
      expected = sprintf("1..%d", length(object$theta_vcov)),
      actual = sprintf("target=%d", target)
    )
  }
  object$theta_vcov[[target]]
}

#' Confidence intervals for `plqr_result` (experimental)
#'
#' @details
#' Coefficient uncertainty requires a fit with `density` or `density_col`.
#' `moment_se` describes score variation and is not a coefficient standard error.
#'
#' Intervals are centered at the coefficient value the object represents:
#' for current fits the supplied `theta_eval`
#' (`meta$theta_source = "supplied"`), not an estimated local root. The
#' source is attached as `attr(., "theta_source")`.
#'
#' @section Lifecycle:
#' **Experimental.** The PLQR API and inference interfaces may change.
#' @param object A `plqr_result` object from [plqr_dnn()].
#' @param parm Optional coefficient names or indices.
#' @param level Confidence level.
#' @param target Optional target indices.
#' @param ... Not used.
#' @return A data frame with confidence intervals and a `theta_source`
#'   attribute.
#' @seealso [plqr_dnn()].
#' @export
confint.plqr_result <- function(object, parm = NULL, level = 0.95,
                                target = NULL, ...) {
  if (is.null(object$theta_vcov) || !identical(object$meta$theta_se_method, "direct_sandwich")) {
    .abort_bad_state(
      "Direct PLQR coefficient intervals are not available. Refit with `density` or `density_col` and `se = TRUE`."
    )
  }
  tab <- .plqr_result_coef_table(object, level = level)

  if (!is.null(target)) {
    target_idx <- .validate_index_vector(
      target,
      upper = ncol(object$theta),
      arg = "target",
      label = "target indices"
    )
    keep_targets <- colnames(object$theta)[target_idx]
    tab <- tab[tab$target %in% keep_targets, , drop = FALSE]
  }

  if (!is.null(parm)) {
    if (is.character(parm)) {
      if (length(parm) == 0L || anyNA(parm) || any(!nzchar(parm)) ||
          any(!(parm %in% rownames(object$theta)))) {
        .abort_bad_arg("parm", "must contain valid coefficient names or indices")
      }
      keep_terms <- parm
    } else if (is.numeric(parm)) {
      parm_idx <- .validate_index_vector(
        parm,
        upper = nrow(object$theta),
        arg = "parm",
        label = "coefficient indices"
      )
      keep_terms <- rownames(object$theta)[parm_idx]
    } else {
      .abort_bad_arg("parm", "must contain valid coefficient names or indices")
    }
    tab <- tab[tab$term %in% keep_terms, , drop = FALSE]
  }

  out <- tab[, c("target", "term", "conf.low", "conf.high"), drop = FALSE]
  rownames(out) <- NULL
  attr(out, "theta_source") <- object$meta$theta_source
  out
}

#' Summarize a `plqr_result` (experimental)
#'
#' @details
#' Coefficient uncertainty requires a fit with `density` or `density_col`.
#' `moment_se` describes score variation and is not a coefficient standard error.
#'
#' @section Lifecycle:
#' **Experimental.** The PLQR API and inference interfaces may change.
#' @param object A `plqr_result` object from [plqr_dnn()].
#' @param level Confidence level.
#' @param ... Not used.
#' @return An object of class `summary.plqr_result`.
#' @seealso [plqr_dnn()].
#' @export
summary.plqr_result <- function(object, level = 0.95, ...) {
  out <- list(
    coefficients = .plqr_result_coef_table(object, level = level),
    moment = as.matrix(object$moment),
    level = level,
    has_inference = identical(object$meta$theta_se_method, "direct_sandwich") &&
      !is.null(object$theta_vcov),
    theta_source = object$meta$theta_source,
    meta = object$meta
  )
  class(out) <- "summary.plqr_result"
  out
}

#' Print a `summary.plqr_result` (experimental)
#'
#' @section Lifecycle:
#' **Experimental.** The PLQR API and inference interfaces may change.
#' @param x A `summary.plqr_result` object.
#' @param digits Number of significant digits for printing.
#' @param ... Not used.
#' @return `x`, invisibly.
#' @seealso [plqr_dnn()].
#' @export
print.summary.plqr_result <- function(x, digits = max(3L, getOption("digits") - 3L), ...) {
  cat("<summary.plqr_result>\n", sep = "")
  cat("  targets: ", length(unique(x$coefficients$target)), "\n", sep = "")
  if (identical(x$theta_source, "supplied")) {
    cat(
      "  theta_source: supplied; no local PLQR root solving was performed\n",
      sep = ""
    )
  } else if (!is.null(x$theta_source)) {
    cat("  theta_source: ", x$theta_source, "\n", sep = "")
  }
  if (!isTRUE(x$has_inference)) {
    cat("  note   : direct coefficient standard errors are not available.\n", sep = "")
  }
  coef_tab <- x$coefficients
  num_cols <- vapply(coef_tab, is.numeric, logical(1L))
  coef_tab[num_cols] <- lapply(coef_tab[num_cols], signif, digits = digits)
  print(coef_tab)
  invisible(x)
}

#' Predict from a `plqr_result` object at new targets (experimental)
#'
#' @details
#' The returned fit matrix repeats the supplied `theta_eval`
#' (`meta$theta_source = "supplied"`). Target-specific variation enters only
#' through the DNN weights and the resulting direct-sandwich standard
#' errors; no target-specific PLQR root is solved.
#'
#' Coefficient uncertainty requires a fit with `density` or `density_col`.
#' `moment_se` describes score variation and is not a coefficient standard error.
#' Prediction standard errors always use the analytic sandwich meat,
#' regardless of the fit's `theta_meat` choice, because the fitted object
#' does not retain the outcome data needed to re-run the jackknife at new
#' targets.
#'
#' @section Lifecycle:
#' **Experimental.** The PLQR API and inference interfaces may change.
#' @param object A `plqr_result` object from [plqr_dnn()].
#' @param newdata Optional matrix/data frame of target covariates.
#' @param X0 Optional numeric matrix of targets (rows).
#' @param x0 Optional numeric vector for a single target.
#' @param se.fit Logical; if `TRUE`, also return direct coefficient standard
#'   errors. The fitted object must retain density values from `density` or
#'   `density_col`.
#' @param interval One of `"none"` or `"confidence"`.
#' @param level Confidence level for interval output.
#' @param ... Not used.
#'
#' @return With `interval = "none"` and `se.fit = FALSE`, a coefficient matrix
#'   with one row per PLQR coefficient and one column per target. With
#'   `se.fit = TRUE`, a list with `fit` and `se.fit` matrices. With
#'   `interval = "confidence"`, a list with `fit`, `lwr`, and `upr` matrices,
#'   and also `se.fit` when requested. The coefficient source is attached to
#'   the fit as `attr(., "theta_source")`.
#' @examples
#' set.seed(3)
#' demo <- plqr_demo_data(n = 80, p = 2, p_d = 2, tau = 0.5, seed = 3)
#' fit <- plqr_dnn(
#'   data = demo$data,
#'   x_cols = c("X1", "X2"),
#'   D_cols = c("D1", "D2"),
#'   g_col = "g_oracle",
#'   mD_cols = c("mD1", "mD2"),
#'   y_col = "Y",
#'   theta_eval = demo$truth$theta_true,
#'   tau = 0.5,
#'   x0 = c(0, 0),
#'   method = "exact",
#'   k = 1L, s = 5L,
#'   se = FALSE
#' )
#' predict(fit, x0 = c(0.2, -0.1))
#' @seealso [plqr_dnn()].
#' @export
predict.plqr_result <- function(
    object,
    newdata = NULL,
    X0 = NULL,
    x0 = NULL,
    se.fit = FALSE,
    interval = c("none", "confidence"),
    level = 0.95,
    ...) {
  interval <- match.arg(interval)
  level <- .validate_scalar_prob(level, arg = "level")
  se.fit <- .validate_scalar_logical(se.fit, arg = "se.fit")

  model <- .plqr_validate_prediction_model(object$.model)
  x0_mat <- .plqr_resolve_targets(model, newdata = newdata, X0 = X0, x0 = x0)

  wt <- weights_dnn(
    X = model$x_mat,
    X0 = x0_mat,
    distance = model$distance,
    method = model$method,
    k = model$k,
    s = model$s,
    normalize = model$normalize,
    ctx = model$ctx,
    return = "W"
  )
  target_names <- paste0("target_", seq_len(ncol(wt$W)))
  theta <- .plqr_prediction_theta(
    theta_eval = model$theta_eval,
    m = ncol(wt$W),
    coef_names = model$coef_names,
    target_names = target_names
  )

  attr(theta, "theta_source") <- object$meta$theta_source

  need_se <- isTRUE(se.fit) || identical(interval, "confidence")
  if (!need_se) {
    return(theta)
  }
  if (is.null(model$density)) {
    .abort_bad_state(paste0(
      "Direct PLQR coefficient standard errors are not available for prediction ",
      "because this object does not retain density values. The jackknife SE on ",
      "the moment is not a validated basis for a coefficient-level interval here ",
      "(see markdown_notes/design/15_plqr_variance_estimation_rationale.md); ",
      "refit with `density` or `density_col` for coefficient-level inference."
    ))
  }

  direct_fit <- .plqr_direct_sandwich(
    d_mat = model$d_mat,
    mD_mat = model$mD_mat,
    density_vec = model$density,
    W = wt$W,
    tau = model$tau,
    coef_names = model$coef_names,
    target_names = target_names
  )
  se_hat <- direct_fit$theta_se

  if (identical(interval, "none")) {
    return(list(fit = theta, se.fit = se_hat))
  }

  crit <- stats::qnorm(0.5 + level / 2)
  out <- list(
    fit = theta,
    lwr = theta - crit * se_hat,
    upr = theta + crit * se_hat
  )
  if (isTRUE(se.fit)) {
    out$se.fit <- se_hat
  }
  out
}
