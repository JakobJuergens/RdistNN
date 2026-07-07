# Methods for plqr_result objects.

.plqr_result_coef_table <- function(object, level = 0.95) {
  if (!is.numeric(level) || length(level) != 1L || !is.finite(level) ||
      level <= 0 || level >= 1) {
    .abort_bad_arg("level", "must be a numeric scalar in (0,1)")
  }

  theta <- as.matrix(object$theta)
  p <- nrow(theta)
  m <- ncol(theta)
  has_inference <- identical(object$meta$se_route, "direct_sandwich") &&
    is.matrix(object$se) && all(dim(object$se) == dim(theta))
  crit <- stats::qnorm(0.5 + level / 2)

  coef_names <- rownames(theta)
  if (is.null(coef_names)) {
    coef_names <- paste0("theta_", seq_len(p))
  }
  target_names <- colnames(theta)
  if (is.null(target_names)) {
    target_names <- paste0("target_", seq_len(m))
  }

  rows <- vector("list", p * m)
  idx <- 1L
  for (target in seq_len(m)) {
    for (coef_idx in seq_len(p)) {
      est <- theta[coef_idx, target]
      se_val <- if (has_inference) object$se[coef_idx, target] else NA_real_
      stat <- if (is.finite(se_val) && se_val > 0) est / se_val else NA_real_
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

#' Print a `plqr_result` object
#'
#' @param x An object returned by [plqr_dnn()].
#' @param ... Not used.
#' @return `x`, invisibly.
#' @export
print.plqr_result <- function(x, ...) {
  cat("<plqr_result>\n", sep = "")
  cat("  data  : n=", x$meta$n, ", p=", x$meta$p, ", p_d=", x$meta$p_d,
    ", m=", x$meta$m, "\n",
    sep = ""
  )
  cat("  model : tau=", x$tau, ", beta=[",
    paste0(round(x$beta, 4), collapse = ", "),
    "]\n",
    sep = ""
  )
  cat("  DNN   : method=", x$meta$method, ", k=", x$meta$k, ", s=", x$meta$s,
    ", distance=", x$meta$distance, "\n",
    sep = ""
  )
  if (!is.null(x$theta)) {
    cat("  est   : theta[,1] = ",
      paste0(round(x$theta[, 1], 4), collapse = ", "),
      if (identical(x$meta$se_route, "direct_sandwich")) {
        "  (coefficient SEs available)\n"
      } else {
        "  (coefficient SEs not available)\n"
      },
      sep = ""
    )
  }
  moment_note <- if (!is.null(x$moment_se)) {
    "  (moment SEs available)\n"
  } else if (identical(x$meta$se_route, "direct_sandwich")) {
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

#' Extract coefficients from a `plqr_result`
#'
#' @param object A `plqr_result` object.
#' @param ... Not used.
#' @return Matrix of PLQR coefficient values.
#' @export
coef.plqr_result <- function(object, ...) {
  object$theta
}

#' Variance-covariance matrices for `plqr_result`
#'
#' @param object A `plqr_result` object.
#' @param target Optional target index. When omitted, returns one matrix for a
#'   single-target fit and a named list for a multi-target fit.
#' @param ... Not used.
#' @return A covariance matrix or list of covariance matrices.
#' @export
vcov.plqr_result <- function(object, target = NULL, ...) {
  if (is.null(object$vcov)) {
    .abort_bad_state(
      "Direct PLQR variance-covariance matrices are not available. Refit with `density` or `density_col` and `se = TRUE`."
    )
  }
  if (is.null(target)) {
    if (length(object$vcov) == 1L) {
      return(object$vcov[[1L]])
    }
    return(object$vcov)
  }
  target <- .validate_scalar_int(target, arg = "target", min = 1L)
  if (target > length(object$vcov)) {
    .abort_bad_arg(
      "target",
      "must select an existing target",
      expected = sprintf("1..%d", length(object$vcov)),
      actual = sprintf("target=%d", target)
    )
  }
  object$vcov[[target]]
}

#' Confidence intervals for `plqr_result`
#'
#' @param object A `plqr_result` object from [plqr_dnn()].
#' @param parm Optional coefficient names or indices.
#' @param level Confidence level.
#' @param target Optional target indices.
#' @param ... Not used.
#' @return A data frame with confidence intervals.
#' @export
confint.plqr_result <- function(object, parm = NULL, level = 0.95,
                                target = NULL, ...) {
  if (is.null(object$vcov) || !identical(object$meta$se_route, "direct_sandwich")) {
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
  out
}

#' Summarize a `plqr_result`
#'
#' @param object A `plqr_result` object from [plqr_dnn()].
#' @param level Confidence level.
#' @param ... Not used.
#' @return An object of class `summary.plqr_result`.
#' @export
summary.plqr_result <- function(object, level = 0.95, ...) {
  out <- list(
    coefficients = .plqr_result_coef_table(object, level = level),
    moment = as.matrix(object$moment),
    level = level,
    has_inference = identical(object$meta$se_route, "direct_sandwich") &&
      !is.null(object$vcov),
    meta = object$meta
  )
  class(out) <- "summary.plqr_result"
  out
}

#' Print a `summary.plqr_result`
#'
#' @param x A `summary.plqr_result` object.
#' @param digits Number of significant digits for printing.
#' @param ... Not used.
#' @return `x`, invisibly.
#' @export
print.summary.plqr_result <- function(x, digits = max(3L, getOption("digits") - 3L), ...) {
  cat("<summary.plqr_result>\n", sep = "")
  cat("  targets: ", length(unique(x$coefficients$target)), "\n", sep = "")
  if (!isTRUE(x$has_inference)) {
    cat("  note   : direct coefficient standard errors are not available.\n", sep = "")
  }
  coef_tab <- x$coefficients
  num_cols <- vapply(coef_tab, is.numeric, logical(1L))
  coef_tab[num_cols] <- lapply(coef_tab[num_cols], signif, digits = digits)
  print(coef_tab)
  invisible(x)
}
