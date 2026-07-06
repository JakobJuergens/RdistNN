# Methods for gmm_result objects.

.gmm_result_coef_table <- function(object, level = 0.95) {
  if (!is.numeric(level) || length(level) != 1L || !is.finite(level) ||
      level <= 0 || level >= 1) {
    .abort_bad_arg("level", "must be a numeric scalar in (0,1)")
  }

  theta <- object$theta
  p <- nrow(theta)
  m <- ncol(theta)
  has_se <- is.matrix(object$se) && all(dim(object$se) == dim(theta))
  crit <- stats::qnorm(0.5 + level / 2)

  rows <- vector("list", p * m)
  idx <- 1L
  for (target in seq_len(m)) {
    for (coef_idx in seq_len(p)) {
      est <- theta[coef_idx, target]
      se_val <- if (has_se) object$se[coef_idx, target] else NA_real_
      stat <- if (is.finite(se_val) && se_val > 0) est / se_val else NA_real_
      pval <- if (is.finite(stat)) 2 * stats::pnorm(-abs(stat)) else NA_real_
      rows[[idx]] <- data.frame(
        target = colnames(theta)[target],
        term = rownames(theta)[coef_idx],
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

#' Print a `gmm_result` object
#'
#' @param x An object returned by [gmm_dnn()].
#' @param ... Not used.
#' @return `x`, invisibly.
#' @export
print.gmm_result <- function(x, ...) {
  cat("<gmm_result>\n", sep = "")
  cat("  data  : n=", x$meta$n, ", p=", x$meta$p, ", m=", x$meta$m,
      "\n", sep = "")
  cat("  model : ", x$meta$model, " (q=", x$meta$q,
      ", p_theta=", x$meta$p_theta, ")\n", sep = "")
  cat("  DNN   : method=", x$meta$method, ", k=", x$meta$k,
      ", s=", x$meta$s, ", distance=", x$meta$distance, "\n", sep = "")
  first_vals <- x$theta[, 1L, drop = TRUE]
  cat("  est   : theta[,1] = ",
      paste0(round(first_vals, 4), collapse = ", "),
      if (!is.null(x$se)) "  (SEs available)\n" else "  (SEs not requested)\n",
      sep = "")
  invisible(x)
}

#' Extract coefficients from a `gmm_result`
#'
#' @param object A `gmm_result` object.
#' @param ... Not used.
#' @return Matrix of estimated `theta` values.
#' @export
coef.gmm_result <- function(object, ...) {
  object$theta
}

#' Variance-covariance matrices for `gmm_result`
#'
#' @param object A `gmm_result` object.
#' @param target Optional target index.  When omitted, returns one matrix for a
#'   single-target fit and a named list for a multi-target fit.
#' @param ... Not used.
#' @return A covariance matrix or list of covariance matrices.
#' @export
vcov.gmm_result <- function(object, target = NULL, ...) {
  if (is.null(object$vcov)) {
    .abort_bad_state("Variance-covariance matrices are not available. Refit with `se = TRUE`.")
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

#' Confidence intervals for `gmm_result`
#'
#' @param object A `gmm_result` object.
#' @param parm Optional coefficient names or indices.
#' @param level Confidence level.
#' @param target Optional target indices.
#' @param ... Not used.
#' @return A data frame with confidence intervals.
#' @export
confint.gmm_result <- function(object, parm = NULL, level = 0.95,
                               target = NULL, ...) {
  tab <- .gmm_result_coef_table(object, level = level)
  if (is.null(object$se)) {
    .abort_bad_state("Standard errors are not available. Refit with `se = TRUE`.")
  }

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

#' Summarize a `gmm_result`
#'
#' @param object A `gmm_result` object from [gmm_dnn()].
#' @param level Confidence level.
#' @param ... Not used.
#' @return An object of class `summary.gmm_result`.
#' @export
summary.gmm_result <- function(object, level = 0.95, ...) {
  out <- list(
    coefficients = .gmm_result_coef_table(object, level = level),
    overid = data.frame(
      target = names(object$J),
      J = as.numeric(object$J),
      df = as.integer(object$df),
      p.value = as.numeric(object$pval_J),
      neff = as.numeric(object$neff),
      stringsAsFactors = FALSE
    ),
    convergence = object$convergence,
    level = level,
    meta = object$meta
  )
  class(out) <- "summary.gmm_result"
  out
}

#' Print a `summary.gmm_result`
#'
#' @param x A `summary.gmm_result` object.
#' @param digits Number of significant digits for printing.
#' @param ... Not used.
#' @return `x`, invisibly.
#' @export
print.summary.gmm_result <- function(x, digits = max(3L, getOption("digits") - 3L), ...) {
  cat("<summary.gmm_result>\n", sep = "")
  cat("  model  : ", x$meta$model, "\n", sep = "")
  cat("  targets: ", length(unique(x$coefficients$target)), "\n", sep = "")
  coef_tab <- x$coefficients
  num_cols <- vapply(coef_tab, is.numeric, logical(1L))
  coef_tab[num_cols] <- lapply(coef_tab[num_cols], signif, digits = digits)
  print(coef_tab)
  if (any(is.finite(x$overid$J))) {
    cat("\nOveridentification diagnostics:\n", sep = "")
    overid_tab <- x$overid
    overid_num <- vapply(overid_tab, is.numeric, logical(1L))
    overid_tab[overid_num] <- lapply(overid_tab[overid_num], signif, digits = digits)
    print(overid_tab)
  }
  invisible(x)
}
