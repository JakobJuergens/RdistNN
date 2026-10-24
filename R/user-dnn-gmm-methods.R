# Methods for rdistnn_gmm_result objects.

.gmm_check_report_arrays <- function(object, n, m, invalid, require_fields, check_array) {
  require_fields(c("theta", "theta_first", "se", "J", "df", "pval_J", "neff", "convergence"))
  meta <- object$meta
  for (field in c("p_theta", "q", "k", "s")) {
    value <- meta[[field, exact = TRUE]]
    if (!is.numeric(value) || is.complex(value) || length(value) != 1L ||
        !is.finite(value) || value < 1 || value > .Machine$integer.max || value != trunc(value)) {
      invalid(paste0("meta$", field), "must be a positive integer scalar")
    }
  }
  p <- meta$p_theta
  q <- meta$q
  targets <- paste0("target_", seq_len(m))
  for (field in c("coef_names", "moment_names", "model", "distance", "method", "optimizer", "jacobian")) {
    value <- meta[[field, exact = TRUE]]
    size <- if (field == "coef_names") p else if (field == "moment_names") q else 1L
    if (!is.character(value) || !is.null(dim(value)) || length(value) != size ||
        anyNA(value) || any(!nzchar(value))) {
      invalid(paste0("meta$", field), "must retain the fitted nonempty names or selector")
    }
  }
  if (!meta$method %in% c("exact", "asymp") ||
      !meta$optimizer %in% c("BFGS", "Nelder-Mead", "L-BFGS-B") ||
      !meta$jacobian %in% c("auto", "analytic", "fd") ||
      meta$k > meta$s || meta$s > n) {
    invalid("meta", "contains invalid fitted GMM controls")
  }
  if (!is.logical(meta$se) || length(meta$se) != 1L || is.na(meta$se)) {
    invalid("meta$se", "must retain TRUE or FALSE")
  }
  check_names <- function(value, field, rows = NULL) {
    aligned <- if (is.null(rows)) identical(names(value), targets) else
      identical(dimnames(value), list(rows, targets))
    if (!aligned) invalid(field, "must retain aligned coefficient and target names")
  }
  for (field in c("theta", "theta_first")) {
    check_array(object[[field]], field, c(p, m), allow_na = TRUE)
    check_names(object[[field]], field, meta$coef_names)
  }
  missing <- colSums(is.na(object$theta))
  failed <- missing == p
  if (any(!missing %in% c(0L, p)) ||
      !identical(is.na(object$theta), is.na(object$theta_first))) {
    invalid("theta/theta_first", "must retain complete coefficient columns or unavailable failed targets")
  }
  if (meta$se) {
    check_array(object$se, "se", c(p, m), allow_na = TRUE, nonnegative = TRUE)
    check_names(object$se, "se", meta$coef_names)
    require_fields("vcov")
    if (!is.list(object$vcov) || length(object$vcov) != m) {
      invalid("vcov", "must retain one covariance matrix or NULL per target")
    }
    check_names(object$vcov, "vcov")
    for (j in seq_len(m)) {
      value <- object$vcov[[j]]
      if (is.null(value)) {
        if (any(!is.na(object$se[, j]))) invalid("se/vcov", "contains inconsistent unavailable uncertainty")
      } else {
        check_array(value, "vcov", c(p, p))
        if (failed[j] || any(diag(value) < 0) ||
            !identical(dimnames(value), list(meta$coef_names, meta$coef_names)) ||
            !isTRUE(all.equal(unname(object$se[, j]), sqrt(unname(diag(value)))))) {
          invalid("se/vcov", "contains inconsistent variances, coefficient names or uncertainty")
        }
      }
    }
  } else if (!is.null(object$se) || !is.null(object$vcov)) {
    invalid("se/vcov/meta$se", "contains inconsistent uncertainty availability")
  }
  for (field in c("J", "df", "pval_J", "neff")) {
    value <- object[[field]]
    check_array(value, field, m, allow_na = TRUE, nonnegative = TRUE)
    check_names(value, field)
    unavailable <- if (field == "neff") failed else failed | q <= p
    if (!identical(unname(is.na(value)), unname(unavailable)) ||
        (field == "df" && any(value[!is.na(value)] != q - p)) ||
        (field == "pval_J" && any(value[!is.na(value)] > 1)) ||
        (field == "neff" && any(value[!is.na(value)] <= 0))) {
      invalid(field, "contains inconsistent GMM diagnostic values or availability")
    }
  }
  conv <- object$convergence
  fields <- c("target", "convergence_first", "convergence_second", "value_first",
    "value_second", "omega_ridge", "omega_condition", "message")
  if (!is.data.frame(conv) || nrow(conv) != m) {
    invalid("convergence", "must retain one diagnostics row per target")
  }
  require_fields(fields, conv, "convergence$")
  if (!identical(conv$target, targets) || !is.character(conv$message) ||
      length(conv$message) != m || any(failed & (is.na(conv$message) | !nzchar(conv$message)))) {
    invalid("convergence", "must retain aligned target identities and numerical failure reasons")
  }
  for (field in fields[2:7]) {
    value <- conv[[field]]
    check_array(value, paste0("convergence$", field), m, allow_na = TRUE, nonnegative = TRUE)
    if (any(!is.na(value[failed])) ||
        (field %in% fields[2:3] && (anyNA(value[!failed]) ||
          any(value[!is.na(value)] != trunc(value[!is.na(value)]))))) {
      invalid(paste0("convergence$", field), "contains inconsistent optimizer diagnostics")
    }
  }
  for (field in intersect(c("omega", "jacobian", "g_obs"), names(object))) {
    values <- object[[field]]
    if (!is.list(values) || length(values) != m) invalid(field, "must retain one array per target")
    check_names(values, field)
    shape <- switch(field, omega = c(q, q), jacobian = c(q, p), g_obs = c(n, q))
    for (j in seq_len(m)) {
      if (failed[j]) {
        if (!is.null(values[[j]])) invalid(field, "must be NULL at failed targets")
      } else if (field == "jacobian" && (!meta$se ||
          (all(is.na(object$se[, j])) && is.null(object$vcov[[j]])))) {
        # Differencing/averaging finite callback values can overflow. Retain
        # that diagnostic only with unavailable or deliberately disabled SEs.
        value <- values[[j]]
        if (!is.matrix(value) || !is.numeric(value) || is.complex(value) ||
            !identical(dim(value), as.integer(shape)) || any(is.na(value) & !is.nan(value))) {
          invalid(field, "must retain an aligned real numeric arithmetic diagnostic")
        }
      } else {
        check_array(values[[j]], field, shape)
      }
    }
  }
  invisible(NULL)
}

.gmm_result_coef_table <- function(object, level = 0.95) {
  if (!is.numeric(level) || length(level) != 1L || !is.finite(level) ||
      level <= 0 || level >= 1) {
    .abort_bad_arg("level", "must be a numeric scalar in (0,1)")
  }

  theta <- object$theta
  p <- nrow(theta)
  m <- ncol(theta)
  has_se <- is.matrix(object$se) && all(dim(object$se) == dim(theta))
  crit <- .dnn_normal_critical(level)

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

#' Print a `rdistnn_gmm_result` object
#'
#' @param x An object returned by [gmm_dnn()].
#' @param ... Not used.
#' @return `x`, invisibly.
#' @export
print.rdistnn_gmm_result <- function(x, ...) {
  .dnn_check_report_state(x)
  cat("<rdistnn_gmm_result>\n", sep = "")
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

#' Extract coefficients from a `rdistnn_gmm_result`
#'
#' @param object A `rdistnn_gmm_result` object.
#' @param ... Not used.
#' @return Matrix of estimated `theta` values.
#' @export
coef.rdistnn_gmm_result <- function(object, ...) {
  .dnn_check_report_state(object)
  object$theta
}

#' Variance-covariance matrices for `rdistnn_gmm_result`
#'
#' @param object A `rdistnn_gmm_result` object.
#' @param target Optional target index.  When omitted, returns one matrix for a
#'   single-target fit and a named list for a multi-target fit.
#' @param ... Not used.
#' @return A covariance matrix or list of covariance matrices.
#' @export
vcov.rdistnn_gmm_result <- function(object, target = NULL, ...) {
  .dnn_reuse_scope(object)
  if (is.null(object$vcov)) {
    .dnn_abort("extra_unavailable", "Variance-covariance matrices are not available. Refit with `se = TRUE`.")
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

#' Confidence intervals for `rdistnn_gmm_result`
#'
#' @param object A `rdistnn_gmm_result` object.
#' @param parm Optional coefficient names or indices.
#' @param level Confidence level.
#' @param target Optional target indices.
#' @param ... Not used.
#' @return A data frame with confidence intervals.
#' @export
confint.rdistnn_gmm_result <- function(object, parm = NULL, level = 0.95,
                               target = NULL, ...) {
  .dnn_check_report_state(object)
  tab <- .gmm_result_coef_table(object, level = level)
  if (is.null(object$se)) {
    .dnn_abort("extra_unavailable", "Standard errors are not available. Refit with `se = TRUE`.")
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

#' Summarize a `rdistnn_gmm_result`
#'
#' @param object A `rdistnn_gmm_result` object from [gmm_dnn()].
#' @param level Confidence level.
#' @param ... Not used.
#' @return An object of class `rdistnn_summary_gmm_result`.
#' @export
summary.rdistnn_gmm_result <- function(object, level = 0.95, ...) {
  .dnn_check_report_state(object)
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
  class(out) <- "rdistnn_summary_gmm_result"
  out
}

#' Print a `rdistnn_summary_gmm_result`
#'
#' @param x A `rdistnn_summary_gmm_result` object.
#' @param digits Number of significant digits for printing.
#' @param ... Not used.
#' @return `x`, invisibly.
#' @export
print.rdistnn_summary_gmm_result <- function(x, digits = max(3L, getOption("digits") - 3L), ...) {
  .dnn_check_tie_policy(x$meta$tie_rule)
  cat("<rdistnn_summary_gmm_result>\n", sep = "")
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
