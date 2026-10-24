# Strict array contracts shared by the affine ratio and IV signal boundaries.
.affine_vector <- function(x, n = NULL, arg) {
  if (is.null(x)) .abort_bad_arg(arg, "must be provided")
  if (!is.numeric(x) || is.complex(x) || !is.null(dim(x))) {
    .abort_bad_arg(arg, "must be a numeric vector")
  }
  if (is.null(n)) {
    if (!length(x)) .abort_bad_arg(arg, "must be nonempty")
  } else if (length(x) != n) {
    .abort_bad_arg(arg, sprintf("must have length %d; got %d", n, length(x)))
  }
  .abort_if_na_or_nonfinite(x, arg = arg, quoted = TRUE)
  as.numeric(x)
}

.affine_row_ids <- function(row_ids, n) {
  if (is.null(row_ids)) {
    return(seq_len(n))
  }
  if ((!is.character(row_ids) && !is.numeric(row_ids)) ||
    is.object(row_ids) || !is.null(dim(row_ids)) || is.complex(row_ids)) {
    .abort_bad_arg("row_ids", "must be a character or numeric vector")
  }
  if (length(row_ids) != n) {
    .abort_bad_arg("row_ids", sprintf("must have length %d; got %d", n, length(row_ids)))
  }
  if (anyNA(row_ids) || (is.numeric(row_ids) && any(!is.finite(row_ids))) ||
    (is.character(row_ids) && any(!nzchar(row_ids))) || anyDuplicated(row_ids)) {
    .abort_bad_arg("row_ids", "must contain unique, nonmissing, nonempty finite identities")
  }
  row_ids
}

.affine_nonnegative_scalar <- function(x, arg) {
  if (!is.numeric(x) || length(x) != 1L || !is.finite(x) || x < 0) {
    .abort_bad_arg(arg, "must be a finite nonnegative numeric scalar")
  }
  as.numeric(x)
}

.affine_integer <- function(x, arg) {
  if (!is.numeric(x) || length(x) != 1L || !is.finite(x) ||
    x < 1 || x > .Machine$integer.max || x != floor(x)) {
    .abort_bad_arg(arg, "must be a positive integer scalar within the R integer range")
  }
  as.integer(x)
}

.affine_column <- function(col, data, n, arg) {
  nm <- .dnn_resolve_data_cols(col, data, arg, allow_multiple = FALSE)
  .affine_vector(data[[nm]], n, arg)
}

.affine_sample_size <- function(x_mat, arg) {
  n <- nrow(x_mat)
  if (n < 1L) .abort_bad_arg(arg, "must contain at least one observation")
  n
}
.affine_rank_weight <- function(model) {
  .dnn_rank_weight_spec(model$method, model$k, model$s,
                        .affine_retained_normalize(model))
}

.affine_retained_normalize <- function(model) {
  normalize <- model$normalize
  if (!is.logical(normalize) || length(normalize) != 1L || is.na(normalize)) {
    .dnn_abort("state_invalid", paste0(
      "Retained `object$.model$normalize` must be TRUE or FALSE. ",
      "Refit with affine_ratio_dnn() or clate_dnn()."
    ))
  }
  normalize
}

.affine_weight_method <- function(method) {
  .dnn_match_arg(method, c("exact", "asymp"))
}

# Fitted and candidate rows share the same paired-ratio availability grammar.
# Callers validate the real numeric estimate/SE arrays before checking statuses.
.affine_validate_retained_status <- function(status, estimate, se, active, invalid, field = "status") {
  required <- c("ratio_status", "covariance_status", "wald_available")
  if (!is.data.frame(status) || nrow(status) != length(estimate) ||
      !all(required %in% names(status))) {
    invalid(field, "must retain aligned ratio and covariance statuses")
  }
  bad <- function(name, reason) invalid(paste0(field, "$", name), reason)
  choices <- list(ratio_status = c("ok", "component_overflow", "zero_denominator", "ratio_overflow"),
    covariance_status = c("ok", "not_computed", "numerical_failure", "zero_score_variance"))
  for (name in names(choices)) {
    value <- status[[name]]
    if (!is.character(value) || !is.null(dim(value)) || anyNA(value) ||
        any(!value %in% choices[[name]])) bad(name, "contains invalid status values")
  }
  wald <- status$wald_available
  if (!is.logical(wald) || !is.null(dim(wald)) || anyNA(wald)) {
    bad("wald_available", "must retain non-missing logical availability flags")
  }
  finite <- is.finite(estimate)
  if (any((status$ratio_status == "ok") != finite)) {
    bad("ratio_status", "is inconsistent with ratio availability")
  }
  if (any((status$covariance_status == "not_computed") != !active)) {
    bad("covariance_status", "is inconsistent with the requested uncertainty policy")
  }
  if (is.null(se)) se <- rep(NA_real_, length(estimate))
  positive <- is.finite(se) & se > 0
  expected <- active & finite & positive & status$covariance_status == "ok"
  if (any(wald != expected) || any(finite & status$covariance_status == "ok" & !positive) ||
      any(status$covariance_status == "zero_score_variance" & (is.na(se) | se != 0)) ||
      any(status$covariance_status == "numerical_failure" & !is.na(se))) {
    bad("wald_available", "is inconsistent with ratio and covariance availability")
  }
  invisible(NULL)
}
