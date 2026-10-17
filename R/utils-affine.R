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
                        if (is.null(model$normalize)) TRUE else model$normalize)
}

.affine_weight_method <- function(method) {
  tryCatch(match.arg(method, c("exact", "asymp")), error = function(e) {
    .abort_bad_arg("method", "must be one of 'exact' or 'asymp'")
  })
}
