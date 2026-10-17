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
  if (is.null(row_ids)) return(seq_len(n))
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

