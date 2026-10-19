.abort_bad_arg <- function(arg, msg, expected = NULL, actual = NULL, quoted = TRUE) {
  arg_lbl <- if (isTRUE(quoted)) sprintf("`%s`", arg) else as.character(arg)
  parts <- c(sprintf("%s %s.", arg_lbl, msg))
  if (!is.null(expected)) parts <- c(parts, sprintf("Expected: %s.", expected))
  if (!is.null(actual)) parts <- c(parts, sprintf("Got: %s.", actual))
  stop(paste(parts, collapse = " "), call. = FALSE)
}

.abort_bad_state <- function(msg) {
  stop(msg, call. = FALSE)
}

.abort_bad_length <- function(arg, relation, expected = NULL, actual = NULL, quoted = TRUE) {
  .abort_bad_arg(
    arg = arg,
    msg = sprintf("length must %s", relation),
    expected = expected,
    actual = actual,
    quoted = quoted
  )
}

.has_na_only <- function(x) {
  any(is.na(x) & !is.nan(x))
}

.has_nonfinite_non_na <- function(x) {
  any(is.nan(x) | is.infinite(x))
}

.abort_if_na_or_nonfinite <- function(x, arg, quoted = FALSE) {
  if (.has_na_only(x)) {
    .abort_bad_arg(arg, "cannot contain `NA` values", quoted = quoted)
  }
  if (.has_nonfinite_non_na(x)) {
    .abort_bad_arg(arg, "cannot contain non-finite values (`Inf`, `-Inf`, or `NaN`)", quoted = quoted)
  }
  invisible(NULL)
}

.validate_scalar_int <- function(x, arg = "x", min = 1L) {
  if (!is.numeric(x) || length(x) != 1L || !is.finite(x) ||
      x < -.Machine$integer.max || x > .Machine$integer.max || x != trunc(x)) {
    .abort_bad_arg(arg, "must be an integer scalar")
  }
  x_int <- as.integer(x)
  if (x_int < as.integer(min)) {
    if (as.integer(min) == 1L) {
      .abort_bad_arg(arg, "must be a positive integer")
    }
    .abort_bad_arg(arg, sprintf("must be >= %d", as.integer(min)))
  }
  x_int
}

.validate_scalar_prob <- function(x, arg = "tau") {
  if (!is.numeric(x) || length(x) != 1L || !is.finite(x) || x <= 0 || x >= 1) {
    .abort_bad_arg(arg, "must be a numeric scalar in (0,1)")
  }
  as.numeric(x)
}

.validate_scalar_logical <- function(x, arg = "x") {
  if (!is.logical(x) || length(x) != 1L || is.na(x)) {
    .abort_bad_arg(arg, "must be TRUE or FALSE")
  }
  x
}

.validate_finite_numeric_vector <- function(x, arg = "x", min_length = 1L) {
  min_length <- .validate_scalar_int(min_length, arg = "min_length", min = 0L)
  if (!is.numeric(x)) {
    .abort_bad_arg(arg, "must be numeric")
  }
  x_vec <- as.numeric(x)
  if (length(x_vec) < min_length) {
    .abort_bad_arg(arg, sprintf("must have length at least %d", min_length))
  }
  .abort_if_na_or_nonfinite(x_vec, arg = arg)
  x_vec
}

.validate_positive_index_matrix <- function(x, arg = "ranks", min_rows = 1L, min_cols = 1L) {
  min_rows <- .validate_scalar_int(min_rows, arg = "min_rows", min = 0L)
  min_cols <- .validate_scalar_int(min_cols, arg = "min_cols", min = 0L)

  if (!is.matrix(x)) {
    .abort_bad_arg(arg, "must be a matrix")
  }
  if (!is.numeric(x)) {
    .abort_bad_arg(arg, "must be numeric or integer")
  }
  x_mat <- as.matrix(x)
  if (nrow(x_mat) < min_rows) {
    .abort_bad_arg(arg, sprintf("must have at least %d row(s)", min_rows))
  }
  if (ncol(x_mat) < min_cols) {
    .abort_bad_arg(arg, sprintf("must have at least %d column(s)", min_cols))
  }
  .abort_if_na_or_nonfinite(x_mat, arg = arg)
  if (any(x_mat != as.integer(x_mat))) {
    .abort_bad_arg(arg, "must contain integer indices")
  }
  x_int <- matrix(as.integer(x_mat), nrow = nrow(x_mat), ncol = ncol(x_mat))
  dimnames(x_int) <- dimnames(x_mat)
  if (any(x_int < 1L)) {
    .abort_bad_arg(arg, "indices must be >= 1")
  }
  x_int
}

.validate_index_vector <- function(x, upper, arg = "idx", label = "indices") {
  upper <- .validate_scalar_int(upper, arg = "upper", min = 1L)
  if (!is.numeric(x)) {
    .abort_bad_arg(arg, sprintf("must contain valid %s in 1..%d", label, upper), quoted = FALSE)
  }
  x_num <- as.numeric(x)
  if (length(x_num) == 0L) {
    .abort_bad_arg(arg, sprintf("must contain valid %s in 1..%d", label, upper), quoted = FALSE)
  }
  .abort_if_na_or_nonfinite(x_num, arg = arg, quoted = FALSE)
  if (any(x_num != floor(x_num))) {
    .abort_bad_arg(arg, sprintf("must contain integer %s", label), quoted = FALSE)
  }
  if (any(x_num < 1L | x_num > upper)) {
    .abort_bad_arg(arg, sprintf("must contain valid %s in 1..%d", label, upper), quoted = FALSE)
  }
  idx <- as.integer(x_num)
  idx
}

.validate_binary01 <- function(x, arg = "x") {
  x_vec <- as.numeric(x)
  .abort_if_na_or_nonfinite(x_vec, arg = arg, quoted = FALSE)
  if (any(!(x_vec %in% c(0, 1)))) {
    .abort_bad_arg(arg, "must be binary in {0,1}", quoted = FALSE)
  }
  as.integer(x_vec)
}

.validate_matrix <- function(
    x,
    arg = "X",
    min_rows = 1L,
    min_cols = 1L,
    numeric = TRUE,
    finite = FALSE,
    strict = FALSE) {
  min_rows <- .validate_scalar_int(min_rows, arg = "min_rows", min = 0L)
  min_cols <- .validate_scalar_int(min_cols, arg = "min_cols", min = 0L)

  if (!is.matrix(x) && isTRUE(strict)) {
    .abort_bad_arg(arg, "must be a matrix", quoted = FALSE)
  }

  x_mat <- as.matrix(x)
  if (nrow(x_mat) < min_rows) {
    .abort_bad_arg(arg, sprintf("must have at least %d row(s)", min_rows), quoted = FALSE)
  }
  if (ncol(x_mat) < min_cols) {
    .abort_bad_arg(arg, sprintf("must have at least %d column(s)", min_cols), quoted = FALSE)
  }
  if (isTRUE(numeric) && !is.numeric(x_mat)) {
    .abort_bad_arg(arg, "must be numeric", quoted = FALSE)
  }
  if (isTRUE(finite)) {
    .abort_if_na_or_nonfinite(x_mat, arg = arg, quoted = FALSE)
  }
  x_mat
}

.validate_matrix_nonempty <- function(x, arg = "X") {
  x_mat <- .validate_matrix(
    x,
    arg = arg,
    min_rows = 0L,
    min_cols = 1L,
    numeric = TRUE,
    finite = TRUE
  )
  if (ncol(x_mat) == 0L) {
    .abort_bad_arg(arg, "must have at least one column", quoted = FALSE)
  }
  x_mat
}

.validate_ncol_relation <- function(x_mat, y, y_arg = "X0", x_arg = "X", relation = c("match", "equal")) {
  relation <- match.arg(relation)
  y_mat <- as.matrix(y)
  if (ncol(x_mat) != ncol(y_mat)) {
    .abort_bad_arg(
      arg = sprintf("ncol(%s)", y_arg),
      msg = sprintf("must %s ncol(%s)", relation, x_arg),
      expected = sprintf("ncol(%s)=%d", x_arg, ncol(x_mat)),
      actual = sprintf("ncol(%s)=%d", y_arg, ncol(y_mat)),
      quoted = FALSE
    )
  }
  y_mat
}

.validate_target_vector <- function(x0, p, x0_arg = "x0", x_arg = "X", strict = FALSE, finite = FALSE) {
  if (isTRUE(strict) && !is.numeric(x0)) {
    .abort_bad_arg(x0_arg, "must be numeric", quoted = FALSE)
  }
  x0_vec <- as.numeric(x0)
  if (length(x0_vec) != p) {
    .abort_bad_length(x0_arg, sprintf("equal ncol(%s)", x_arg), quoted = FALSE)
  }
  if (isTRUE(finite)) {
    .abort_if_na_or_nonfinite(x0_vec, arg = x0_arg, quoted = FALSE)
  }
  x0_vec
}

.validate_x0_inv_s <- function(x0, inv_s, p, x_arg = "X", strict = FALSE, finite = FALSE) {
  if (isTRUE(strict) && !is.numeric(x0)) {
    .abort_bad_arg("x0", "must be numeric", quoted = FALSE)
  }
  if (isTRUE(strict) && !is.numeric(inv_s)) {
    .abort_bad_arg("inv_s", "must be numeric", quoted = FALSE)
  }
  x0_vec <- as.numeric(x0)
  inv_s_vec <- as.numeric(inv_s)
  if (length(x0_vec) != p || length(inv_s_vec) != p) {
    .abort_bad_arg("x0/inv_s", sprintf("must both have length ncol(%s)", x_arg), quoted = FALSE)
  }
  if (isTRUE(finite)) {
    if (.has_na_only(x0_vec) || .has_na_only(inv_s_vec)) {
      .abort_bad_arg("x0/inv_s", "cannot contain `NA` values", quoted = FALSE)
    }
    if (.has_nonfinite_non_na(x0_vec) || .has_nonfinite_non_na(inv_s_vec)) {
      .abort_bad_arg("x0/inv_s", "cannot contain non-finite values (`Inf`, `-Inf`, or `NaN`)", quoted = FALSE)
    }
  }
  list(x0 = x0_vec, inv_s = inv_s_vec)
}

.validate_vector_nrow <- function(x, n, arg = "x", ref = "data", strict = FALSE, finite = FALSE) {
  if (isTRUE(strict) && !is.numeric(x)) {
    .abort_bad_arg(arg, "must be numeric", quoted = FALSE)
  }
  x_vec <- as.numeric(x)
  if (length(x_vec) != n) {
    .abort_bad_length(arg, sprintf("equal nrow(%s)", ref), quoted = FALSE)
  }
  if (isTRUE(finite)) {
    .abort_if_na_or_nonfinite(x_vec, arg = arg, quoted = FALSE)
  }
  x_vec
}

.validate_same_n <- function(..., n, args = NULL, ref = "n") {
  vals <- list(...)
  if (is.null(args)) {
    args <- names(vals)
  }
  if (is.null(args) || any(!nzchar(args))) {
    args <- paste0("arg", seq_along(vals))
  }
  if (length(args) != length(vals)) {
    .abort_bad_arg("args", "length must match number of vectors")
  }
  for (i in seq_along(vals)) {
    if (length(vals[[i]]) != n) {
      .abort_bad_arg(args[[i]], sprintf("length must equal %s (%d)", ref, n), quoted = FALSE)
    }
  }
  invisible(TRUE)
}

.validate_subset_indices <- function(idxs_or_subs, n, s, unique_within_subset = TRUE, allow_empty = FALSE) {
  n <- .validate_scalar_int(n, arg = "n", min = 1L)
  s <- .validate_scalar_int(s, arg = "s", min = 1L)

  if (!is.matrix(idxs_or_subs)) {
    if (!is.numeric(idxs_or_subs)) {
      .abort_bad_arg("idxs", "must be a numeric or integer vector of indices")
    }
    idxs_num <- as.numeric(idxs_or_subs)
    if (length(idxs_num) != s) {
      .abort_bad_arg("idxs", sprintf("must have %d elements (order s)", s))
    }
    .abort_if_na_or_nonfinite(idxs_num, arg = "idxs")
    if (any(idxs_num != as.integer(idxs_num))) {
      .abort_bad_arg("idxs", "must contain integer indices")
    }
    idxs <- as.integer(idxs_num)
    if (length(idxs) != s) {
      .abort_bad_arg("idxs", sprintf("must have %d elements (order s)", s))
    }
    if (any(idxs < 1L | idxs > n)) {
      .abort_bad_arg("idxs", "must lie in 1..n")
    }
    if (isTRUE(unique_within_subset) && anyDuplicated(idxs)) {
      .abort_bad_arg("idxs", "must not contain duplicated indices")
    }
    return(matrix(idxs, nrow = 1L))
  }

  subs <- idxs_or_subs
  if (ncol(subs) != s) {
    .abort_bad_arg("subs", sprintf("must have %d columns (order s)", s))
  }
  if (nrow(subs) == 0L && !isTRUE(allow_empty)) {
    .abort_bad_arg("subs", "must have at least one row")
  }
  if (!is.numeric(subs)) {
    .abort_bad_arg("subs", "must be a numeric or integer matrix of indices")
  }
  .abort_if_na_or_nonfinite(subs, arg = "subs")
  if (any(subs != as.integer(subs))) {
    .abort_bad_arg("subs", "must contain integer indices")
  }
  if (!is.integer(subs)) {
    storage.mode(subs) <- "integer"
  }
  if (any(subs < 1L | subs > n)) {
    .abort_bad_arg("subs", "indices must lie in 1..n")
  }
  if (isTRUE(unique_within_subset) && any(apply(subs, 1L, function(row) anyDuplicated(row) > 0L))) {
    .abort_bad_arg("subs", "each row must not contain duplicated indices")
  }
  subs
}

.validate_permutation_vector <- function(ranks, arg = "ranks") {
  if (!is.numeric(ranks)) {
    .abort_bad_arg(arg, "must be numeric or integer ranks")
  }
  ranks_num <- as.numeric(ranks)
  .abort_if_na_or_nonfinite(ranks_num, arg = arg)
  if (any(ranks_num != as.integer(ranks_num))) {
    .abort_bad_arg(arg, "must contain integer ranks")
  }
  ranks <- as.integer(ranks_num)
  n <- length(ranks)
  if (n == 0L) .abort_bad_arg(arg, "must have positive length")
  if (!setequal(ranks, seq_len(n))) {
    .abort_bad_arg(arg, "must be a permutation of 1..n")
  }
  ranks
}

.validate_permutation_matrix <- function(rank_matrix, arg = "R", col_arg = "ranks") {
  rank_matrix <- as.matrix(rank_matrix)
  if (nrow(rank_matrix) == 0L || ncol(rank_matrix) == 0L) {
    .abort_bad_arg(arg, "must have positive nrow and ncol")
  }
  for (j in seq_len(ncol(rank_matrix))) {
    .validate_permutation_vector(rank_matrix[, j], arg = col_arg)
  }
  if (!is.integer(rank_matrix)) storage.mode(rank_matrix) <- "integer"
  rank_matrix
}

.warn_unknown_return <- function(return, valid) {
  unknown <- setdiff(return, valid)
  if (length(unknown) > 0L) {
    warning(sprintf(
      "Unknown `return` value(s) ignored: %s.\nValid options: %s.",
      paste(unknown, collapse = ", "),
      paste(valid, collapse = ", ")
    ), call. = FALSE)
  }
  invisible(NULL)
}

.validate_score_idxs <- function(idxs, n, s) {
  .validate_subset_indices(idxs, n, s, unique_within_subset = FALSE, allow_empty = FALSE)[1L, ]
}

.validate_score_subs <- function(subs, n, s) {
  .validate_subset_indices(subs, n, s, unique_within_subset = FALSE, allow_empty = FALSE)
}
