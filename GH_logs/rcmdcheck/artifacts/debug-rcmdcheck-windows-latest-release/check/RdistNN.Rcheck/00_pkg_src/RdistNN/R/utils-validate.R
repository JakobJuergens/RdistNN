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

.validate_scalar_int <- function(x, arg = "x", min = 1L) {
  if (!is.numeric(x) || length(x) != 1L || !is.finite(x) || x != as.integer(x)) {
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

.validate_binary01 <- function(x, arg = "x") {
  x_vec <- as.numeric(x)
  if (anyNA(x_vec) || any(!is.finite(x_vec)) || any(!(x_vec %in% c(0, 1)))) {
    .abort_bad_arg(arg, "must be binary in {0,1}")
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
  if (isTRUE(finite) && any(!is.finite(x_mat))) {
    .abort_bad_arg(arg, "must be finite", quoted = FALSE)
  }
  x_mat
}

.validate_matrix_nonempty <- function(x, arg = "X") {
  x_mat <- .validate_matrix(x, arg = arg, min_rows = 0L, min_cols = 1L, numeric = FALSE)
  if (ncol(x_mat) == 0L) {
    .abort_bad_arg(arg, "must have at least one column", quoted = FALSE)
  }
  x_mat
}

.validate_ncol_relation <- function(x_mat, y, y_arg = "X0", x_arg = "X", relation = c("match", "equal")) {
  relation <- match.arg(relation)
  y_mat <- as.matrix(y)
  if (ncol(x_mat) != ncol(y_mat)) {
    .abort_bad_state(sprintf("ncol(%s) must %s ncol(%s).", y_arg, relation, x_arg))
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
  if (isTRUE(finite) && any(!is.finite(x0_vec))) {
    .abort_bad_arg(x0_arg, "must be finite", quoted = FALSE)
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
    .abort_bad_state(sprintf("x0 and inv_s must be length p == ncol(%s).", x_arg))
  }
  if (isTRUE(finite) && (any(!is.finite(x0_vec)) || any(!is.finite(inv_s_vec)))) {
    .abort_bad_state("x0 and inv_s must be finite.")
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
  if (isTRUE(finite) && any(!is.finite(x_vec))) {
    .abort_bad_arg(arg, "must be finite", quoted = FALSE)
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
    .abort_bad_state("`args` length must match number of vectors.")
  }
  for (i in seq_along(vals)) {
    if (length(vals[[i]]) != n) {
      .abort_bad_state(sprintf("%s length must equal %s (%d).", args[[i]], ref, n))
    }
  }
  invisible(TRUE)
}

.validate_subset_indices <- function(idxs_or_subs, n, s, unique_within_subset = TRUE, allow_empty = FALSE) {
  n <- .validate_scalar_int(n, arg = "n", min = 1L)
  s <- .validate_scalar_int(s, arg = "s", min = 1L)

  if (!is.matrix(idxs_or_subs)) {
    idxs <- as.integer(idxs_or_subs)
    if (length(idxs) != s) {
      .abort_bad_state(sprintf("`idxs` must have %d elements (order s).", s))
    }
    if (anyNA(idxs)) {
      .abort_bad_state("`idxs` cannot contain NA.")
    }
    if (any(idxs < 1L | idxs > n)) {
      .abort_bad_state("`idxs` must lie in 1..n.")
    }
    if (isTRUE(unique_within_subset) && anyDuplicated(idxs)) {
      .abort_bad_state("`idxs` must not contain duplicated indices.")
    }
    return(matrix(idxs, nrow = 1L))
  }

  subs <- idxs_or_subs
  if (ncol(subs) != s) {
    .abort_bad_state(sprintf("`subs` must have %d columns (order s).", s))
  }
  if (nrow(subs) == 0L && !isTRUE(allow_empty)) {
    .abort_bad_state("`subs` must have at least one row.")
  }
  if (!is.integer(subs)) {
    storage.mode(subs) <- "integer"
  }
  if (any(!is.finite(subs))) {
    .abort_bad_state("`subs` contains non-finite entries.")
  }
  if (any(subs < 1L | subs > n)) {
    .abort_bad_state("`subs` indices must lie in 1..n.")
  }
  if (isTRUE(unique_within_subset) && any(apply(subs, 1L, function(row) anyDuplicated(row) > 0L))) {
    .abort_bad_state("Each row of `subs` must not contain duplicated indices.")
  }
  subs
}

.validate_permutation_vector <- function(ranks, arg = "ranks") {
  ranks <- as.integer(ranks)
  n <- length(ranks)
  if (n == 0L) .abort_bad_arg(arg, "must have positive length")
  if (anyNA(ranks)) .abort_bad_arg(arg, "cannot contain NA")
  if (!setequal(ranks, seq_len(n))) {
    .abort_bad_arg(arg, "must be a permutation of 1..n")
  }
  ranks
}

.validate_permutation_matrix <- function(rank_matrix, arg = "R", col_arg = "ranks") {
  rank_matrix <- as.matrix(rank_matrix)
  if (!is.integer(rank_matrix)) storage.mode(rank_matrix) <- "integer"
  if (nrow(rank_matrix) == 0L || ncol(rank_matrix) == 0L) {
    .abort_bad_arg(arg, "must have positive nrow and ncol")
  }
  for (j in seq_len(ncol(rank_matrix))) {
    .validate_permutation_vector(rank_matrix[, j], arg = col_arg)
  }
  rank_matrix
}

.validate_score_idxs <- function(idxs, n, s) {
  .validate_subset_indices(idxs, n, s, unique_within_subset = FALSE, allow_empty = FALSE)[1L, ]
}

.validate_score_subs <- function(subs, n, s) {
  .validate_subset_indices(subs, n, s, unique_within_subset = FALSE, allow_empty = FALSE)
}
