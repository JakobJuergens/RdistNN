# Fit input routes, column identities and canonical target matrices.

.dnn_validate_data_column_names <- function(data, arg = "data") {
  col_names <- names(data)
  if (is.null(col_names) || length(col_names) != ncol(data) ||
      anyNA(col_names) || any(!nzchar(col_names))) {
    .abort_bad_arg(arg, "must have non-missing, non-empty column names")
  }

  duplicated_names <- unique(col_names[duplicated(col_names)])
  if (length(duplicated_names) > 0L) {
    .abort_bad_arg(
      arg,
      "must have unique column names",
      actual = sprintf(
        "duplicate column name(s): %s",
        paste(sprintf("'%s'", duplicated_names), collapse = ", ")
      )
    )
  }

  col_names
}

.dnn_resolve_data_cols <- function(cols, data, arg, allow_multiple = TRUE) {
  if (is.null(cols)) {
    .abort_bad_arg(arg, "must be provided")
  }

  data_names <- .dnn_validate_data_column_names(data)

  if (is.character(cols)) {
    if (length(cols) == 0L || anyNA(cols) || any(!nzchar(cols))) {
      .abort_bad_arg(arg, "must be a non-empty character vector of column names")
    }
    missing <- setdiff(cols, data_names)
    if (length(missing) > 0L) {
      .dnn_abort("argument_invalid", sprintf(
        "Column(s) %s (specified by `%s`) not found in `data`. Available columns: %s.",
        .dnn_message_examples(sprintf("'%s'", missing), width = 40L)$text,
        .dnn_brief(arg),
        .dnn_message_examples(data_names, width = 40L)$text
      ), argument = arg,
        expected = .dnn_message_examples(data_names)$text,
        actual = .dnn_message_examples(missing)$text,
        hint = "Check the selector against the available column names.")
    }
    resolved <- cols
  } else if (is.numeric(cols)) {
    if (length(cols) == 0L) {
      .abort_bad_arg(arg, "must be a non-empty integer vector of column indices")
    }
    if (.has_na_only(cols)) {
      .abort_bad_arg(arg, "cannot contain `NA` indices")
    }
    if (.has_nonfinite_non_na(cols)) {
      .abort_bad_arg(arg, "cannot contain non-finite indices (`Inf`, `-Inf`, or `NaN`)")
    }
    if (any(cols != trunc(cols))) {
      .abort_bad_arg(arg, "must be a non-empty integer vector of column indices")
    }
    p <- ncol(data)
    if (any(cols < 1L | cols > p | cols > .Machine$integer.max)) {
      .abort_bad_arg(arg, sprintf("contains out-of-range indices; valid range is 1..%d", p))
    }
    idx <- as.integer(cols)
    resolved <- data_names[idx]
  } else {
    .abort_bad_arg(arg, "must be a character vector of column names or an integer vector of indices")
  }

  if (!isTRUE(allow_multiple) && length(resolved) != 1L) {
    .abort_bad_arg(arg, "must reference exactly one column")
  }
  resolved
}

# Extract and validate a single data-frame column by name or index.
.dnn_col_vec <- function(col, data, n, arg, ref = "data", finite = FALSE) {
  col_name <- .dnn_resolve_data_cols(col, data = data, arg = arg, allow_multiple = FALSE)
  .validate_vector_nrow(data[[col_name]], n, arg = arg, ref = ref, strict = TRUE, finite = finite)
}

# Extract and validate multiple data-frame columns by names/indices.
.dnn_col_mat <- function(cols, data, arg, out_arg = "X") {
  col_names <- .dnn_resolve_data_cols(cols, data = data, arg = arg, allow_multiple = TRUE)
  .validate_matrix_nonempty(data[, col_names, drop = FALSE], arg = out_arg)
}

.dnn_resolve_targets <- function(x_mat, X0 = NULL, x0 = NULL) {
  if (!is.null(X0) && !is.null(x0)) {
    .abort_bad_arg("X0/x0", "must provide either `X0` or `x0`, not both", quoted = FALSE)
  }
  if (is.null(X0)) {
    if (is.null(x0)) {
      .abort_bad_arg("X0/x0", "must provide `X0` (m x p) or a single `x0` (length p)", quoted = FALSE)
    }
    x0_vec <- .validate_target_vector(
      x0,
      ncol(x_mat),
      x0_arg = "x0",
      x_arg = "X",
      strict = TRUE,
      finite = TRUE
    )
    return(matrix(x0_vec, nrow = 1L))
  }
  x0_mat <- .validate_matrix(
    X0,
    arg = "X0",
    min_rows = 1L,
    min_cols = 1L,
    numeric = TRUE,
    finite = TRUE
  )
  .validate_ncol_relation(x_mat, x0_mat, y_arg = "X0", x_arg = "X", relation = "equal")
}

.dnn_x0_with_names <- function(x0_mat, x_mat) {
  # Result objects store the evaluation points so downstream consumers
  # (plotting, reporting) can reconstruct target profiles. Propagate
  # covariate names from the training matrix when the targets have none.
  if (is.null(colnames(x0_mat)) && !is.null(colnames(x_mat)) &&
      ncol(x0_mat) == ncol(x_mat)) {
    colnames(x0_mat) <- colnames(x_mat)
  }
  x0_mat
}

# Reject ignored input-route arguments before canonicalization. Allowed direct
# signals (data + y or data + pseudo_outcome) are deliberately excluded.
.dnn_validate_input_route <- function(data, direct, selectors) {
  supplied <- function(x) names(x)[!vapply(x, is.null, logical(1L))]
  if (!is.null(data)) {
    bad <- supplied(direct)
    if (length(bad)) .abort_bad_arg(bad[1L], "cannot be combined with `data`; use one input interface")
  } else if (length(supplied(selectors))) {
    .abort_bad_arg("data", "must be supplied when column selectors are used")
  }
  invisible(NULL)
}
