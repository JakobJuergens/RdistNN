# Only bounded descriptions belong in conditions. Never retain a failing array,
# fitted object or call environment merely to format a diagnostic.
.dnn_brief <- function(x, limit = 80L) {
  x <- encodeString(as.character(x), quote = "")
  long <- nchar(x, type = "chars") > limit
  x[long] <- paste0(substr(x[long], 1L, limit - 3L), "...")
  x
}

.dnn_message_examples <- function(examples, width = 80L) {
  if (is.null(examples)) return(list(values = NULL, text = NULL))
  values <- .dnn_brief(utils::head(examples, 5L), width)
  omitted <- length(examples) - length(values)
  text <- paste(values, collapse = ", ")
  if (omitted > 0L) text <- sprintf("%s (%d more)", text, omitted)
  list(values = values, text = text)
}

# The dynamic scope is installed at a public orchestration boundary and dies
# with that operation. Numerical loops do not look it up or signal conditions.
.dnn_reporting_scope <- function() {
  for (frame in rev(sys.frames())) {
    if (exists("..dnn_distance_warning_scope", envir = frame, inherits = FALSE)) {
      scope <- get("..dnn_distance_warning_scope", envir = frame, inherits = FALSE)
      if (is.environment(scope)) return(scope)
    }
  }
  NULL
}

# Resolve a public function only while signaling. Retain its bounded name,
# never the call's arguments, closure, frame or traceback.
.dnn_public_operation <- function() {
  ns <- environment(.dnn_abort)
  if (!isNamespace(ns)) return(NULL)
  public <- getNamespaceExports(ns)
  methods <- getNamespaceInfo(ns, "S3methods")
  if (length(methods)) public <- unique(c(public, methods[, 3L]))
  calls <- sys.calls()
  for (i in seq_along(calls)) {
    fun <- sys.function(i)
    if (is.primitive(fun) || !identical(environment(fun), ns)) next
    head <- calls[[i]][[1L]]
    name <- if (is.symbol(head)) as.character(head) else if (is.call(head) &&
      as.character(head[[1L]]) %in% c("::", ":::")) as.character(head[[3L]]) else ""
    if (name %in% public && identical(fun, get(name, ns, inherits = FALSE))) return(.dnn_brief(name))
    # A user can bind a public function to a local alias or call it via do.call.
    for (candidate in public) {
      if (exists(candidate, ns, inherits = FALSE) &&
          identical(fun, get(candidate, ns, inherits = FALSE))) return(.dnn_brief(candidate))
    }
  }
  NULL
}

.dnn_condition_message <- function(message, operation) {
  if (is.null(operation) || startsWith(message, paste0(operation, "():"))) return(message)
  paste0(operation, "(): ", message)
}

.dnn_abort <- function(code, message, argument = NULL, expected = NULL, actual = NULL,
                       hint = NULL, parent = NULL) {
  families <- c(argument_invalid = "rdistnn_input_error", state_invalid = "rdistnn_state_error",
    extra_unavailable = "rdistnn_extra_unavailable", dependency_missing = "rdistnn_dependency_error",
    engine_unavailable = "rdistnn_engine_error", numerical_range = "rdistnn_numerical_error",
    candidate_failure = "rdistnn_candidate_failure", callback_invalid = "rdistnn_callback_error",
    internal_invariant = "rdistnn_internal_error")
  family <- unname(families[code])
  family <- family[!is.na(family)]
  operation <- .dnn_public_operation()
  hint <- if (is.null(hint)) NULL else .dnn_brief(hint, 240L)
  cnd <- structure(list(message = .dnn_condition_message(paste(c(message, hint), collapse = " "), operation),
    call = NULL, operation = operation,
    code = code, argument = if (is.null(argument)) NULL else .dnn_brief(argument),
    expected = if (is.null(expected)) NULL else .dnn_brief(expected),
    actual = if (is.null(actual)) NULL else .dnn_brief(actual), hint = hint,
    parent = if (is.null(parent)) NULL else list(message = .dnn_brief(conditionMessage(parent), 240L),
      class = .dnn_brief(utils::head(class(parent), 5L)))),
    class = c(family, "rdistnn_error", "error", "condition"))
  stop(cnd)
}

.dnn_warn <- function(code, message, affected = NULL, examples = NULL, hint = NULL, location = NULL) {
  scope <- .dnn_reporting_scope()
  if (!is.null(scope)) {
    if (code %in% scope$emitted) return(invisible(NULL))
    scope$emitted <- c(scope$emitted, code)
  }
  details <- .dnn_message_examples(examples)
  if (is.null(location)) {
    location <- switch(code, partial_result = "target", weak_localization = "target",
      infeasible_candidates = "candidate", NULL)
  }
  locations <- if (is.null(location)) {
    details$text
  } else {
    labeled <- .dnn_message_examples(paste(location, utils::head(examples, 5L)))$text
    omitted <- length(examples) - length(details$values)
    if (omitted > 0L) paste0(labeled, " (", omitted, " more)") else labeled
  }
  operation <- .dnn_public_operation()
  hint <- if (is.null(hint)) NULL else .dnn_brief(hint, 240L)
  cnd <- structure(list(message = .dnn_condition_message(paste(c(message,
    if (length(details$values)) paste0("First affected: ", locations, "."), hint), collapse = " "), operation),
    call = NULL, operation = operation, code = code, affected = affected,
    examples = details$values, location = location, hint = hint),
    class = c(paste0("rdistnn_", code, "_warning"), "rdistnn_warning", "warning", "condition"))
  warning(cnd)
  invisible(NULL)
}

.abort_bad_arg <- function(arg, msg, expected = NULL, actual = NULL, quoted = TRUE) {
  arg_lbl <- if (isTRUE(quoted)) sprintf("`%s`", .dnn_brief(arg)) else .dnn_brief(arg)
  parts <- c(sprintf("%s %s.", arg_lbl, msg))
  if (!is.null(expected)) parts <- c(parts, sprintf("Expected: %s.", .dnn_brief(expected)))
  if (!is.null(actual)) parts <- c(parts, sprintf("Got: %s.", .dnn_brief(actual)))
  .dnn_abort("argument_invalid", paste(parts, collapse = " "), argument = arg,
    expected = expected, actual = actual)
}

# Match public choices with base R's NULL/default/partial/multiple semantics,
# while identifying the actual argument. No broad error handler can relabel an
# interrupt or resource failure from matching or from a caller's promise.
.dnn_match_arg <- function(arg, choices, several.ok = FALSE,
                           argument = as.character(substitute(arg))) {
  if (missing(choices)) {
    caller <- sys.parent()
    choices <- eval(formals(sys.function(caller))[[argument]], envir = sys.frame(caller))
  }
  if (is.null(arg)) return(choices[1L])
  invalid <- function(problem) {
    valid <- .dnn_message_examples(choices)$text
    .abort_bad_arg(argument, problem, expected = paste0("one of: ", valid),
      actual = if (is.character(arg)) .dnn_message_examples(arg)$text else
        sprintf("%s of length %d", typeof(arg), length(arg)))
  }
  if (!is.character(arg)) invalid("must be NULL or a character vector")
  all_match <- is.character(several.ok) && length(several.ok) == 1L &&
    !is.na(several.ok) && startsWith(several.ok, "all")
  if (all_match) several.ok <- TRUE
  if (!is.logical(several.ok) || length(several.ok) != 1L || is.na(several.ok)) {
    .abort_bad_arg("several.ok", "must be TRUE, FALSE or a string starting with 'all'")
  }
  if (!several.ok) {
    if (identical(arg, choices)) return(arg[1L])
    if (length(arg) != 1L) invalid("must have length 1")
  } else if (!length(arg)) {
    invalid("must have length at least 1")
  }
  matched <- pmatch(arg, choices, nomatch = 0L, duplicates.ok = TRUE)
  if (all(matched == 0L) || (all_match && any(matched == 0L))) {
    invalid("must match an available choice")
  }
  choices[matched[matched > 0L]]
}

# Reuse the ordinary validators at a retained-state owner, while distinguishing
# corruption from a new bad argument. Only owned input validation is translated;
# interruptions, resource failures and unexpected errors keep their identity.
.dnn_validate_retained <- function(code) {
  tryCatch(code, rdistnn_input_error = function(e) {
    .dnn_abort("state_invalid", conditionMessage(e),
      argument = e$argument, expected = e$expected, actual = e$actual,
      hint = "Refit with the current estimator before reusing this fit.", parent = e)
  })
}

# Evaluate the callback first. This boundary translates only known validators
# applied to its returned value, never exceptions raised by user execution.
.dnn_validate_callback <- function(code, argument = NULL) {
  tryCatch(code, rdistnn_input_error = function(e) {
    .dnn_abort("callback_invalid", conditionMessage(e),
      argument = if (is.null(argument)) e$argument else argument,
      expected = e$expected, actual = e$actual,
      hint = "Correct the callback to return the documented result.", parent = e)
  })
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
      .abort_bad_arg(arg, "must be a positive integer", expected = "an integer >= 1", actual = x_int)
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
  if (any(x_mat != trunc(x_mat)) || any(x_mat > .Machine$integer.max)) {
    .abort_bad_arg(arg, "must contain integer indices")
  }
  if (any(x_mat < 1)) {
    .abort_bad_arg(arg, "indices must be >= 1")
  }
  x_int <- matrix(as.integer(x_mat), nrow = nrow(x_mat), ncol = ncol(x_mat))
  dimnames(x_int) <- dimnames(x_mat)
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
    if (any(idxs_num != trunc(idxs_num))) {
      .abort_bad_arg("idxs", "must contain integer indices")
    }
    if (any(idxs_num < 1 | idxs_num > n)) {
      .abort_bad_arg("idxs", "must lie in 1..n")
    }
    idxs <- as.integer(idxs_num)
    if (length(idxs) != s) {
      .abort_bad_arg("idxs", sprintf("must have %d elements (order s)", s))
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
  if (any(subs != trunc(subs))) {
    .abort_bad_arg("subs", "must contain integer indices")
  }
  if (any(subs < 1 | subs > n)) {
    .abort_bad_arg("subs", "indices must lie in 1..n")
  }
  if (!is.integer(subs)) {
    storage.mode(subs) <- "integer"
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
  if (any(ranks_num != trunc(ranks_num))) {
    .abort_bad_arg(arg, "must contain integer ranks")
  }
  n <- length(ranks_num)
  if (n == 0L) .abort_bad_arg(arg, "must have positive length")
  if (!setequal(ranks_num, seq_len(n))) {
    .abort_bad_arg(arg, "must be a permutation of 1..n")
  }
  as.integer(ranks_num)
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
  if (!is.null(return) && (!is.character(return) || !is.null(dim(return)) || anyNA(return))) {
    .abort_bad_arg("return", "must be NULL or a character vector without missing values",
      expected = paste("a subset of", paste(valid, collapse = ", ")))
  }
  unknown <- setdiff(return, valid)
  if (length(unknown) > 0L) {
    .dnn_warn("unknown_return", "Unknown `return` value(s) ignored.",
      affected = length(unknown), examples = unknown,
      hint = sprintf("Valid options: %s.", paste(valid, collapse = ", ")))
  }
  invisible(NULL)
}

.validate_score_idxs <- function(idxs, n, s) {
  .validate_subset_indices(idxs, n, s, unique_within_subset = FALSE, allow_empty = FALSE)[1L, ]
}

.validate_score_subs <- function(subs, n, s) {
  .validate_subset_indices(subs, n, s, unique_within_subset = FALSE, allow_empty = TRUE)
}
