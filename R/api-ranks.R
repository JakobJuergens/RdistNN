# Lightweight rank-weight API + tiny registry (non-S3), consistent with distance API

.rank_registry <- new.env(parent = emptyenv())

# ---- predicates & adapters ---------------------------------------------------

#' Is this a rdistnn_rank_weight spec?
#'
#' A rank-weight spec is a list with `eval(ranks, ctx)` and/or
#' `eval_batch(rank_matrix, ctx)`, where `rank_matrix` is an n x m rank matrix
#' (each column a permutation of 1..n). Specs may also define `table(n)` as an
#' optional rank-to-weight lookup table and `setup(data)` for manual or
#' wrapper-level context construction.
#'
#' @param x Any object.
#' @return Logical.
#' @examples
#' is_rank_weight(rank_weight_dnn_exact(k = 1, s = 1))
#' @export
is_rank_weight <- function(x) {
  inherits(x, "rdistnn_rank_weight") ||
    (is.list(x) && identical(x$.__kind, "rank_weight"))
}

#' Construct a rank-weight spec
#'
#' @param name Character scalar.
#' @param eval Optional function `(ranks, ctx) -> numeric(n)`.
#' @param eval_batch Optional function `(rank_matrix, ctx) -> numeric(n x m)`.
#' @param setup Optional function to build context. At the top level,
#'   [weights_matrix()] passes any supplied `ctx` through to the rank-weight
#'   evaluator, but it does not automatically call `setup()` for you.
#' @return A validated `rdistnn_rank_weight` specification.
#' @examples
#' rw <- new_rank_weight(
#'   name = "inv_rank",
#'   eval = function(ranks, ctx = NULL) {
#'     w <- 1 / as.numeric(ranks)
#'     w / sum(w)
#'   }
#' )
#' is_rank_weight(rw)
#' @export
new_rank_weight <- function(name, eval = NULL, eval_batch = NULL, setup = NULL) {
  .rank_validate_spec(list(
    name = name,
    eval = eval,
    eval_batch = eval_batch,
    setup = setup
  ))
}

# Internal shape checks for rank-weight specs and outputs
.rank_validate_spec <- function(spec) {
  if (!is.list(spec)) .abort_bad_arg("spec", "must be a list")
  if (!is.null(spec$name) &&
      (!is.character(spec$name) || length(spec$name) != 1L ||
       is.na(spec$name) || !nzchar(spec$name))) {
    .abort_bad_arg("name", "must be a non-empty character scalar", quoted = FALSE)
  }
  if (!is.null(spec$eval) && !is.function(spec$eval)) {
    .abort_bad_arg("eval", "must be NULL or a function")
  }
  if (!is.null(spec$eval_batch) && !is.function(spec$eval_batch)) {
    .abort_bad_arg("eval_batch", "must be NULL or a function")
  }
  if (!is.null(spec$setup) && !is.function(spec$setup)) {
    .abort_bad_arg("setup", "must be NULL or a function")
  }
  if (!is.null(spec$table) && !is.function(spec$table)) {
    .abort_bad_arg("spec$table", "must be a function", quoted = FALSE)
  }
  has_eval <- !is.null(spec$eval) && is.function(spec$eval)
  has_eval_batch <- !is.null(spec$eval_batch) && is.function(spec$eval_batch)
  has_table <- !is.null(spec$table) && is.function(spec$table)
  if (!(has_eval || has_eval_batch || has_table)) {
    .abort_bad_arg("spec", "must provide at least one of `eval` or `eval_batch`, or `table`")
  }
  if (is.null(spec$name) || !nzchar(as.character(spec$name))) {
    spec$name <- "unnamed_rank_weight"
  }
  spec$.__kind <- "rank_weight"
  class(spec) <- union("rdistnn_rank_weight", class(spec))
  spec
}
.rank_validate_output_vec <- function(w, n) {
  if (!is.numeric(w) || length(w) != n || any(!is.finite(w))) {
    .dnn_abort("callback_invalid",
      "rank_weight$eval() must return finite numeric length(ranks).",
      argument = "rank_weight$eval()",
      expected = sprintf("finite numeric vector of length %d", n),
      actual = sprintf("%s of length %d", class(w)[1L], length(w)),
      hint = "Correct the callback to return the documented result."
    )
  }
  as.numeric(w)
}

.rank_weight_table <- function(spec, n) {
  if (!is_rank_weight(spec)) .abort_bad_arg("spec", "must be a rdistnn_rank_weight")
  n <- .validate_scalar_int(n, "n")
  if (is.null(spec$table)) {
    return(NULL)
  }
  if (!is.function(spec$table)) {
    .dnn_abort("argument_invalid", "`rank_weight$table` must be a function.", argument = "rank_weight$table")
  }

  tab <- spec$table(n)
  if (!is.numeric(tab) || length(tab) != n || any(!is.finite(tab))) {
    actual <- if (is.numeric(tab)) {
      sprintf("numeric vector of length %d", length(tab))
    } else {
      sprintf("%s of length %d", class(tab)[1L], length(tab))
    }
    .dnn_abort("callback_invalid", sprintf(
      "`rank_weight$table(n)` must return a finite numeric vector of length %d; got %s.",
      n, actual
    ), argument = "rank_weight$table", hint = "Correct the callback to return the documented finite rank-weight table.")
  }

  tab <- as.numeric(tab)
  tab_sum <- .compensated_sum(tab)
  if (!is.finite(tab_sum)) {
    # Bound every normalized partial sum without dividing tiny terms into
    # underflow. Keep those terms in original units to retain their sign when
    # the larger signed terms cancel. Only mass validation uses this fallback.
    scale <- 2^ceiling(log2(length(tab)))
    if (scale < length(tab)) scale <- scale * 2
    large <- abs(tab) >= .Machine$double.xmin * scale
    tab_sum <- .compensated_sum(c(
      .compensated_sum(tab[large] / scale) * scale,
      .compensated_sum(tab[!large])))
  }
  if (!is.finite(tab_sum) || tab_sum <= 0) {
    .dnn_abort("callback_invalid",
      "`rank_weight$table(n)` must have positive finite sum.",
      argument = "rank_weight$table",
      hint = "Correct the callback to return the documented finite rank-weight table."
    )
  }
  tab
}

# -----------------------------------------------------------------------------
.as_fun_accepts_ctx <- function(fun) {
  fmls <- formals(fun)
  if (is.null(fmls)) {
    return(FALSE)
  }
  # Accept either explicitly named 'ctx' or any >= 2-arg signature
  "ctx" %in% names(fmls) || length(fmls) >= 2L
}

#' Wrap a function as a rdistnn_rank_weight spec
#'
#' @param fun Function `(ranks, ctx) -> numeric(n)` OR just `(ranks) -> numeric(n)`.
#' @param name Character name.
#' @return A rank-weight spec (list) with fields `name`, `eval`.
#' @examples
#' rw <- as_rank_weight(function(ranks) {
#'   w <- 1 / as.numeric(ranks)
#'   w / sum(w)
#' })
#' rank_weight_eval(rw, c(3L, 1L, 2L))
#' @export
as_rank_weight <- function(fun, name = "user_rank_weight") {
  if (!is.function(fun)) {
    .abort_bad_arg("fun", "must be a function")
  }
  takes_ctx <- .as_fun_accepts_ctx(fun)

  new_rank_weight(
    name = name,
    eval = function(ranks, ctx = NULL) {
      r <- .rank_check_vec(ranks)
      w <- if (takes_ctx) fun(r, ctx) else fun(r)
      .rank_validate_output_vec(w, length(r))
    }
  )
}

# ---- registry ---------------------------------------------------------------

#' Register a rdistnn_rank_weight under a name
#'
#' If `spec` is a function, it is wrapped with [as_rank_weight()].
#'
#' @param name Character scalar registry key.
#' @param spec A rdistnn_rank_weight spec or function `(ranks, ctx) -> numeric(n)`.
#' @return Invisibly, the registered spec.
#' @examples
#' nm <- "tmp_rank_demo"
#' register_rank_weight(nm, function(ranks) rep(1 / length(ranks), length(ranks)))
#' nm %in% rank_weight_names()
#' @export
register_rank_weight <- function(name, spec) {
  name <- .validate_registry_name(name)
  if (is.function(spec)) spec <- as_rank_weight(spec, name = name)
  if (is.list(spec)) {
    missing_name <- is.null(spec$name) ||
      (is.character(spec$name) && length(spec$name) == 1L &&
       (is.na(spec$name) || !nzchar(spec$name)))
    if (missing_name) {
      spec$name <- name
    }
  }
  if (is.list(spec)) spec <- .rank_validate_spec(spec)
  if (!is_rank_weight(spec)) .abort_bad_arg("spec", "is not a rdistnn_rank_weight")
  assign(name, spec, envir = .rank_registry)
  invisible(spec)
}

#' Retrieve a registered rdistnn_rank_weight
#' @param name Character scalar.
#' @return A rdistnn_rank_weight spec.
#' @examples
#' register_rank_weight("tmp_rank_get", function(ranks) rep(1 / length(ranks), length(ranks)))
#' is_rank_weight(get_rank_weight("tmp_rank_get"))
#' @export
get_rank_weight <- function(name) {
  name <- .validate_registry_name(name)
  if (!exists(name, envir = .rank_registry, inherits = FALSE)) {
    known <- .registry_names(.rank_registry)
    .abort_bad_arg(
      "name",
      "Unknown rdistnn_rank_weight",
      expected = if (length(known) > 0L) paste(known, collapse = ", ") else "<none registered>",
      actual = as.character(name)
    )
  }
  get(name, envir = .rank_registry, inherits = FALSE)
}

#' List registered rdistnn_rank_weight names
#' @return Character vector.
#' @examples
#' head(rank_weight_names(), 3)
#' @export
rank_weight_names <- function() {
  .registry_names(.rank_registry)
}

# ---- evaluators -------------------------------------------------------------

# internal validators used by the evaluators
.rank_check_vec <- function(ranks) {
  .validate_permutation_vector(ranks, arg = "ranks")
}
.rank_check_mat <- function(rank_matrix) {
  .validate_permutation_matrix(rank_matrix, arg = "R", col_arg = "ranks")
}

#' Evaluate rank-weight for a single target
#'
#' @param spec Rank-weight spec (see [is_rank_weight()]).
#' @param ranks Integer permutation of `1..n` (length `n`).
#' @param ctx Optional context (if the spec uses it).
#' @return Numeric vector of length `n`.
#' @examples
#' rw <- as_rank_weight(function(r) rep(1 / length(r), length(r)))
#' rank_weight_eval(rw, c(2L, 1L, 3L))
#' @export
rank_weight_eval <- function(spec, ranks, ctx = NULL) {
  if (!is_rank_weight(spec)) .abort_bad_arg("spec", "must be a rdistnn_rank_weight")
  r <- .rank_check_vec(ranks)

  if (is.function(spec$table) || is.function(spec$eval_batch)) {
    # Prefer batch path whenever available.
    rank_matrix_one <- matrix(r, nrow = length(r), ncol = 1L)
    weight_matrix_one <- rank_weight_eval_batch(spec, rank_matrix_one, ctx)
    w <- weight_matrix_one[, 1L]
  } else if (is.function(spec$eval)) {
    w <- spec$eval(r, ctx)
  } else {
    .abort_bad_arg("spec", "rdistnn_rank_weight spec has neither eval nor eval_batch nor table")
  }

  if (!is.numeric(w) || length(w) != length(r) || any(!is.finite(w))) {
    .dnn_abort("callback_invalid",
      "rank_weight$eval() must return finite numeric of length(ranks).",
      argument = "rank_weight$eval()",
      expected = sprintf("finite numeric vector of length %d", length(r)),
      actual = sprintf("%s of length %d", class(w)[1L], length(w)),
      hint = "Correct the callback to return the documented result."
    )
  }
  as.numeric(w)
}

#' Evaluate rank-weight for many targets (batch)
#'
#' A scalar-only rule is evaluated and validated separately for each target.
#' The assembled result remains an `n x m` matrix when `n = 1` or `m = 1`.
#' A supplied `eval_batch` must itself return that matrix shape.
#'
#' @param spec Rank-weight spec (see [is_rank_weight()]).
#' @param R Integer matrix (n x m), each column a permutation of `1..n`.
#' @param ctx Optional context (if the spec uses it).
#' @return Numeric matrix (n x m).
#' @examples
#' rw <- as_rank_weight(function(r) rep(1 / length(r), length(r)))
#' rank_matrix <- cbind(c(1L, 2L, 3L), c(2L, 1L, 3L))
#' rank_weight_eval_batch(rw, rank_matrix)
#' dim(rank_weight_eval_batch(rw, matrix(1L, nrow = 1L, ncol = 2L)))
#' @export
rank_weight_eval_batch <- function(spec, R, ctx = NULL) {
  if (!is_rank_weight(spec)) .abort_bad_arg("spec", "must be a rdistnn_rank_weight")
  rank_matrix <- .rank_check_mat(R)

  tab <- .rank_weight_table(spec, nrow(rank_matrix))
  if (!is.null(tab)) {
    weight_matrix <- matrix(tab[rank_matrix], nrow = nrow(rank_matrix), ncol = ncol(rank_matrix))
  } else if (is.function(spec$eval_batch)) {
    weight_matrix <- spec$eval_batch(rank_matrix, ctx)
  } else if (is.function(spec$eval)) {
    # Validate each scalar-target result before assembling the batch.
    weight_matrix <- apply(rank_matrix, 2L, function(r) {
      .rank_validate_output_vec(spec$eval(r, ctx), length(r))
    })
    if (!is.matrix(weight_matrix)) {
      weight_matrix <- matrix(weight_matrix, nrow = nrow(rank_matrix), ncol = ncol(rank_matrix))
    }
  } else {
    .abort_bad_arg("spec", "rdistnn_rank_weight spec has neither eval nor eval_batch nor table")
  }

  if (!is.matrix(weight_matrix) || !is.numeric(weight_matrix) ||
    any(!is.finite(weight_matrix)) ||
    nrow(weight_matrix) != nrow(rank_matrix) || ncol(weight_matrix) != ncol(rank_matrix)) {
    got_dim <- if (is.matrix(weight_matrix)) {
      sprintf("%d x %d", nrow(weight_matrix), ncol(weight_matrix))
    } else {
      sprintf("non-matrix (%s)", class(weight_matrix)[1L])
    }
    .dnn_abort(
      "callback_invalid",
      sprintf(
      "rank_weight$eval_batch must return a finite numeric matrix of shape %d x %d; got %s%s.",
      nrow(rank_matrix), ncol(rank_matrix), got_dim,
      if (is.matrix(weight_matrix) && nrow(weight_matrix) == nrow(rank_matrix) &&
          ncol(weight_matrix) == ncol(rank_matrix)) " (contains non-finite values)" else ""
    ),
      argument = "rank_weight",
      hint = "Correct the callback to return the documented type, dimensions and finite values."
    )
  }
  weight_matrix
}
