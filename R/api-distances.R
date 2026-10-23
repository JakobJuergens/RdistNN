# -----------------------------------------------------------------------------
# Unified API + registry for distance specifications

# lightweight registry (name -> spec)
.dist_registry <- new.env(parent = emptyenv())

# --- type + adapters ----------------------------------------------------------
.distance_has_eval <- function(spec) !is.null(spec$eval) && is.function(spec$eval)
.distance_has_eval_batch <- function(spec) !is.null(spec$eval_batch) && is.function(spec$eval_batch)
.distance_eval_is_single <- function(fun) {
  fmls <- tryCatch(names(formals(fun)), error = function(e) character())
  length(fmls) > 0L && ("x0" %in% fmls) && !("X0" %in% fmls)
}
.distance_call_user <- function(fun, X, x0, ctx) {
  args <- list(X = X, x0 = x0, ctx = ctx)
  fmls <- tryCatch(names(formals(fun)), error = function(e) character())
  if ("..." %in% fmls) {
    return(do.call(fun, args))
  }
  do.call(fun, args[intersect(names(args), fmls)])
}
.distance_validate_output_vec <- function(d, n, src = "Distance function") {
  if (!is.numeric(d) || length(d) != n || any(!is.finite(d))) {
    stop(sprintf(
      "%s must return a finite numeric vector of length n=%d; got %s of length %d%s.",
      src, n, class(d)[1L], length(d),
      if (is.numeric(d) && length(d) == n) " (contains non-finite values)" else ""
    ), call. = FALSE)
  }
  as.numeric(d)
}
.distance_validate_output_mat <- function(d_mat, n, m, src = "Distance function") {
  if (!is.matrix(d_mat) || !is.numeric(d_mat) || nrow(d_mat) != n || ncol(d_mat) != m || any(!is.finite(d_mat))) {
    got_dim <- if (is.matrix(d_mat)) {
      sprintf("%d x %d", nrow(d_mat), ncol(d_mat))
    } else {
      sprintf("non-matrix (%s)", class(d_mat)[1L])
    }
    stop(sprintf(
      "%s must return a finite numeric matrix of shape %d x %d; got %s%s.",
      src, n, m, got_dim,
      if (is.matrix(d_mat) && is.numeric(d_mat) && nrow(d_mat) == n && ncol(d_mat) == m)
        " (contains non-finite values)" else ""
    ), call. = FALSE)
  }
  d_mat
}
.distance_validate_target_matrix <- function(X0, x_mat) {
  x0_mat <- .validate_matrix(
    X0,
    arg = "X0",
    min_rows = 0L,
    min_cols = 1L,
    numeric = TRUE,
    finite = TRUE
  )
  .validate_ncol_relation(x_mat, x0_mat, y_arg = "X0", x_arg = "X", relation = "match")
}
.distance_validate_spec <- function(spec) {
  if (!is.list(spec)) .abort_bad_arg("spec", "must be a list")
  if (!is.null(spec$name) &&
    (!is.character(spec$name) || length(spec$name) != 1L || is.na(spec$name) || !nzchar(spec$name))) {
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
  if (!(.distance_has_eval(spec) || .distance_has_eval_batch(spec))) {
    .abort_bad_arg("spec", "must provide at least one of `eval` or `eval_batch`")
  }
  if (is.null(spec$name) || !nzchar(as.character(spec$name))) {
    spec$name <- "unnamed_distance"
  }
  class(spec) <- union("dnn_distance", class(spec))
  spec
}

#' Is this a distance spec?
#'
#' A distance spec is a list with `eval(X, x0, ctx)` and/or
#' `eval_batch(X, X0, ctx)`. It may optionally define `setup(X)` to precompute
#' a context (e.g., scaling).
#'
#' @param x Any R object.
#' @return Logical.
#' @examples
#' is_distance(distance_l2())
#' is_distance(list(eval = function(X, x0, ctx) rep(0, nrow(X))))
#' @export
is_distance <- function(x) {
  inherits(x, "dnn_distance") ||
    (is.list(x) && (is.function(x$eval) || is.function(x$eval_batch)))
}

#' Construct a distance spec
#'
#' @param name Character scalar.
#' @param eval Optional function `(X, x0, ctx) -> numeric(n)`.
#' @param eval_batch Optional function `(X, X0, ctx) -> numeric(n x m)`.
#' @param setup Optional function `(X) -> ctx`.
#' @details A batch evaluator must return an actual `n x m` numeric matrix,
#' including when either dimension is one. With only `eval`, the shared batch
#' evaluator assembles those dimensions for you and validates each result.
#' When `ctx` is `NULL`, `setup(X)` runs once per evaluation call, not once
#' per target in a batch. [weights_matrix()] uses this same evaluator.
#' @return A validated `dnn_distance` spec.
#' @examples
#' spec <- new_distance(
#'   name = "sq_l2",
#'   eval = function(X, x0, ctx = NULL) rowSums((X - matrix(x0, nrow(X), ncol(X), TRUE))^2)
#' )
#' is_distance(spec)
#' @export
new_distance <- function(name, eval = NULL, eval_batch = NULL, setup = NULL) {
  .distance_validate_spec(list(
    name = name,
    eval = eval,
    eval_batch = eval_batch,
    setup = setup
  ))
}

#' Wrap a function as a distance spec
#'
#' @param fun Function `(X, x0, ctx) -> numeric(n)`.
#' @param name Character name for the spec.
#' @return A distance spec with class `dnn_distance`.
#' @examples
#' spec <- as_distance(
#'   function(X, x0, ctx = NULL) {
#'     rowSums(abs(X - matrix(x0, nrow(X), ncol(X), TRUE)))
#'   }
#' )
#' x_mat <- matrix(1:9, nrow = 3)
#' distance_eval(spec, x_mat, c(0, 0, 0))
#' @export
as_distance <- function(fun, name = "user_distance") {
  if (!is.function(fun)) {
    .abort_bad_arg("fun", "must be a function")
  }
  new_distance(
    name = name,
    eval = function(X, x0, ctx = NULL) {
      x_mat <- .validate_matrix_nonempty(X, arg = "X")
      x0_vec <- .validate_target_vector(x0, ncol(x_mat), x0_arg = "x0", x_arg = "X")

      d <- .distance_call_user(fun, x_mat, x0_vec, ctx)
      .distance_validate_output_vec(d, nrow(x_mat))
    }
  )
}

# --- registry -----------------------------------------------------------------

#' Register a distance spec under a name
#'
#' @param name Character scalar registry key.
#' @param spec Either a distance spec or a function `(X, x0, ctx) -> numeric(n)`.
#' @return (Invisibly) the registered spec.
#' @examples
#' nm <- "tmp_l1_demo"
#' register_distance(
#'   nm,
#'   function(X, x0, ctx = NULL) {
#'     rowSums(abs(X - matrix(x0, nrow(X), ncol(X), TRUE)))
#'   }
#' )
#' nm %in% distance_names()
#' @export
register_distance <- function(name, spec) {
  name <- .validate_registry_name(name)
  if (is.function(spec)) spec <- as_distance(spec, name = name)
  if (is.list(spec) && (is.null(spec$name) ||
    (is.character(spec$name) && length(spec$name) == 1L &&
      (is.na(spec$name) || !nzchar(spec$name))))) {
    spec$name <- name
  }
  spec <- .distance_validate_spec(spec)
  assign(name, spec, envir = .dist_registry)
  invisible(spec)
}

#' Retrieve a registered distance spec
#'
#' @param name Character scalar.
#' @return A distance spec.
#' @examples
#' d <- get_distance("l2")
#' is_distance(d)
#' @export
get_distance <- function(name) {
  name <- .validate_registry_name(name)
  if (!exists(name, envir = .dist_registry, inherits = FALSE)) {
    known <- ls(.dist_registry, all.names = FALSE)
    .abort_bad_arg(
      "name",
      "Unknown distance",
      expected = if (length(known) > 0L) paste(known, collapse = ", ") else "<none registered>",
      actual = as.character(name)
    )
  }
  get(name, envir = .dist_registry, inherits = FALSE)
}

#' List registered distances
#' @return Character vector of names.
#' @examples
#' head(distance_names(), 3)
#' @export
distance_names <- function() {
  sort(ls(.dist_registry, all.names = FALSE))
}

# --- evaluators ---------------------------------------------------------------

#' Evaluate a distance spec at a single target
#'
#' @param spec Distance spec (see [is_distance()]).
#' @param X Numeric matrix \(n  x  p\).
#' @param x0 Numeric vector length \(p\).
#' @param ctx Optional context; if `NULL` and `spec$setup` exists, it will be
#'   called as `spec$setup(X)` to create a context.
#' @return Numeric vector length `n`.
#' @examples
#' set.seed(1)
#' x_mat <- matrix(rnorm(12), 4, 3)
#' x0_vec <- c(0, 0, 0)
#' distance_eval(distance_l2(), x_mat, x0_vec)
#' @export
distance_eval <- function(spec, X, x0, ctx = NULL) {
  if (!is_distance(spec)) .abort_bad_arg("spec", "must be a distance")
  spec <- .distance_validate_spec(spec)
  x_mat <- .validate_matrix_nonempty(X, arg = "X")
  x0_vec <- .validate_target_vector(x0, ncol(x_mat), x0_arg = "x0", x_arg = "X", strict = TRUE, finite = TRUE)
  if (is.null(ctx) && is.function(spec$setup)) ctx <- spec$setup(x_mat)

  has_eval <- .distance_has_eval(spec)
  has_batch <- .distance_has_eval_batch(spec)

  # Use eval() only if it clearly looks like (X, x0, ctx) and not a batch shape
  if (has_eval && .distance_eval_is_single(spec$eval)) {
    out <- spec$eval(x_mat, x0_vec, ctx)
    return(.distance_validate_output_vec(out, nrow(x_mat)))
  }

  if (has_batch) {
    d_mat <- spec$eval_batch(x_mat, matrix(x0_vec, nrow = 1L), ctx)
    d_mat <- .distance_validate_output_mat(d_mat, nrow(x_mat), 1L, src = "distance eval_batch")
    return(as.numeric(d_mat[, 1L]))
  }

  if (has_eval) {
    out <- spec$eval(x_mat, x0_vec, ctx)
    return(.distance_validate_output_vec(out, nrow(x_mat)))
  }

  .abort_bad_arg("spec", "No R `eval` or `eval_batch` provided")
}

#' Evaluate a distance spec for many targets (batch)
#'
#' Scalar-only specifications are supported through the shared fallback.
#' The result retains both matrix dimensions for singleton samples/targets.
#'
#' @param spec Distance spec (see [is_distance()]).
#' @param X Numeric matrix \(n  x  p\).
#' @param X0 Numeric matrix \(m  x  p\) (targets as rows).
#' @param ctx Optional context; if `NULL` and `spec$setup` exists, it will be
#'   called as `spec$setup(X)` to create a context.
#' @return Numeric matrix \(n  x  m\).
#' @examples
#' set.seed(1)
#' x_mat <- matrix(rnorm(12), 4, 3)
#' x0_mat <- matrix(rnorm(6), 2, 3)
#' distance_eval_batch(distance_l2(), x_mat, x0_mat)
#' scalar_l1 <- as_distance(function(X, x0) {
#'   rowSums(abs(sweep(X, 2L, x0, "-")))
#' })
#' dim(distance_eval_batch(scalar_l1, x_mat[1L, , drop = FALSE], x0_mat))
#' @export
distance_eval_batch <- function(spec, X, X0, ctx = NULL) {
  if (!is_distance(spec)) .abort_bad_arg("spec", "must be a distance")
  spec <- .distance_validate_spec(spec)
  x_mat <- .validate_matrix_nonempty(X, arg = "X")
  x0_mat <- .distance_validate_target_matrix(X0, x_mat)

  if (is.null(ctx) && is.function(spec$setup)) ctx <- spec$setup(x_mat)

  has_eval <- .distance_has_eval(spec)
  has_batch <- .distance_has_eval_batch(spec)
  n <- nrow(x_mat)
  m <- nrow(x0_mat)

  if (has_batch) {
    out <- spec$eval_batch(x_mat, x0_mat, ctx)
    return(.distance_validate_output_mat(out, n, m, src = "distance eval_batch"))
  }

  if (has_eval) {
    # Only loop if eval looks like (X, x0, ctx); otherwise, we don't know how to call it
    if (!.distance_eval_is_single(spec$eval)) {
      .abort_bad_arg("spec", "No R `eval_batch` or `eval` provided")
    }
    out <- matrix(NA_real_, n, m)
    for (j in seq_len(m)) {
      out[, j] <- .distance_validate_output_vec(spec$eval(x_mat, x0_mat[j, ], ctx), n, src = "distance eval")
    }
    return(out)
  }

  .abort_bad_arg("spec", "No R `eval_batch` or `eval` provided")
}

# --- built-ins registration (call from .onLoad) -------------------------------

#' @noRd
.register_builtin_distances <- function() {
  # Unscaled
  register_distance("l2", distance_l2())
  register_distance("l2_sq", distance_l2(squared = TRUE))
  register_distance("l1", distance_l1())

  # Optionally: scaled variants later (kept commented until needed)
  register_distance("l1_scaled_sd", distance_l1_scaled("sd"))
  register_distance("l2_scaled_sd", distance_l2_scaled("sd"))

  # Mahalanobis (full sample covariance; setup-then-eval)
  register_distance("mahalanobis", distance_mahalanobis())
}
