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
    stop(sprintf("%s must return numeric length nrow(X) with finite entries.", src), call. = FALSE)
  }
  as.numeric(d)
}
.distance_validate_output_mat <- function(d_mat, n, m, src = "Distance function") {
  if (!is.matrix(d_mat) || !is.numeric(d_mat) || nrow(d_mat) != n || ncol(d_mat) != m || any(!is.finite(d_mat))) {
    stop(sprintf("%s must return a finite numeric matrix with shape nrow(X) x nrow(X0).", src), call. = FALSE)
  }
  d_mat
}
.distance_validate_spec <- function(spec) {
  if (!is.list(spec)) stop("`spec` must be a list.", call. = FALSE)
  if (!(.distance_has_eval(spec) || .distance_has_eval_batch(spec))) {
    stop("Distance spec must provide at least one of `eval` or `eval_batch`.", call. = FALSE)
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
    name = as.character(name),
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
#' Accepts either `(name, spec)` or legacy single-argument form `(spec)`,
#' where `spec$name` is used as the key.
#'
#' @param name Character scalar, or a spec in legacy single-argument form.
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
register_distance <- function(name, spec = NULL) {
  if (missing(spec)) {
    spec <- name
    if (is.function(spec)) {
      name <- "user_distance"
      spec <- as_distance(spec, name = name)
    } else if (is.list(spec) && is.character(spec$name) && length(spec$name) == 1L) {
      name <- spec$name
    } else {
      stop("Legacy register_distance(spec): spec must have a non-empty $name.", call. = FALSE)
    }
  } else {
    if (!is.character(name) || length(name) != 1L) {
      .abort_bad_arg("name", "must be a character scalar")
    }
    if (is.function(spec)) spec <- as_distance(spec, name = name)
    if (is.list(spec) && (is.null(spec$name) || !nzchar(as.character(spec$name)))) {
      spec$name <- name
    }
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
  if (!is.character(name) || length(name) != 1L) {
    .abort_bad_arg("name", "must be a character scalar")
  }
  if (!exists(name, envir = .dist_registry, inherits = FALSE)) {
    known <- ls(.dist_registry, all.names = FALSE)
    stop(sprintf(
      "Unknown distance '%s'. Known: %s",
      name, paste(known, collapse = ", ")
    ), call. = FALSE)
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

#' Back-compat alias for [distance_names()]
#' @return Character vector of names.
#' @examples
#' identical(distance_names(), list_distances())
#' @export
list_distances <- function() distance_names()

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
#' x_mat <- matrix(rnorm(12), 4, 3)
#' x0_vec <- c(0, 0, 0)
#' distance_eval(distance_l2(), x_mat, x0_vec)
#' @export
distance_eval <- function(spec, X, x0, ctx = NULL) {
  if (!is_distance(spec)) stop("`spec` must be a distance.", call. = FALSE)
  x_mat <- .validate_matrix_nonempty(X, arg = "X")
  x0_vec <- .validate_target_vector(x0, ncol(x_mat), x0_arg = "x0", x_arg = "X")
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

  stop("No R `eval` or `eval_batch` provided.", call. = FALSE)
}

#' Evaluate a distance spec for many targets (batch)
#'
#' @param spec Distance spec (see [is_distance()]).
#' @param X Numeric matrix \(n  x  p\).
#' @param X0 Numeric matrix \(m  x  p\) (targets as rows).
#' @param ctx Optional context; if `NULL` and `spec$setup` exists, it will be
#'   called as `spec$setup(X)` to create a context.
#' @return Numeric matrix \(n  x  m\).
#' @examples
#' x_mat <- matrix(rnorm(12), 4, 3)
#' x0_mat <- matrix(rnorm(6), 2, 3)
#' distance_eval_batch(distance_l2(), x_mat, x0_mat)
#' @export
distance_eval_batch <- function(spec, X, X0, ctx = NULL) {
  if (!is_distance(spec)) stop("`spec` must be a distance.", call. = FALSE)
  x_mat <- .validate_matrix_nonempty(X, arg = "X")
  x0_mat <- .validate_ncol_relation(x_mat, X0, y_arg = "X0", x_arg = "X", relation = "match")

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
      stop("No R `eval_batch` or `eval` provided.", call. = FALSE)
    }
    out <- matrix(NA_real_, n, m)
    for (j in seq_len(m)) {
      out[, j] <- .distance_validate_output_vec(spec$eval(x_mat, x0_mat[j, ], ctx), n, src = "distance eval")
    }
    return(out)
  }

  stop("No R `eval_batch` or `eval` provided.", call. = FALSE)
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
