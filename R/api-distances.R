# -----------------------------------------------------------------------------
# Unified API + registry for distance specifications

# lightweight registry (name -> spec)
.dist_registry <- new.env(parent = emptyenv())

# --- type + adapters ----------------------------------------------------------

#' Is this a distance spec?
#'
#' A distance spec is a list with `eval(X, x0, ctx)` and/or
#' `eval_batch(X, X0, ctx)`. It may optionally define `setup(X)` to precompute
#' a context (e.g., scaling).
#'
#' @param x Any R object.
#' @return Logical.
#' @export
is_distance <- function(x) {
  inherits(x, "dnn_distance") ||
    (is.list(x) && (is.function(x$eval) || is.function(x$eval_batch)))
}

#' Wrap a function as a distance spec
#'
#' @param fun Function `(X, x0, ctx) -> numeric(n)`.
#' @param name Character name for the spec.
#' @return A distance spec with class `dnn_distance`.
#' @export
as_distance <- function(fun, name = "user_distance") {
  stopifnot(is.function(fun))
  spec <- list(
    name = name,
    eval = function(X, x0, ctx = NULL) {
      X <- as.matrix(X)
      x0 <- as.numeric(x0)
      if (ncol(X) == 0L) stop("X must have at least one column.", call. = FALSE)
      if (length(x0) != ncol(X)) stop("x0 length must equal ncol(X).", call. = FALSE)

      d <- fun(X, x0, ctx)
      if (!is.numeric(d) || length(d) != nrow(X) || any(!is.finite(d))) {
        stop("Distance function must return numeric length nrow(X) with finite entries.", call. = FALSE)
      }
      as.numeric(d)
    }
  )
  class(spec) <- union("dnn_distance", class(spec))
  spec
}

# --- registry -----------------------------------------------------------------

#' Register a distance spec under a name
#'
#' @param name Character scalar.
#' @param spec Either a distance spec or a function `(X, x0, ctx) -> numeric(n)`.
#' @return (Invisibly) the registered spec.
#' @export
register_distance <- function(name, spec) {
  stopifnot(is.character(name), length(name) == 1L)
  if (is.function(spec)) spec <- as_distance(spec, name = name)
  if (!is_distance(spec)) stop("`spec` is not a distance.", call. = FALSE)
  if (is.null(spec$name)) spec$name <- name
  class(spec) <- union("dnn_distance", class(spec))
  assign(name, spec, envir = .dist_registry)
  invisible(spec)
}

#' Retrieve a registered distance spec
#'
#' @param name Character scalar.
#' @return A distance spec.
#' @export
get_distance <- function(name) {
  stopifnot(is.character(name), length(name) == 1L)
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
#' @export
distance_eval <- function(spec, X, x0, ctx = NULL) {
  if (!is_distance(spec)) stop("`spec` must be a distance.", call. = FALSE)
  X <- as.matrix(X)
  x0 <- as.numeric(x0)
  if (ncol(X) == 0L) stop("X must have at least one column.", call. = FALSE)
  if (length(x0) != ncol(X)) stop("x0 length must equal ncol(X).", call. = FALSE)
  if (is.null(ctx) && is.function(spec$setup)) ctx <- spec$setup(X)

  has_eval <- !is.null(spec$eval) && is.function(spec$eval)
  has_batch <- !is.null(spec$eval_batch) && is.function(spec$eval_batch)

  # Use eval() only if it clearly looks like (X, x0, ctx) and not a batch shape
  if (has_eval) {
    fmls <- tryCatch(names(formals(spec$eval)), error = function(e) character())
    single_ok <- length(fmls) > 0L && ("x0" %in% fmls) && !("X0" %in% fmls)
    if (single_ok) {
      return(spec$eval(X, x0, ctx))
    }
    # else: fall through and try batch
  }

  if (has_batch) {
    D <- spec$eval_batch(X, matrix(x0, nrow = 1L), ctx)
    return(as.numeric(D[, 1L]))
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
#' @export
distance_eval_batch <- function(spec, X, X0, ctx = NULL) {
  if (!is_distance(spec)) stop("`spec` must be a distance.", call. = FALSE)
  X <- as.matrix(X)
  X0 <- as.matrix(X0)
  if (ncol(X) == 0L) stop("X must have at least one column.", call. = FALSE)
  if (ncol(X) != ncol(X0)) stop("ncol(X0) must match ncol(X).", call. = FALSE)

  if (is.null(ctx) && is.function(spec$setup)) ctx <- spec$setup(X)

  has_eval <- !is.null(spec$eval) && is.function(spec$eval)
  has_batch <- !is.null(spec$eval_batch) && is.function(spec$eval_batch)

  if (has_batch) {
    return(spec$eval_batch(X, X0, ctx))
  }

  if (has_eval) {
    # Only loop if eval looks like (X, x0, ctx); otherwise, we don't know how to call it
    fmls <- tryCatch(names(formals(spec$eval)), error = function(e) character())
    single_ok <- length(fmls) > 0L && ("x0" %in% fmls) && !("X0" %in% fmls)
    if (!single_ok) {
      stop("No R `eval_batch` or `eval` provided.", call. = FALSE)
    }
    m <- nrow(X0)
    out <- matrix(NA_real_, nrow(X), m)
    for (j in seq_len(m)) out[, j] <- spec$eval(X, X0[j, ], ctx)
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
}
