#' Distance specification (S3) for DNN weights
#'
#' A distance is evaluated from a data matrix `X` (n x p) to a target point `x0` (length p).
#' Provide either `ptr` (compiled fast path), or an R `eval`, or an R `eval_batch`.
#'
#' @param name Character scalar.
#' @param setup Optional function(data, ...) -> ctx object (list or externalptr).
#'   Use this to precompute scaling, standardization, etc., from `data`.
#' @param eval R fallback: function(X, x0, ctx) -> numeric(n).
#' @param eval_batch R fast path for many targets: function(X, X0, ctx) -> matrix(n, m),
#'   where X0 is an m x p matrix of target points (rows = targets).
#' @param ptr External pointer to compiled distance; if present, used by the engine.
#'   The engine will call a compiled wrapper that conforms to the `eval`/`eval_batch` contract.
#' @param randomized Logical: does the distance use randomness (rare; default FALSE).
#' @return An object of class `dnn_distance`.
#' @export
new_distance <- function(
  name,
  setup = NULL,
  eval = NULL,
  eval_batch = NULL,
  ptr = NULL,
  randomized = FALSE
) {
  dist <- list(
    name = as.character(name),
    setup = setup,
    eval = eval,
    eval_batch = eval_batch,
    ptr = ptr,
    randomized = isTRUE(randomized)
  )
  class(dist) <- "dnn_distance"
  validate_distance(dist)
}

# --- minimal internal helper (no extra deps) ---
.is_fun <- function(f, n = 0L) is.function(f) && (length(formals(f)) >= n)

#' Validate a distance object
#' @noRd
validate_distance <- function(dist) {
  stopifnot(inherits(dist, "dnn_distance"))
  if (!nzchar(dist$name)) stop("distance name must be non-empty")
  if (is.null(dist$ptr) && is.null(dist$eval) && is.null(dist$eval_batch)) {
    stop("Provide at least one of: ptr, eval, eval_batch")
  }
  dist
}

#' @export
print.dnn_distance <- function(x, ...) {
  cat("<dnn_distance: ", x$name, ">\n", sep = "")
  cat("  impl : ",
    if (!is.null(x$ptr)) {
      "compiled(ptr)"
    } else if (!is.null(x$eval_batch)) {
      "R eval_batch"
    } else if (!is.null(x$eval)) {
      "R eval"
    } else {
      "NONE"
    },
    if (isTRUE(x$randomized)) " (randomized)" else "",
    "\n",
    sep = ""
  )
  invisible(x)
}


# ---------- S3 generics (parallel to score.*) ----------

#' Set up distance context (precomputation)
#' @param distance dnn_distance
#' @param data data.frame/matrix used to derive any scaling/normalization
#' @param ... passed to distance$setup
#' @export
distance_setup <- function(distance, data, ...) UseMethod("distance_setup")
#' @export
distance_setup.dnn_distance <- function(distance, data, ...) {
  if (is.null(distance$setup)) {
    return(NULL)
  }
  distance$setup(data, ...)
}

#' Evaluate distance from X to x0 (R fallback or compiled)
#' @param distance dnn_distance
#' @param X numeric matrix (n x p)
#' @param x0 numeric length-p vector
#' @param ctx optional context from distance_setup()
#' @export
distance_eval <- function(distance, X, x0, ctx = NULL) UseMethod("distance_eval")
#' @export
distance_eval.dnn_distance <- function(distance, X, x0, ctx = NULL) {
  if (is.null(distance$eval)) {
    stop("No R `eval` provided; use compiled `ptr` or R `eval_batch`.")
  }
  X <- as.matrix(X)
  x0 <- as.numeric(x0)
  if (ncol(X) == 0L) stop("X must have at least one column.")
  if (length(x0) != ncol(X)) stop("x0 length must equal ncol(X).")
  out <- distance$eval(X, x0, ctx)
  if (!is.numeric(out) || length(out) != nrow(X) || any(!is.finite(out))) {
    stop("distance$eval must return numeric(nrow(X)) with finite entries.")
  }
  out
}

#' Evaluate distance from X to many targets X0 (R fast-path or compiled)
#' @param X0 numeric matrix (m x p); each row is a target
#' @return numeric matrix (n x m), column j = distances to X0(j,)
#' @export
distance_eval_batch <- function(distance, X, X0, ctx = NULL) UseMethod("distance_eval_batch")
#' @export
distance_eval_batch.dnn_distance <- function(distance, X, X0, ctx = NULL) {
  if (is.null(distance$eval_batch)) {
    stop("No R `eval_batch` provided; use compiled `ptr` or R `eval`.")
  }
  X <- as.matrix(X)
  X0 <- as.matrix(X0)
  if (ncol(X) == 0L) stop("X must have at least one column.")
  if (ncol(X0) != ncol(X)) stop("ncol(X0) must match ncol(X).")
  out <- distance$eval_batch(X, X0, ctx)
  if (!is.matrix(out) || nrow(out) != nrow(X) || ncol(out) != nrow(X0) || any(!is.finite(out))) {
    stop("distance$eval_batch must return a finite numeric matrix of size nrow(X) x nrow(X0).")
  }
  out
}

# ---------- Convenience & predicates ----------

#' Coerce a plain function to a distance
#' @param f function(X, x0, ctx) -> numeric(n)
#' @export
as_distance <- function(f, name = deparse(substitute(f))) {
  stopifnot(.is_fun(f, 2L))
  new_distance(name = name, eval = f)
}

#' @export
is_distance <- function(x) inherits(x, "dnn_distance")

# compiled-pointer predicate the engine can use
#' @export
is_compiled_distance <- function(distance) is_distance(distance) && !is.null(distance$ptr)
