#' Score specification (S3) for U-stat kernels
#'
#' A score (kernel) is a symmetric function evaluated on s-sized subsets.
#' Provide either `ptr` (compiled fast path), or an R `eval_batch`, or an R `eval`.
#'
#' @param name Character scalar.
#' @param order Integer order s, or function(n, data, ctx) -> s.
#' @param setup Optional function(data, ...) -> ctx object (list or externalptr).
#' @param eval R fallback: function(idxs, data, ctx, rng) -> numeric(1).
#' @param eval_batch R fast path: function(subs, data, ctx, rng) -> numeric(nrow(subs)).
#' @param ptr External pointer to compiled score (see vignette); if present, used by engine.
#' @param randomized Logical: does the kernel draw randomness?
#' @param symmetric Logical: should be TRUE (the engine assumes permutation symmetry).
#' @return An object of class `ustat_score`.
#' @export
new_score <- function(
  name,
  order,
  setup = NULL,
  eval = NULL,
  eval_batch = NULL,
  ptr = NULL,
  randomized = FALSE,
  symmetric = TRUE
) {
  sc <- list(
    name = as.character(name),
    order = order,
    setup = setup,
    eval = eval,
    eval_batch = eval_batch,
    ptr = ptr,
    randomized = isTRUE(randomized),
    symmetric = isTRUE(symmetric)
  )
  class(sc) <- "ustat_score"
  validate_score(sc)
}

#' @export
print.ustat_score <- function(x, ...) {
  cat("<ustat_score:", x$name, ">\n", sep = "")
  cat("  order: ", if (is.function(x$order)) "function(n, data, ctx)" else x$order, "\n", sep = "")
  cat("  impl : ",
    if (!is.null(x$ptr)) {
      "compiled(ptr) "
    } else if (!is.null(x$eval_batch)) {
      "R eval_batch "
    } else if (!is.null(x$eval)) {
      "R eval "
    } else {
      "NONE"
    },
    if (isTRUE(x$randomized)) "(randomized)" else "",
    "\n",
    sep = ""
  )
  invisible(x)
}

# minimal internal helpers (no extra deps)
.is_scalar_int <- function(x) is.numeric(x) && length(x) == 1L && is.finite(x) && as.integer(x) == x
.is_fun <- function(f, n = 0L) is.function(f) && (length(formals(f)) >= n)

#' Validate a score object
#' @noRd
validate_score <- function(sc) {
  stopifnot(inherits(sc, "ustat_score"))
  if (!nzchar(sc$name)) stop("score name must be non-empty")
  if (!(.is_scalar_int(sc$order) || .is_fun(sc$order, 1L))) {
    stop("order must be integer or function(n, data, ctx)")
  }
  if (is.null(sc$ptr) && is.null(sc$eval_batch) && is.null(sc$eval)) {
    stop("Provide at least one of: ptr, eval_batch, eval")
  }
  if (!isTRUE(sc$symmetric)) warning("score marked non-symmetric; engine assumes symmetry")
  sc
}

# ---------- S3 generics (namespaced to avoid clashing with base::order) ----------

#' Set up score context (precomputation)
#' @param score ustat_score
#' @param data data.frame/matrix
#' @param ... passed to score$setup
#' @export
score_setup <- function(score, data, ...) UseMethod("score_setup")
#' @export
score_setup.ustat_score <- function(score, data, ...) {
  if (is.null(score$setup)) {
    return(NULL)
  }
  score$setup(data, ...)
}

#' Resolve kernel order s
#' @export
score_order <- function(score, n, data, ctx = NULL) UseMethod("score_order")
#' @export
score_order.ustat_score <- function(score, n, data, ctx = NULL) {
  if (.is_scalar_int(score$order)) {
    return(as.integer(score$order))
  }
  score$order(n, data, ctx)
}

#' Evaluate kernel on one subset (R fallback)
#' @param idxs integer vector of length s
#' @param rng optional integer or raw key for RNG (deterministic if provided)
#' @export
score_eval <- function(score, idxs, data, ctx = NULL, rng = NULL) UseMethod("score_eval")
#' @export
score_eval.ustat_score <- function(score, idxs, data, ctx = NULL, rng = NULL) {
  if (is.null(score$eval)) {
    stop("No R `eval` provided; use compiled `ptr` or R `eval_batch`.")
  }
  out <- score$eval(idxs, data, ctx, rng)
  if (!is.numeric(out) || length(out) != 1L || !is.finite(out)) {
    stop("score$eval must return a finite scalar numeric.")
  }
  out
}

#' Evaluate kernel on many subsets (R fast-path)
#' @param subs integer matrix with nrow=K, ncol=s (each row is a subset)
#' @export
score_eval_batch <- function(score, subs, data, ctx = NULL, rng = NULL) UseMethod("score_eval_batch")
#' @export
score_eval_batch.ustat_score <- function(score, subs, data, ctx = NULL, rng = NULL) {
  if (is.null(score$eval_batch)) {
    stop("No R `eval_batch` provided; use compiled `ptr` or R `eval`.")
  }
  out <- score$eval_batch(subs, data, ctx, rng)
  if (!is.numeric(out) || length(out) != nrow(subs) || any(!is.finite(out))) {
    stop("score$eval_batch must return numeric(K) with finite entries.")
  }
  out
}

# ---------- Convenience & predicates ----------

#' Coerce a plain function to a score (simple cases)
#' @param f function(idxs, data, ctx, rng) -> scalar
#' @param s integer order
#' @export
as_score <- function(f, s, name = deparse(substitute(f))) {
  stopifnot(.is_fun(f, 2L), .is_scalar_int(s))
  new_score(name = name, order = as.integer(s), eval = f)
}

#' @export
is_score <- function(x) inherits(x, "ustat_score")

# compiled-pointer predicate the engine can use
#' @export
is_compiled <- function(score) is_score(score) && !is.null(score$ptr)
