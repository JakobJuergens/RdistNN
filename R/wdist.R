# R/wdist.R
# -----------------------------------------------------------------------------
# Distance-weight specification (S3)
# A distance-weight consumes distances dx (length n) or a distance matrix DX
# (n x m; n = reference rows, m = targets) and returns weights that typically
# sum to 1 per target.
# -----------------------------------------------------------------------------

#' Distance-weight specification (S3)
#'
#' A distance-weight takes a numeric distance vector \code{dx} (length \code{n})
#' or a distance matrix \code{DX} (\code{n x m}) and returns weights \code{w}
#' (length \code{n}) or a weight matrix (\code{n x m}). Provide either a compiled
#' \code{ptr} (fast path), an R \code{eval_dx} for single-target evaluation, or an
#' R \code{eval_dx_batch} for many targets.
#'
#' @param name Character scalar.
#' @param setup Optional function(data, ...) -> ctx object (list or externalptr),
#'   used to precompute bandwidths, scaling, etc., from training data.
#' @param eval_dx R fallback: function(dx, ctx) -> numeric(n).
#' @param eval_dx_batch R fast path: function(DX, ctx) -> matrix(n, m).
#' @param ptr External pointer to compiled implementation (optional).
#' @param randomized Logical: does the weighting rule use randomness?
#' @return An object of class \code{dist_weight}.
#' @export
new_dist_weight <- function(
  name,
  setup = NULL,
  eval_dx = NULL,
  eval_dx_batch = NULL,
  ptr = NULL,
  randomized = FALSE
) {
  dw <- list(
    name = as.character(name),
    setup = setup,
    eval_dx = eval_dx,
    eval_dx_batch = eval_dx_batch,
    ptr = ptr,
    randomized = isTRUE(randomized)
  )
  class(dw) <- "dist_weight"
  validate_dist_weight(dw)
}

# --- internal helpers (no extra deps) -----------------------------------------

.is_fun <- function(f, n = 0L) is.function(f) && (length(formals(f)) >= n)

.check_dx_vec <- function(dx) {
  # NEW: forbid matrices/arrays on the single-target path
  if (!is.null(dim(dx))) {
    stop("dx must be a numeric vector; use dist_weight_eval_batch() for matrices.")
  }
  x <- as.numeric(dx)
  n <- length(x)
  if (n == 0L) stop("dx must have positive length.")
  if (anyNA(x)) stop("dx cannot contain NA.")
  if (any(!is.finite(x))) stop("dx must be finite.")
  x
}

.check_dx_mat <- function(DX) {
  M <- as.matrix(DX)
  if (!is.numeric(M)) storage.mode(M) <- "double"
  if (nrow(M) == 0L || ncol(M) == 0L) stop("DX must have positive nrow and ncol.")
  if (any(!is.finite(M))) stop("DX must be finite.")
  M
}

# --- validation & printing ----------------------------------------------------

#' @noRd
validate_dist_weight <- function(dw) {
  stopifnot(inherits(dw, "dist_weight"))
  if (!nzchar(dw$name)) stop("distance-weight name must be non-empty")
  if (is.null(dw$ptr) && is.null(dw$eval_dx) && is.null(dw$eval_dx_batch)) {
    stop("Provide at least one of: ptr, eval_dx, eval_dx_batch")
  }
  dw
}

#' @export
print.dist_weight <- function(x, ...) {
  cat("<dist_weight: ", x$name, ">\n", sep = "")
  cat("  impl : ",
    if (!is.null(x$ptr)) {
      "compiled(ptr)"
    } else if (!is.null(x$eval_dx_batch)) {
      "R eval_dx_batch"
    } else if (!is.null(x$eval_dx)) {
      "R eval_dx"
    } else {
      "NONE"
    },
    if (isTRUE(x$randomized)) " (randomized)" else "",
    "\n",
    sep = ""
  )
  invisible(x)
}

# --- S3 generics --------------------------------------------------------------

#' Set up distance-weight context (precomputation)
#' @param dist_weight dist_weight
#' @param data data.frame/matrix used to derive bandwidths/scaling
#' @param ... passed to dist_weight$setup
#' @export
dist_weight_setup <- function(dist_weight, data, ...) UseMethod("dist_weight_setup")
#' @export
dist_weight_setup.dist_weight <- function(dist_weight, data, ...) {
  if (is.null(dist_weight$setup)) {
    return(NULL)
  }
  dist_weight$setup(data, ...)
}

#' Evaluate distance-weights for one target
#' @param dist_weight dist_weight
#' @param dx numeric vector of distances (length n)
#' @param ctx optional context
#' @return numeric(n)
#' @export
dist_weight_eval <- function(dist_weight, dx, ctx = NULL) UseMethod("dist_weight_eval")
#' @export
dist_weight_eval.dist_weight <- function(dist_weight, dx, ctx = NULL) {
  if (is.null(dist_weight$eval_dx)) {
    stop("No R `eval_dx` provided; use compiled `ptr` or R `eval_dx_batch`.")
  }
  x <- .check_dx_vec(dx)
  w <- dist_weight$eval_dx(x, ctx)
  if (!is.numeric(w) || length(w) != length(x) || any(!is.finite(w))) {
    stop("dist_weight$eval_dx must return numeric(length(dx)) with finite entries.")
  }
  w
}

#' Evaluate distance-weights for many targets
#' @param dist_weight dist_weight
#' @param DX numeric matrix (n x m) of distances
#' @param ctx optional context
#' @return numeric matrix (n x m)
#' @export
dist_weight_eval_batch <- function(dist_weight, DX, ctx = NULL) UseMethod("dist_weight_eval_batch")
#' @export
dist_weight_eval_batch.dist_weight <- function(dist_weight, DX, ctx = NULL) {
  if (is.null(dist_weight$eval_dx_batch)) {
    stop("No R `eval_dx_batch` provided; use compiled `ptr` or R `eval_dx`.")
  }
  M <- .check_dx_mat(DX)
  W <- dist_weight$eval_dx_batch(M, ctx)
  if (!is.matrix(W) || nrow(W) != nrow(M) || ncol(W) != ncol(M) || any(!is.finite(W))) {
    stop("dist_weight$eval_dx_batch must return a finite numeric matrix of size nrow(DX) x ncol(DX).")
  }
  W
}

# --- convenience & predicates -------------------------------------------------

#' Wrap a plain function into a dist_weight (simple cases)
#' @param f function(dx, ctx) -> numeric(n)
#' @param name character scalar used as the dist-weight name (optional)
#' @return A \code{dist_weight}
#' @export
as_dist_weight <- function(f, name = deparse(substitute(f))) {
  stopifnot(.is_fun(f, 1L))
  new_dist_weight(name = name, eval_dx = f)
}

#' Is object a dist_weight?
#' Is object a \code{dist_weight}?
#' @param x any object
#' @return logical
#' @export
is_dist_weight <- function(x) inherits(x, "dist_weight")

#' Does a dist_weight have a compiled implementation?
#' @param dw a \code{dist_weight}
#' @return logical
#' @export
is_compiled_dist_weight <- function(dw) is_dist_weight(dw) && !is.null(dw$ptr)
