# R/rank.R

#' Rank-weight specification (S3)
#'
#' A rank-weight takes a *permutation of ranks* `ranks` (length n; values 1..n)
#' and returns weights `w` (length n, typically nonnegative and summing to 1).
#' Provide either `ptr` (compiled fast path), or an R `eval`, or an R `eval_batch`.
#'
#' @param name Character scalar.
#' @param setup Optional function(data, ...) -> ctx object (list or externalptr).
#'   Use this for precomputation (e.g., hyperparameters) based on `data`.
#' @param eval R fallback: function(ranks, ctx) -> numeric(n).
#' @param eval_batch R fast path: function(R, ctx) -> matrix(n, m),
#'   where `R` is an n x m integer matrix of ranks (one column per target).
#' @param ptr External pointer to compiled path (optional).
#' @param randomized Logical: does the weighting use randomness?
#' @return An object of class `rank_weight`.
#' @export
new_rank_weight <- function(
    name,
    setup = NULL,
    eval = NULL,
    eval_batch = NULL,
    ptr = NULL,
    randomized = FALSE
) {
  rw <- list(
    name = as.character(name),
    setup = setup,
    eval = eval,
    eval_batch = eval_batch,
    ptr = ptr,
    randomized = isTRUE(randomized)
  )
  class(rw) <- "rank_weight"
  validate_rank_weight(rw)
}

# Internal helpers -------------------------------------------------------------

.is_fun <- function(f, n = 0L) is.function(f) && (length(formals(f)) >= n)

.check_ranks_vec <- function(ranks) {
  r <- as.integer(ranks)
  n <- length(r)
  if (n == 0L) stop("ranks must have positive length.")
  if (anyNA(r)) stop("ranks cannot contain NA.")
  if (!setequal(r, seq_len(n))) stop("ranks must be a permutation of 1..n.")
  r
}

.check_ranks_mat <- function(R) {
  R <- as.matrix(R)
  if (!is.integer(R)) storage.mode(R) <- "integer"
  if (nrow(R) == 0L || ncol(R) == 0L) stop("R must have positive nrow and ncol.")
  for (j in seq_len(ncol(R))) .check_ranks_vec(R[, j])
  R
}

# Validation & printing --------------------------------------------------------

#' @noRd
validate_rank_weight <- function(rw) {
  stopifnot(inherits(rw, "rank_weight"))
  if (!nzchar(rw$name)) stop("rank-weight name must be non-empty")
  if (is.null(rw$ptr) && is.null(rw$eval) && is.null(rw$eval_batch))
    stop("Provide at least one of: ptr, eval, eval_batch")
  rw
}

#' @export
print.rank_weight <- function(x, ...) {
  cat("<rank_weight: ", x$name, ">\n", sep = "")
  cat("  impl : ",
      if (!is.null(x$ptr))           "compiled(ptr)"
      else if (!is.null(x$eval_batch)) "R eval_batch"
      else if (!is.null(x$eval))     "R eval"
      else                           "NONE",
      if (isTRUE(x$randomized)) " (randomized)" else "",
      "\n", sep = "")
  invisible(x)
}

# Generics ---------------------------------------------------------------------

#' @export
rank_weight_setup <- function(rank_weight, data, ...) UseMethod("rank_weight_setup")
#' @export
rank_weight_setup.rank_weight <- function(rank_weight, data, ...) {
  if (is.null(rank_weight$setup)) return(NULL)
  rank_weight$setup(data, ...)
}

#' Evaluate weights from ranks (single target)
#' @param ranks integer permutation 1..n
#' @export
rank_weight_eval <- function(rank_weight, ranks, ctx = NULL) UseMethod("rank_weight_eval")
#' @export
rank_weight_eval.rank_weight <- function(rank_weight, ranks, ctx = NULL) {
  if (is.null(rank_weight$eval))
    stop("No R `eval` provided; use compiled `ptr` or R `eval_batch`.")
  r <- .check_ranks_vec(ranks)
  w <- rank_weight$eval(r, ctx)
  if (!is.numeric(w) || length(w) != length(r) || any(!is.finite(w)))
    stop("rank_weight$eval must return numeric(length(ranks)) with finite entries.")
  w
}

#' Evaluate weights from ranks (many targets)
#' @param R integer matrix (n x m), each column a permutation 1..n
#' @return matrix (n x m)
#' @export
rank_weight_eval_batch <- function(rank_weight, R, ctx = NULL) UseMethod("rank_weight_eval_batch")
#' @export
rank_weight_eval_batch.rank_weight <- function(rank_weight, R, ctx = NULL) {
  if (is.null(rank_weight$eval_batch))
    stop("No R `eval_batch` provided; use compiled `ptr` or R `eval`.")
  RR <- .check_ranks_mat(R)
  W  <- rank_weight$eval_batch(RR, ctx)
  if (!is.matrix(W) || any(!is.finite(W)) ||
      nrow(W) != nrow(RR) || ncol(W) != ncol(RR))
    stop("rank_weight$eval_batch must return a finite numeric matrix of size nrow(R) x ncol(R).")
  W
}

# Convenience ------------------------------------------------------------------

#' Wrap a plain function into a rank_weight
#' @param f function(ranks, ctx) -> numeric(n)
#' @export
as_rank_weight <- function(f, name = deparse(substitute(f))) {
  stopifnot(.is_fun(f, 1L))
  new_rank_weight(name = name, eval = f)
}

#' Is object a rank_weight?
#' @export
is_rank_weight <- function(x) inherits(x, "rank_weight")

# Utilities --------------------------------------------------------------------

#' Convert an index order to ranks (inverse permutation)
#' @param ord integer vector (1..n) giving increasing order
#' @return integer ranks vector r where r[i] = rank of i
#' @export
ord_to_ranks <- function(ord) {
  ord <- as.integer(ord)
  n <- length(ord)
  if (!setequal(ord, seq_len(n))) stop("ord must be a permutation of 1..n.")
  r <- integer(n)
  r[ord] <- seq_len(n)
  r
}
