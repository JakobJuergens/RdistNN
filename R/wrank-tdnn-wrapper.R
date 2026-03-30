# Rank-weight spec constructors for TDNN (Demirkaya two-scale DNN, k = 1).
#
# These return rank_weight specs compatible with the rank API:
#   rank_weight_eval(), rank_weight_eval_batch(), weights_matrix(), etc.
#
# They use a table-based evaluation strategy: for a given n, the TDNN
# rank-order weight table is computed once and cached in the spec closure,
# then W[, t] = w_rank[R[, t]] for each target column t.  This avoids
# recomputing the table for every column of a rank matrix.
#
# Both exact and asymptotic variants share the same interface and the same
# caching/dispatch pattern; they differ only in the underlying rank table.
#
# They rely on:
#   tdnn_rank_table_exact(n, s1, s2, dim)    -> w_rank[1:n]
#   tdnn_rank_table_asymp(n, s1, s2, dim, normalize) -> w_rank[1:n]

#' Rank-weight spec: TDNN exact (Demirkaya two-scale DNN, k = 1)
#'
#' Returns a `rank_weight` spec whose `eval_batch` uses the exact TDNN rank
#' table from [tdnn_rank_table_exact()].  The rank table is computed once per
#' sample size `n` and cached inside the spec closure, so calling
#' [rank_weight_eval_batch()] repeatedly on rank matrices of the same size
#' does not recompute the table.
#'
#' **Signed weights.** TDNN weights are a signed linear combination of two
#' DNN tables and can be negative.  They sum to 1 across all `n`
#' observations but the negative entries are part of the bias-correction
#' mechanism and should not be interpreted as a standard probability mass.
#'
#' @param s1 Positive integer; smaller U-statistic order (\code{s1 < s2}).
#' @param s2 Positive integer; larger U-statistic order.
#' @param dim Positive integer; localization dimension.  For standard Euclidean
#'   ranking this is \code{ncol(X)}; for custom distances it should be set
#'   explicitly.
#'
#' @return A \code{rank_weight} spec (list with class \code{"rank_weight"}).
#'
#' @seealso [rank_weight_tdnn_asymp()], [tdnn_rank_table_exact()],
#'   [weights_tdnn_exact()]
#'
#' @examples
#' rw <- rank_weight_tdnn_exact(s1 = 3L, s2 = 6L, dim = 2L)
#' R <- matrix(c(2L, 1L, 3L, 4L, 5L,
#'               5L, 2L, 1L, 3L, 4L), nrow = 5)
#' W <- rank_weight_eval_batch(rw, R)
#' colSums(W)  # c(1, 1)
#' @export
rank_weight_tdnn_exact <- function(s1, s2, dim) {
  s1  <- as.integer(s1)
  s2  <- as.integer(s2)
  dim <- as.integer(dim)

  # Per-n table cache (keyed by character n; lives inside the spec closure)
  .cache <- new.env(parent = emptyenv())

  .get_table <- function(n) {
    key <- as.character(n)
    if (!exists(key, envir = .cache, inherits = FALSE)) {
      assign(key, tdnn_rank_table_exact(n, s1, s2, dim), envir = .cache)
    }
    get(key, envir = .cache, inherits = FALSE)
  }

  eval_batch <- function(R, ctx = NULL) {
    n <- nrow(R)
    m <- ncol(R)
    w_rank <- .get_table(n)
    # Table-based: index into w_rank once for the entire rank matrix
    W <- w_rank[R]
    dim(W) <- c(n, m)
    W
  }

  spec <- list(
    name       = sprintf("tdnn_exact[s1=%d,s2=%d,dim=%d]", s1, s2, dim),
    eval_batch = eval_batch
  )
  spec$.__kind <- "rank_weight"
  class(spec)  <- c("rank_weight", class(spec))
  spec
}

#' Rank-weight spec: TDNN asymptotic (Demirkaya two-scale DNN, k = 1)
#'
#' Returns a `rank_weight` spec whose `eval_batch` uses the asymptotic TDNN
#' rank table from [tdnn_rank_table_asymp()].  Structurally identical to
#' [rank_weight_tdnn_exact()] but uses the Beta order-statistic approximation,
#' which gives full rank support (all n entries are non-zero).
#'
#' @param s1 Positive integer; smaller U-statistic order (\code{s1 < s2}).
#' @param s2 Positive integer; larger U-statistic order.
#' @param dim Positive integer; localization dimension.
#' @param normalize Logical; renormalize the rank table to sum 1
#'   (default \code{TRUE}; guards floating-point drift).
#'
#' @return A \code{rank_weight} spec (list with class \code{"rank_weight"}).
#'
#' @seealso [rank_weight_tdnn_exact()], [tdnn_rank_table_asymp()],
#'   [weights_tdnn_asymp()]
#'
#' @examples
#' rw <- rank_weight_tdnn_asymp(s1 = 3L, s2 = 6L, dim = 2L)
#' R <- matrix(c(2L, 1L, 3L, 4L, 5L,
#'               5L, 2L, 1L, 3L, 4L), nrow = 5)
#' W <- rank_weight_eval_batch(rw, R)
#' colSums(W)  # c(1, 1)
#' @export
rank_weight_tdnn_asymp <- function(s1, s2, dim, normalize = TRUE) {
  s1        <- as.integer(s1)
  s2        <- as.integer(s2)
  dim       <- as.integer(dim)
  normalize <- isTRUE(normalize)

  .cache <- new.env(parent = emptyenv())

  .get_table <- function(n) {
    key <- as.character(n)
    if (!exists(key, envir = .cache, inherits = FALSE)) {
      assign(key, tdnn_rank_table_asymp(n, s1, s2, dim, normalize = normalize),
             envir = .cache)
    }
    get(key, envir = .cache, inherits = FALSE)
  }

  eval_batch <- function(R, ctx = NULL) {
    n <- nrow(R)
    m <- ncol(R)
    w_rank <- .get_table(n)
    W <- w_rank[R]
    dim(W) <- c(n, m)
    W
  }

  spec <- list(
    name = sprintf("tdnn_asymp[s1=%d,s2=%d,dim=%d,norm=%s]",
                   s1, s2, dim, if (normalize) "T" else "F"),
    eval_batch = eval_batch
  )
  spec$.__kind <- "rank_weight"
  class(spec)  <- c("rank_weight", class(spec))
  spec
}
