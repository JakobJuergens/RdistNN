# R/wrank-DNN-wrapper.R

#' Rank-weight: DNN exact (full U, WOR)
#' @param k Integer, 1 <= k <= s
#' @param s Integer, 1 <= s <= n
#' @return A `rank_weight` object
#' @family rank-weight
#' @export
rank_weight_dnn_exact <- function(k, s) {
  k <- as.integer(k); s <- as.integer(s)
  new_rank_weight(
    name = sprintf("dnn_exact[k=%d,s=%d]", k, s),
    eval = function(ranks, ctx=NULL)
      dnn_weights_exact_full_from_ranks(ranks, k, s)$w,
    eval_batch = function(R, ctx=NULL) {
      n <- nrow(R); m <- ncol(R)
      W <- matrix(NA_real_, n, m)
      for (j in seq_len(m)) W[, j] <- dnn_weights_exact_full_from_ranks(R[, j], k, s)$w
      W
    }
  )
}

#' Rank-weight: DNN asymptotic (full U, WOR)
#' @inheritParams rank_weight_dnn_exact
#' @param normalize logical; renormalize to sum 1
#' @return A `rank_weight` object
#' @family rank-weight
#' @export
rank_weight_dnn_asymp <- function(k, s, normalize = TRUE) {
  k <- as.integer(k); s <- as.integer(s)
  new_rank_weight(
    name = sprintf("dnn_asymp[k=%d,s=%d]", k, s),
    eval = function(ranks, ctx=NULL)
      dnn_weights_asymp_full_from_ranks(ranks, k, s, normalize)$w,
    eval_batch = function(R, ctx=NULL) {
      n <- nrow(R); m <- ncol(R)
      W <- matrix(NA_real_, n, m)
      for (j in seq_len(m)) W[, j] <- dnn_weights_asymp_full_from_ranks(R[, j], k, s, normalize)$w
      W
    }
  )
}
