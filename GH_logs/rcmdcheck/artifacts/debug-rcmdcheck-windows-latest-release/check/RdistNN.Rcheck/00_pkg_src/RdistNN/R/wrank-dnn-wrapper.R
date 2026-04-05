# Rank-weight wrappers for DNN weights compatible with the new rank API.
# These return a rank_weight spec (list) with fields:
#   name, eval(ranks, ctx), eval_batch(R, ctx),
#   and the marker: spec$.__kind <- "rank_weight"; class(spec) <- "rank_weight".
#
# They rely on:
#   dnn_weights_exact_full_from_ranks(ranks, k, s)$w
#   dnn_weights_asymp_full_from_ranks(ranks, k, s, normalize)$w

#' Rank-weight: DNN exact (full U, WOR)
#' @param k Integer; 1 <= k <= s
#' @param s Integer; U-statistic order
#' @return rank_weight spec (list)
#' @examples
#' rw <- rank_weight_dnn_exact(k = 2, s = 3)
#' rank_weight_eval(rw, c(2L, 1L, 3L, 4L))
#' @export
rank_weight_dnn_exact <- function(k, s) {
  k <- as.integer(k)
  s <- as.integer(s)

  eval_one <- function(ranks, ctx = NULL) {
    # ranks is already validated upstream by rank_weight_eval(_batch)
    w <- dnn_weights_exact_full_from_ranks(ranks, k, s)$w
    # be defensive about normalization
    w <- as.numeric(w)
    sw <- sum(w)
    if (!is.finite(sw) || sw <= 0) stop("DNN exact produced nonpositive/invalid sum.", call. = FALSE)
    w / sw
  }

  eval_many <- function(R, ctx = NULL) {
    n <- nrow(R)
    m <- ncol(R)
    W <- matrix(NA_real_, n, m)
    for (j in seq_len(m)) W[, j] <- eval_one(R[, j], ctx)
    W
  }

  spec <- list(
    name = sprintf("dnn_exact[k=%d,s=%d]", k, s),
    eval = eval_one,
    eval_batch = eval_many
  )
  spec$.__kind <- "rank_weight"
  class(spec) <- c("rank_weight", class(spec))
  spec
}

#' Rank-weight: DNN asymptotic (full U, WOR)
#' @inheritParams rank_weight_dnn_exact
#' @param normalize logical; renormalize to sum 1 (default TRUE)
#' @return rank_weight spec (list)
#' @examples
#' rw <- rank_weight_dnn_asymp(k = 2, s = 3, normalize = TRUE)
#' rank_weight_eval(rw, c(2L, 1L, 3L, 4L))
#' @export
rank_weight_dnn_asymp <- function(k, s, normalize = TRUE) {
  k <- as.integer(k)
  s <- as.integer(s)
  normalize <- isTRUE(normalize)

  eval_one <- function(ranks, ctx = NULL) {
    w <- dnn_weights_asymp_full_from_ranks(ranks, k, s, normalize)$w
    w <- as.numeric(w)
    sw <- sum(w)
    if (!is.finite(sw) || sw <= 0) stop("DNN asymp produced nonpositive/invalid sum.", call. = FALSE)
    w / sw
  }

  eval_many <- function(R, ctx = NULL) {
    n <- nrow(R)
    m <- ncol(R)
    W <- matrix(NA_real_, n, m)
    for (j in seq_len(m)) W[, j] <- eval_one(R[, j], ctx)
    W
  }

  spec <- list(
    name = sprintf("dnn_asymp[k=%d,s=%d,norm=%s]", k, s, if (normalize) "T" else "F"),
    eval = eval_one,
    eval_batch = eval_many
  )
  spec$.__kind <- "rank_weight"
  class(spec) <- c("rank_weight", class(spec))
  spec
}
