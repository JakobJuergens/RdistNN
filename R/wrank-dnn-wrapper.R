# Rank-weight wrappers for DNN weights compatible with the new rank API.
# These return a rdistnn_rank_weight spec (list) with fields:
#   name, table(n), eval(ranks, ctx), eval_batch(R, ctx),
#   and the marker: spec$.__kind <- "rank_weight"; class(spec) <- "rdistnn_rank_weight".
#
# They rely on:
#   dnn_weights_exact_full_from_ranks(ranks, k, s)$w
#   dnn_weights_asymp_full_from_ranks(ranks, k, s, normalize)$w

#' Rank-weight: DNN exact (full U, WOR)
#' @param k Integer; 1 <= k <= s
#' @param s Integer; U-statistic order
#' @return rdistnn_rank_weight spec (list)
#' @examples
#' rw <- rank_weight_dnn_exact(k = 1, s = 3)
#' rank_weight_eval(rw, c(2L, 1L, 3L, 4L))
#' @export
rank_weight_dnn_exact <- function(k, s) {
  k <- .validate_scalar_int(k, "k")
  s <- .validate_scalar_int(s, "s")
  if (k > s) .abort_bad_arg("k", sprintf("must satisfy k <= s (got k=%d, s=%d)", k, s))

  table_one <- function(n) {
    n <- .validate_scalar_int(n, "n")
    dnn_weights_exact_full_from_ranks(seq_len(n), k, s)$w_rank
  }

  eval_one <- function(ranks, ctx = NULL) {
    # ranks is already validated upstream by rank_weight_eval(_batch)
    table_one(length(ranks))[ranks]
  }

  eval_many <- function(R, ctx = NULL) {
    n <- nrow(R)
    m <- ncol(R)
    matrix(table_one(n)[R], nrow = n, ncol = m)
  }

  spec <- list(
    name = sprintf("dnn_exact[k=%d,s=%d]", k, s),
    table = table_one,
    eval = eval_one,
    eval_batch = eval_many
  )
  spec$.__kind <- "rank_weight"
  class(spec) <- c("rdistnn_rank_weight", class(spec))
  spec
}

#' Rank-weight: DNN asymptotic (full U, WOR)
#' @inheritSection weights_dnn_asymp Numerical limits
#' @inheritParams rank_weight_dnn_exact
#' @param normalize logical; renormalize to sum 1 (default TRUE)
#' @return rdistnn_rank_weight spec (list)
#' @examples
#' rw <- rank_weight_dnn_asymp(k = 1, s = 3, normalize = TRUE)
#' rank_weight_eval(rw, c(2L, 1L, 3L, 4L))
#' @export
rank_weight_dnn_asymp <- function(k, s, normalize = TRUE) {
  k <- .validate_scalar_int(k, "k")
  s <- .validate_scalar_int(s, "s")
  if (k > s) .abort_bad_arg("k", sprintf("must satisfy k <= s (got k=%d, s=%d)", k, s))
  normalize <- .validate_scalar_logical(normalize, "normalize")

  table_one <- function(n) {
    n <- .validate_scalar_int(n, "n")
    w <- dnn_weights_asymp_full_from_ranks(seq_len(n), k, s, normalize)$w_rank
    w <- as.numeric(w)
    sw <- .compensated_sum(w)
    if (!is.finite(sw) || sw <= 0) .dnn_abort("numerical_range",
      "Asymptotic DNN weights have a nonpositive or nonfinite computed total.",
      hint = "Use less extreme feasible orders, or report a reproducible example if this persists."
    )
    if (normalize) w / sw else w
  }

  eval_one <- function(ranks, ctx = NULL) {
    table_one(length(ranks))[ranks]
  }

  eval_many <- function(R, ctx = NULL) {
    n <- nrow(R)
    m <- ncol(R)
    matrix(table_one(n)[R], nrow = n, ncol = m)
  }

  spec <- list(
    name = sprintf("dnn_asymp[k=%d,s=%d,norm=%s]", k, s, if (normalize) "T" else "F"),
    table = table_one,
    eval = eval_one,
    eval_batch = eval_many
  )
  spec$.__kind <- "rank_weight"
  class(spec) <- c("rdistnn_rank_weight", class(spec))
  spec
}
