# Asymptotic TDNN rank-weight table (Demirkaya two-scale DNN, k = 1)
#
# Identical TDNN construction as the exact variant but using the Beta
# order-statistic approximation for the underlying single-scale tables:
#
#   w_rank^DNN_asymp[j] = F_{Beta(1, s)}(j/n) - F_{Beta(1, s)}((j-1)/n)
#                       = pbeta(j/n, 1, s) - pbeta((j-1)/n, 1, s)
#
# For k = 1 the inner sum over i collapses to a single Beta(1, s) term.
# The TDNN asymptotic table is then:
#
#   w_rank^TDNN[j] = w1* * w_rank^DNN_asymp(s1)[j]
#                  + w2* * w_rank^DNN_asymp(s2)[j]
#
# Like the exact variant these weights are signed and sum to 1.
# Unlike the exact variant the asymptotic table has full support (all
# n ranks receive non-zero weight, though many will be very small).
# --------------------------------------------------------------------------

# ---- Internal: single-scale asymptotic rank table for k = 1 -------------

.dnn_rank_table_asymp_k1 <- function(n, s, normalize = TRUE) {
  # Compute the asymptotic DNN rank-order weight table for k = 1.
  # For k=1 the Beta-increment sum reduces to a single Beta(1, s) term:
  #   w_rank[j] = pbeta(j/n, 1, s) - pbeta((j-1)/n, 1, s)
  n <- as.integer(n)
  s <- as.integer(s)
  j <- seq_len(n)
  w_rank <- .beta_interval_mass((j - 1L) / n, j / n, 1L, s)
  w_rank <- as.numeric(w_rank)
  if (isTRUE(normalize)) {
    sw <- .compensated_sum(w_rank)
    if (!is.finite(sw) || sw <= 0) {
      stop("Internal error: non-positive total weight in .dnn_rank_table_asymp_k1.", call. = FALSE)
    }
    if (abs(sw - 1) > 1e-12) w_rank <- w_rank / sw
  }
  w_rank
}

# ---- Public: TDNN asymptotic rank table ----------------------------------

#' Asymptotic TDNN rank-weight table (Demirkaya two-scale DNN, k = 1)
#'
#' Computes the rank-order weight sequence for the two-scale DNN (TDNN)
#' estimator using the Beta order-statistic approximation for the underlying
#' single-scale tables.  Structurally identical to `tdnn_rank_table_exact()`
#' but evaluates the asymptotic approximation across all `n` ranks.
#' Combined entries can be zero through signed cancellation or underflow.
#'
#' \deqn{w^{\text{TDNN,asymp}}_j = w_1^* \, \Delta_j\text{Beta}(1, s_1)
#'   + w_2^* \, \Delta_j\text{Beta}(1, s_2)}
#'
#' where \eqn{\Delta_j\text{Beta}(1, s) = F_{\text{Beta}(1,s)}(j/n) -
#' F_{\text{Beta}(1,s)}((j-1)/n)} and the bias-correction scalars
#' \eqn{w_1^*, w_2^*} are the same as in `tdnn_rank_table_exact()`.
#'
#' @param n Positive integer; sample size.
#' @param s1 Positive integer; smaller U-statistic order (\code{s1 < s2 <= n}).
#' @param s2 Positive integer; larger U-statistic order (\code{s1 < s2 <= n}).
#' @param dim Positive integer; localization dimension (typically \code{ncol(X)}).
#' @param normalize Logical; renormalize so entries sum exactly to 1
#'   (default \code{TRUE}; guards floating-point drift).
#'
#' @return A numeric vector of length \code{n} giving the rank-order TDNN
#'   asymptotic weights \code{w_rank[1:n]}.  Entries sum to 1 but some can
#'   be negative.
#'
#' @seealso `tdnn_rank_table_exact()`, `dnn_weights_tdnn_asymp_full_from_ranks()`,
#'   [rank_weight_tdnn_asymp()]
#'
#' @examples
#' n <- 20; s1 <- 3L; s2 <- 6L; d <- 2L
#' w <- tdnn_rank_table_asymp(n, s1, s2, d)
#' sum(w)      # 1
#' any(w < 0)  # TRUE
#' @keywords internal
#' @noRd
tdnn_rank_table_asymp <- function(n, s1, s2, dim, normalize = TRUE) {
  n   <- as.integer(n)
  s1  <- as.integer(s1)
  s2  <- as.integer(s2)
  dim <- as.integer(dim)

  if (n < 1L) stop("`n` must be a positive integer.", call. = FALSE)
  if (s1 < 1L || s1 >= s2) stop("`s1` must satisfy 1 <= s1 < s2.", call. = FALSE)
  if (s2 > n)  stop("`s2` must satisfy s2 <= n.", call. = FALSE)
  if (dim < 1L) stop("`dim` must be a positive integer.", call. = FALSE)
  normalize <- .validate_scalar_logical(normalize, "normalize")

  sc <- .tdnn_scalars(s1, s2, dim)
  t1 <- .dnn_rank_table_asymp_k1(n, s1, normalize = FALSE)
  t2 <- .dnn_rank_table_asymp_k1(n, s2, normalize = FALSE)
  w_rank <- sc["w1"] * t1 + sc["w2"] * t2

  if (normalize) {
    sw <- .compensated_sum(w_rank)
    if (!is.finite(sw) || abs(sw) < .Machine$double.eps * 100) {
      stop("Internal error: TDNN asymptotic rank table has zero/non-finite sum.", call. = FALSE)
    }
    if (abs(sw - 1) > 1e-10) w_rank <- w_rank / sw
  }
  w_rank
}

# ---- Public: computation function (mirrors dnn_weights_asymp_full_from_ranks) ----

#' Asymptotic TDNN weights from ranks (Demirkaya two-scale DNN, k = 1)
#'
#' Given a rank permutation `ranks`, computes the TDNN weight for each
#' observation using the asymptotic rank-weight table from
#' `tdnn_rank_table_asymp()`.
#'
#' @param ranks Integer permutation of \code{1:n}; \code{ranks[i]} is the
#'   rank (1 = nearest) of observation \code{i} relative to the target.
#' @param s1 Positive integer; smaller U-statistic order.
#' @param s2 Positive integer; larger U-statistic order.
#' @param dim Positive integer; localization dimension.
#' @param normalize Logical; renormalize to sum 1 (default \code{TRUE}).
#'
#' @return A list with
#' \itemize{
#'   \item \code{w}: numeric length-\code{n} weights in data order (sum to 1,
#'     some entries may be negative).
#'   \item \code{w_rank}: numeric length-\code{n} weights in rank order (1..n).
#'   \item \code{meta}: list with \code{n}, \code{s1}, \code{s2}, \code{dim},
#'     \code{normalize}, and \code{scalars}.
#' }
#'
#' @seealso `tdnn_rank_table_asymp()`, [rank_weight_tdnn_asymp()],
#'   `dnn_weights_tdnn_exact_full_from_ranks()`
#'
#' @examples
#' set.seed(1)
#' n <- 15L; ranks <- sample.int(n)
#' out <- dnn_weights_tdnn_asymp_full_from_ranks(ranks, s1 = 3L, s2 = 6L, dim = 2L)
#' sum(out$w)        # 1
#' any(out$w < 0)    # TRUE
#' @keywords internal
#' @noRd
dnn_weights_tdnn_asymp_full_from_ranks <- function(ranks, s1, s2, dim,
                                                    normalize = TRUE) {
  if (is.null(ranks)) stop("`ranks` must be provided.", call. = FALSE)
  if (!is.integer(ranks)) ranks <- as.integer(ranks)
  n <- length(ranks)
  if (n < 1L) stop("`ranks` must have positive length.", call. = FALSE)
  if (anyNA(ranks)) stop("`ranks` cannot contain NA.", call. = FALSE)
  if (!setequal(ranks, seq_len(n))) stop("`ranks` must be a permutation of 1..n.", call. = FALSE)

  s1  <- as.integer(s1)
  s2  <- as.integer(s2)
  dim <- as.integer(dim)
  if (s1 < 1L || s1 >= s2) stop("`s1` must satisfy 1 <= s1 < s2.", call. = FALSE)
  if (s2 > n)  stop("`s2` must satisfy s2 <= n.", call. = FALSE)
  if (dim < 1L) stop("`dim` must be a positive integer.", call. = FALSE)
  normalize <- .validate_scalar_logical(normalize, "normalize")

  w_rank <- tdnn_rank_table_asymp(n, s1, s2, dim, normalize = normalize)
  w <- w_rank[ranks]

  list(
    w      = w,
    w_rank = w_rank,
    meta   = list(
      n         = n,
      s1        = s1,
      s2        = s2,
      dim       = dim,
      normalize = normalize,
      scalars   = .tdnn_scalars(s1, s2, dim)
    )
  )
}
