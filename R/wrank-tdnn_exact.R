# Exact TDNN rank-weight table (Demirkaya two-scale DNN, k = 1)
#
# TDNN is a signed linear combination of two single-scale exact DNN
# rank-weight tables:
#
#   w_rank^TDNN = w1* * w_rank^DNN(s1) + w2* * w_rank^DNN(s2)
#
# The bias-correction scalars (Demirkaya et al.) are:
#
#   w1* = 1 / (1 - (s1 / s2)^(-2 / d))
#   w2* = 1 - w1*
#
# where d is the localization dimension (typically ncol(X)).
# Since s1 < s2 and d > 0, we have (s1/s2)^(-2/d) > 1, so w1* < 0
# and w2* > 1.  TDNN weights are therefore signed; they sum to 1 but
# individual entries can be negative.
#
# For k = 1 the single-scale DNN rank coefficient at rank j is:
#
#   a_j = C(n - j, s - 1) / C(n, s)    for j = 1, ..., n - s + 1
#         0                              otherwise
#
# The TDNN support cutoff is n - s1 + 1 (not n - s2 + 1) because the
# smaller s1 piece can still assign mass to ranks that the larger s2
# piece cannot reach.
# --------------------------------------------------------------------------

# ---- Internal: single-scale rank table for k = 1 -------------------------

.dnn_rank_table_exact_k1 <- function(n, s) {
  .dnn_exact_table(as.integer(n), 1L, as.integer(s))
}

# ---- Bias-correction scalars ---------------------------------------------

.tdnn_scalars <- function(s1, s2, dim) {
  # Returns c(w1 = w1*, w2 = w2*) for the Demirkaya TDNN construction.
  # s1 < s2, dim > 0 (localization dimension).
  # Retain the small order difference before taking logarithms, then avoid
  # cancellation when the exponential is close to one.
  w1 <- -1 / expm1((2 / dim) * log1p((s2 - s1) / s1))
  c(w1 = w1, w2 = 1 - w1)         # w2 > 1, w1 + w2 = 1
}

# ---- Public: TDNN exact rank table ---------------------------------------

#' Exact TDNN rank-weight table (Demirkaya two-scale DNN, k = 1)
#'
#' Computes the rank-order weight sequence for the two-scale DNN (TDNN)
#' estimator of Demirkaya et al. as a signed linear combination of two
#' exact single-scale DNN rank-weight tables evaluated at `(k=1, s1)` and
#' `(k=1, s2)`:
#'
#' \deqn{w^{\text{TDNN}}_j = w_1^* \, w^{\text{DNN}}_j(s_1) + w_2^* \, w^{\text{DNN}}_j(s_2)}
#'
#' with bias-correction scalars
#' \deqn{w_1^* = \frac{1}{1 - (s_1/s_2)^{-2/d}}, \quad w_2^* = 1 - w_1^*,}
#' where \eqn{d} is the localization dimension.
#'
#' **Note:** TDNN weights are signed.  They sum to 1 but individual
#' entries can be negative (specifically, ranks near 1 receive negative
#' weight from the \eqn{s_1} term).  This is by design: the signed
#' combination cancels the leading bias term in the DNN MSE expansion.
#'
#' @param n Positive integer; sample size.
#' @param s1 Positive integer; smaller U-statistic order (\code{s1 < s2 <= n}).
#' @param s2 Positive integer; larger U-statistic order (\code{s1 < s2 <= n}).
#' @param dim Positive integer; localization dimension (typically \code{ncol(X)}).
#'
#' @return A numeric vector of length \code{n} giving the rank-order TDNN
#'   weights \code{w_rank[1:n]} (rank 1 = nearest).  Entries sum to 1 but
#'   some can be negative.
#'
#' @seealso `tdnn_rank_table_asymp()`, `dnn_weights_tdnn_exact_full_from_ranks()`,
#'   [rank_weight_tdnn_exact()]
#'
#' @examples
#' n <- 20; s1 <- 3L; s2 <- 6L; d <- 2L
#' w <- tdnn_rank_table_exact(n, s1, s2, d)
#' sum(w)          # 1
#' any(w < 0)      # TRUE (bias-correction weights)
#' @keywords internal
#' @noRd
tdnn_rank_table_exact <- function(n, s1, s2, dim) {
  n   <- as.integer(n)
  s1  <- as.integer(s1)
  s2  <- as.integer(s2)
  dim <- as.integer(dim)

  if (n < 1L) stop("`n` must be a positive integer.", call. = FALSE)
  if (s1 < 1L || s1 >= s2) stop("`s1` must satisfy 1 <= s1 < s2.", call. = FALSE)
  if (s2 > n)  stop("`s2` must satisfy s2 <= n.", call. = FALSE)
  if (dim < 1L) stop("`dim` must be a positive integer.", call. = FALSE)

  sc <- .tdnn_scalars(s1, s2, dim)
  t1 <- .dnn_rank_table_exact_k1(n, s1)
  t2 <- .dnn_rank_table_exact_k1(n, s2)
  w_rank <- sc["w1"] * t1 + sc["w2"] * t2

  # Guard: sum should equal 1 exactly (w1* + w2* = 1 and sum(t1) = sum(t2) = 1)
  sw <- .compensated_sum(w_rank)
  if (!is.finite(sw) || abs(sw) < .Machine$double.eps * 100) {
    stop("Internal error: TDNN rank table has zero/non-finite sum.", call. = FALSE)
  }
  if (abs(sw - 1) > 1e-10) w_rank <- w_rank / sw
  w_rank
}

# ---- Public: computation function (mirrors dnn_weights_exact_full_from_ranks) ----

#' Exact TDNN weights from ranks (Demirkaya two-scale DNN, k = 1)
#'
#' Given a rank permutation `ranks`, computes the TDNN weight for each
#' observation using the exact rank-weight table from `tdnn_rank_table_exact()`.
#'
#' @param ranks Integer permutation of \code{1:n}; \code{ranks[i]} is the
#'   rank (1 = nearest) of observation \code{i} relative to the target.
#' @param s1 Positive integer; smaller U-statistic order.
#' @param s2 Positive integer; larger U-statistic order.
#' @param dim Positive integer; localization dimension.
#'
#' @return A list with
#' \itemize{
#'   \item \code{w}: numeric length-\code{n} weights in data order (sum to 1,
#'     some entries may be negative).
#'   \item \code{w_rank}: numeric length-\code{n} weights in rank order (1..n).
#'   \item \code{meta}: list with \code{n}, \code{s1}, \code{s2}, \code{dim},
#'     and \code{scalars} (the bias-correction coefficients \eqn{w_1^*, w_2^*}).
#' }
#'
#' @seealso `tdnn_rank_table_exact()`, [rank_weight_tdnn_exact()],
#'   `dnn_weights_tdnn_asymp_full_from_ranks()`
#'
#' @examples
#' set.seed(1)
#' n <- 15L; ranks <- sample.int(n)
#' out <- dnn_weights_tdnn_exact_full_from_ranks(ranks, s1 = 3L, s2 = 6L, dim = 2L)
#' sum(out$w)        # 1
#' any(out$w < 0)    # TRUE
#' @keywords internal
#' @noRd
dnn_weights_tdnn_exact_full_from_ranks <- function(ranks, s1, s2, dim) {
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

  w_rank <- tdnn_rank_table_exact(n, s1, s2, dim)
  w <- w_rank[ranks]

  list(
    w      = w,
    w_rank = w_rank,
    meta   = list(
      n       = n,
      s1      = s1,
      s2      = s2,
      dim     = dim,
      scalars = .tdnn_scalars(s1, s2, dim)
    )
  )
}
