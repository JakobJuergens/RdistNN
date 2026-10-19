# Asymptotic DNN L-statistic weights (FULL U-statistic, WOR) from ranks
#
# Assumes: averaging over ALL subsets of size s without replacement.
# Input :  ranks (length n), where ranks[i] is the global rank (1..n) of obs i.
# Output:  weights in data order and in rank order, using the asymptotic L-stat:
#
#   For rank j, approximate weight
#     w_rank[j] ≈ (1/k) * Σ_{i=1}^k [ F_{Beta(i, s-i+1)}(j/n) - F_{Beta(i, s-i+1)}((j-1)/n) ].
#
# Justification: when n is large, the ordered positions of an s-subset sampled
#               without replacement from {1,..,n} behave like the order stats of
#               s i.i.d. Uniform(0,1), whose i-th order stat has Beta(i, s-i+1).
#               See Steele (2009) and your conditional Z estimation paper.
# -------------------------------------------------------------------

#' Asymptotic DNN weights from ranks (full U-statistic, WOR)
#'
#' Compute **asymptotic** DNN (k-NN) weights via the L-statistic representation
#' over the full U-statistic of order \code{s}, using Beta order-statistic
#' increments as an approximation to the exact finite-population law.
#'
#' Unlike the exact finite-population rule, this asymptotic approximation is
#' evaluated on the full rank grid `1:n` and is therefore typically treated as a
#' full-support rank-weight sequence.
#'
#' @param ranks Integer permutation of \code{1:n}; \code{ranks[i]} is the global rank of obs \code{i}.
#' @param k Integer, number of nearest neighbors used inside each s-subset; must satisfy \code{1 <= k <= s}.
#' @param s Integer, U-statistic order (subset size); must satisfy \code{1 <= s <= n}.
#' @param normalize Logical; if TRUE (default) renormalizes weights to sum exactly 1 (guards roundoff).
#'
#' @return A list with
#' \itemize{
#'   \item \code{w}: numeric length-\code{n} weights in data order (sum to 1).
#'   \item \code{w_rank}: numeric length-\code{n} weights in rank order (1..n).
#'   \item \code{meta}: list with \code{n}, \code{k}, \code{s}, \code{approx = "beta(OrderStat)"}.
#' }
#' @examples
#' set.seed(1)
#' n <- 50
#' s <- 10
#' k <- 3
#' ranks <- sample.int(n)
#' out <- dnn_weights_asymp_full_from_ranks(ranks, k = k, s = s)
#' sum(out$w) # 1
#' sum(out$w_rank) # 1
#' @keywords internal
#' @noRd
dnn_weights_asymp_full_from_ranks <- function(ranks, k, s, normalize = TRUE) {
  # ---- validate ranks ----
  if (is.null(ranks)) stop("`ranks` must be provided.", call. = FALSE)
  if (!is.integer(ranks)) ranks <- as.integer(ranks)
  n <- length(ranks)
  if (n < 1L) stop("`ranks` must have positive length.", call. = FALSE)
  if (anyNA(ranks)) stop("`ranks` cannot contain NA.", call. = FALSE)
  if (!setequal(ranks, seq_len(n))) stop("`ranks` must be a permutation of 1..n.", call. = FALSE)

  # ---- validate k, s ----
  k <- .validate_scalar_int(k, "k")
  s <- .validate_scalar_int(s, "s")
  if (s > n) .abort_bad_arg("s", sprintf("must satisfy s <= n (got s=%d, n=%d)", s, n))
  if (k > s) .abort_bad_arg("k", sprintf("must satisfy k <= s (got k=%d, s=%d)", k, s))
  normalize <- .validate_scalar_logical(normalize, "normalize")

  # ---- asymptotic Beta-increment weights over the j/n grid ----
  j <- seq_len(n)
  ii <- seq_len(k)
  # n x k matrix of CDF increments on intervals ((j-1)/n, j/n]
  inc <- outer(
    j, ii,
    function(J, I) {
      .beta_interval_mass((J - 1) / n, J / n, I, s - I + 1)
    }
  )
  w_rank <- as.numeric(rowMeans(inc))

  # numerical guard: renormalize if requested
  if (normalize) {
    sw <- .compensated_sum(w_rank)
    if (!is.finite(sw) || sw <= 0) stop("Internal error: nonpositive total weight.", call. = FALSE)
    if (abs(sw - 1) > 1e-12) w_rank <- w_rank / sw
  }

  # map back to data order
  w <- w_rank[ranks]

  list(
    w = w,
    w_rank = w_rank,
    meta = list(
      n = n, k = as.integer(k), s = as.integer(s),
      approx = "beta(OrderStat)"
    )
  )
}
