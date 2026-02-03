# Exact DNN L-statistic weights (FULL U-statistic, WOR) from ranks
#
# Assumes: we average over ALL subsets of size s without replacement.
# Input:   ranks (length n), where ranks[i] is the rank (1..n) of obs i.
# Output:  weights in data order, computed via the exact L-statistic law:
#
#   For rank j, the exact weight is
#     w_rank[j] = (1/k) * sum_{i=1}^k P{ the i-th NN in a random s-subset is the j-th global NN }
#   with
#     P = C(j-1, i-1) * C(n-j, s-i) / C(n, s),
#   valid when 1 <= i <= min(k, s), i <= j, and s - i <= n - j; else 0.
#
# The returned weights sum to 1 (up to numerical tolerance).
# -------------------------------------------------------------------

#' Exact DNN weights from ranks (full U-statistic, WOR)
#'
#' Compute exact DNN (k-NN) weights for a target point using the L-statistic
#' representation over the **full** U-statistic of order \code{s} (all \eqn{\binom{n}{s}}
#' subsets sampled \emph{without} replacement). The function is always given the
#' permutation of ranks \code{ranks}, where \code{ranks[i]} is the global rank
#' (1 = nearest) of observation \code{i} relative to the target.
#'
#' @param ranks integer permutation of \code{1:n}; \code{ranks[i]} is the rank of obs \code{i}.
#' @param k integer, number of nearest neighbors used inside each s-subset; must satisfy \code{1 <= k <= s}.
#' @param s integer, U-statistic order (subset size); must satisfy \code{1 <= s <= n}.
#' @return A list with
#' \itemize{
#'   \item \code{w}: numeric length-\code{n} weights in data order (sum to 1).
#'   \item \code{w_rank}: numeric length-\code{n} weights in rank order (1..n).
#'   \item \code{meta}: list with \code{n}, \code{k}, \code{s}.
#' }
#' @examples
#' set.seed(1)
#' n <- 10
#' s <- 4
#' k <- 2
#' ranks <- sample.int(n) # ranks[i] in 1..n
#' out <- dnn_weights_exact_full_from_ranks(ranks, k = k, s = s)
#' sum(out$w) # 1
#' sum(out$w_rank) # 1
#' @export
dnn_weights_exact_full_from_ranks <- function(ranks, k, s) {
  # ---- validate ranks ----
  if (is.null(ranks)) stop("`ranks` must be provided.")
  if (!is.integer(ranks)) ranks <- as.integer(ranks)
  n <- length(ranks)
  if (n < 1L) stop("`ranks` must have positive length.")
  if (anyNA(ranks)) stop("`ranks` cannot contain NA.")
  if (!setequal(ranks, seq_len(n))) stop("`ranks` must be a permutation of 1..n.")

  # ---- validate k, s ----
  if (!is.numeric(s) || length(s) != 1L || !is.finite(s) || s < 1L || s > n || as.integer(s) != s) {
    stop("`s` must be an integer with 1 <= s <= n.")
  }
  if (!is.numeric(k) || length(k) != 1L || !is.finite(k) || k < 1L || as.integer(k) != k) {
    stop("`k` must be a positive integer.")
  }
  if (k > s) stop("`k` must satisfy k <= s in the full U-statistic setting.")

  # ---- exact rank-weights via hypergeometric law ----
  # w_rank[j] = (1/k) * sum_{i=1}^k C(j-1, i-1) * C(n-j, s-i) / C(n, s)
  # with the convention that out-of-range binomials contribute 0.
  L_den <- lchoose(n, s)
  w_rank <- numeric(n)

  for (j in seq_len(n)) {
    # i must satisfy: 1 <= i <= k, i <= j, and s - i <= n - j  <=>  i >= s - (n - j)
    i_lo <- max(1L, s - (n - j))
    i_hi <- min(k, j, s)
    if (i_lo <= i_hi) {
      acc <- 0.0
      for (i in i_lo:i_hi) {
        # log-probability; exp() is safe for moderate n and exact enough in double
        lp <- lchoose(j - 1L, i - 1L) + lchoose(n - j, s - i) - L_den
        acc <- acc + exp(lp)
      }
      w_rank[j] <- acc / k
    } # else remains 0
  }

  # normalize (guard against tiny numerical drift)
  sw <- sum(w_rank)
  if (!is.finite(sw) || sw <= 0) stop("Internal error: nonpositive total weight.")
  if (abs(sw - 1) > 1e-12) w_rank <- w_rank / sw

  # map back to data order: w[i] = w_rank[ ranks[i] ]
  w <- w_rank[ranks]

  list(
    w      = w,
    w_rank = w_rank,
    meta   = list(n = n, k = as.integer(k), s = as.integer(s))
  )
}
