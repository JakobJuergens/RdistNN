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
#' For the exact rule, weights are computed from the full rank permutation, not
#' by truncating to the first `s` observations. The resulting exact support in
#' rank order is:
#'
#' - `j <= n - s + k` in the general `k`-NN-within-subset case;
#' - in particular, `j <= n - s + 1` when `k = 1`.
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
#' @keywords internal
#' @noRd
dnn_weights_exact_full_from_ranks <- function(ranks, k, s) {
  # ---- validate ranks ----
  if (is.null(ranks)) .dnn_abort("argument_invalid", "`ranks` must be provided.", argument = "ranks")
  if (!is.integer(ranks)) ranks <- as.integer(ranks)
  n <- length(ranks)
  if (n < 1L) .dnn_abort("argument_invalid", "`ranks` must have positive length.", argument = "ranks")
  if (anyNA(ranks)) .dnn_abort("argument_invalid", "`ranks` cannot contain NA.", argument = "ranks")
  if (!setequal(ranks, seq_len(n))) .dnn_abort(
    "argument_invalid",
    "`ranks` must be a permutation of 1..n.",
    argument = "ranks"
  )

  # ---- validate k, s ----
  k <- .validate_scalar_int(k, "k")
  s <- .validate_scalar_int(s, "s")
  if (s > n) .abort_bad_arg("s", sprintf("must satisfy s <= n (got s=%d, n=%d)", s, n))
  if (k > s) .abort_bad_arg("k", sprintf("must satisfy k <= s (got k=%d, s=%d)", k, s))

  w_rank <- .dnn_exact_table(n, k, s)

  # map back to data order: w[i] = w_rank[ ranks[i] ]
  w <- w_rank[ranks]

  list(
    w      = w,
    w_rank = w_rank,
    meta   = list(n = n, k = as.integer(k), s = as.integer(s))
  )
}
