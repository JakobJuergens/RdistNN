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
  if (s == 1L) {
    w_rank <- rep(1 / n, n)
  } else {
    j <- seq_len(n)
    w_rank <- .beta_interval_mass((j - 1L) / n, j / n, 1L, s)
    w_rank <- as.numeric(w_rank)
  }
  if (isTRUE(normalize)) {
    sw <- .compensated_sum(w_rank)
    if (!is.finite(sw) || sw <= 0) {
      .dnn_abort("numerical_range",
        "Asymptotic DNN weights have a nonpositive or nonfinite computed total.",
        hint = "Use less extreme feasible orders, or report a reproducible example if this persists."
      )
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
#' @seealso `tdnn_rank_table_exact()`, [rank_weight_tdnn_asymp()],
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

  if (n < 1L) .dnn_abort("argument_invalid", "`n` must be a positive integer.", argument = "n")
  if (s1 < 1L || s1 >= s2) .dnn_abort("argument_invalid", "`s1` must satisfy 1 <= s1 < s2.", argument = "s1")
  if (s2 > n)  .dnn_abort("argument_invalid", "`s2` must satisfy s2 <= n.", argument = "s2")
  if (dim < 1L) .dnn_abort("argument_invalid", "`dim` must be a positive integer.", argument = "dim")
  normalize <- .validate_scalar_logical(normalize, "normalize")

  sc <- .tdnn_scalars(s1, s2, dim)
  t1 <- .dnn_rank_table_asymp_k1(n, s1, normalize = FALSE)
  t2 <- .dnn_rank_table_asymp_k1(n, s2, normalize = FALSE)
  w_rank <- sc["w1"] * t1 + sc["w2"] * t2

  if (normalize) {
    sw <- .compensated_sum(w_rank)
    if (!is.finite(sw) || abs(sw) < .Machine$double.eps * 100) {
      .dnn_abort("numerical_range",
        "Asymptotic TDNN weights have a zero or nonfinite computed total.",
        hint = "Use less extreme feasible orders, or report a reproducible example if this persists."
      )
    }
    if (abs(sw - 1) > 1e-10) w_rank <- w_rank / sw
  }
  w_rank
}

# ---- Public: computation function (mirrors dnn_weights_asymp_full_from_ranks) ----
