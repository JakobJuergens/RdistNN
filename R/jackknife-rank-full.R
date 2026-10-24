# Generic rank-order delete-one jackknife via prefix/suffix sums.
#
# jackknife_rank_full() is the public entry point: it dispatches to the
# compiled C++ engine (rustat_jackknife_rank_full) when available, and falls
# back to the pure-R implementation (.jackknife_rank_full_R) otherwise.
#
# This function is intentionally rank-weight agnostic: the caller supplies
# precomputed rank-weight tables wn (at sample size n) and wnm1 (at n-1).
# The TDNN-specific wrapper that prepares those tables lives in
# user-tdnn-core.R.  The same engine can be reused for single-scale DNN or
# any other rule whose rank-weight sequence can be precomputed.

# ---- Pure-R fallback ---------------------------------------------------------

.jackknife_rank_full_R <- function(ranks, PHI, wn, wnm1) {
  # Implements the same O(n * m * q) prefix/suffix algorithm as the C++ engine.
  #
  # For each target t and signal dimension qq:
  #   1. Build signal in rank order: z[r] = PHI[ord[r], qq].
  #   2. Forward prefix:   P[r] = sum_{j=1}^{r-1} wnm1[j] * z[j]
  #   3. Backward suffix:  S[r] = sum_{j=r+1}^{n} wnm1[j-1] * z[j]
  #   4. theta_loo[i, qq, t] = P[ranks[i,t]] + S[ranks[i,t]]
  #
  # All index arithmetic uses 1-based ranks (R convention).

  n   <- nrow(ranks)
  m   <- ncol(ranks)
  q   <- ncol(PHI)
  table_mass <- .compensated_sum_R(wnm1)
  weight_amplitude <- if (length(wnm1)) max(abs(wnm1)) else 0

  theta0  <- matrix(NA_real_, q, m)
  loo_arr <- array(NA_real_, dim = c(n, q, m))

  for (t in seq_len(m)) {
    r_col <- ranks[, t]          # 1-based rank vector (length n)
    ord   <- order(r_col)        # ord[r]: 1-based obs index at rank r

    # Baseline: theta0[qq, t] = sum_i wn[r_i] * PHI[i, qq]
    for (qq in seq_len(q))
      theta0[qq, t] <- .compensated_sum_R(wn[r_col] * PHI[, qq])

    for (qq in seq_len(q)) {
      z <- PHI[ord, qq]          # signal in rank order (length n)
      active <- .jackknife_active_ranks(wnm1)
      level <- .numeric_signal_level(matrix(z[active], ncol = 1L),
        weights = weight_amplitude, weight_mass = table_mass)
      z[!active] <- 0
      z[active] <- z[active] - level

      value <- .compensated_deletions_R(z, wnm1,
        level * table_mass)
      loo_arr[, qq, t] <- value[r_col]
    }
  }

  list(
    theta0    = theta0,
    theta_loo = loo_arr,
    info      = list(n = n, m = m, q = q)
  )
}

# ---- Public entry point ------------------------------------------------------

#' Generic rank-order delete-one jackknife (prefix/suffix engine)
#'
#' Computes all \eqn{n} leave-one-out estimates for a rank-weighted linear
#' estimator \eqn{\hat\theta_t = \sum_i w_n[r_{it}]\,\phi_i} at each of
#' \eqn{m} target points, using an \eqn{O(n \cdot m \cdot q)} prefix/suffix
#' algorithm.
#'
#' The rank-weight rule is supplied via precomputed tables `wn` and `wnm1`
#' rather than a spec object, so the function is completely agnostic to the
#' weight construction.  Callers are responsible for providing consistent
#' tables:
#' \itemize{
#'   \item `wn` should be the rank-weight sequence evaluated at sample size
#'     \eqn{n} (i.e. for the full-sample estimator).
#'   \item `wnm1` should be the **same rule** evaluated at sample size
#'     \eqn{n-1} (i.e. for each leave-one-out sub-sample).
#' }
#'
#' With compiled execution enabled and the native routine available, this uses
#' the C++ engine. Otherwise it uses the same algorithm in R. See [RdistNN] for
#' engine controls. Native execution errors propagate without an automatic retry.
#'
#' @param ranks Integer matrix (\eqn{n \times m}); column \eqn{t} is a
#'   permutation of \code{1:n} giving the rank of each observation for
#'   target \eqn{t} (rank 1 = nearest).
#' @param PHI Finite numeric matrix (\eqn{n \times q}); row \eqn{i} is the
#'   \eqn{q}-dimensional signal contribution of observation \eqn{i}.
#'   For a scalar signal, pass \code{cbind(y)}.
#' @param wn Finite numeric vector of length \eqn{n}; rank-weight table at full
#'   sample size \eqn{n}.  Entry \code{wn[r]} is the weight assigned to
#'   rank \eqn{r}.
#' @param wnm1 Finite numeric vector of length \eqn{n-1}; rank-weight table at
#'   leave-one-out sample size \eqn{n-1}.
#'
#' @return A list with:
#' \describe{
#'   \item{`theta0`}{Numeric matrix (\eqn{q \times m}); baseline estimates
#'     at full sample size.}
#'   \item{`theta_loo`}{Array with \code{dim = c(n, q, m)}; element
#'     \code{[i, qq, t]} is the leave-one-out estimate after deleting
#'     observation \eqn{i} for signal dimension \eqn{qq} and target \eqn{t}.}
#'   \item{`info`}{List with fields \code{n}, \code{m}, \code{q}.}
#' }
#'
#' @seealso [jackknife_fast_rank()] for the finite-support L-stat engine;
#'   [jackknife_slow_rank()] for the reference implementation.
#'
#' @examples
#' set.seed(1)
#' n <- 20L; m <- 3L; q <- 2L
#' ranks <- apply(matrix(runif(n * m), n, m), 2L, rank)
#' storage.mode(ranks) <- "integer"
#' PHI <- matrix(rnorm(n * q), n, q)
#' # toy uniform rank weights
#' wn   <- rep(1 / n, n)
#' wnm1 <- rep(1 / (n - 1L), n - 1L)
#' out <- jackknife_rank_full(ranks, PHI, wn, wnm1)
#' dim(out$theta_loo)   # c(n, q, m)
#' dim(out$theta0)      # c(q, m)
#' @export
jackknife_rank_full <- function(ranks, PHI, wn, wnm1) {
  ranks <- .validate_permutation_matrix(ranks, arg = "ranks", col_arg = "ranks")
  PHI <- .validate_matrix(
    PHI,
    arg = "PHI",
    min_rows = 1L,
    min_cols = 1L,
    numeric = TRUE,
    finite = TRUE,
    strict = TRUE
  )
  wn <- .validate_finite_numeric_vector(wn, arg = "wn", min_length = 1L)
  wnm1 <- .validate_finite_numeric_vector(wnm1, arg = "wnm1", min_length = 0L)

  n <- nrow(ranks)
  if (length(wn) != n)
    .dnn_abort("argument_invalid", "`wn` must have length nrow(ranks).", argument = "wn")
  if (length(wnm1) != n - 1L)
    .dnn_abort("argument_invalid", "`wnm1` must have length nrow(ranks) - 1.", argument = "wnm1")
  if (nrow(PHI) != n)
    .dnn_abort("argument_invalid", "`PHI` must have nrow(ranks) rows.", argument = "PHI")

  cpp <- .rdistnn_compiled_fun("rustat_jackknife_rank_full")
  if (!is.null(cpp)) {
    return(cpp(ranks, PHI, wn, wnm1))
  }

  .jackknife_rank_full_R(ranks, PHI, wn, wnm1)
}
