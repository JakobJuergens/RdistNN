#' Slow delete-one jackknife for rank-weighted L-stat rules (fixed ranks)
#'
#' A reference implementation that recomputes *aggregates* per deletion
#' using the provided neighbor ranks. This **does not** recompute ranks
#' after deletion (useful as a correctness oracle and for testing).
#'
#' Like [jackknife_fast_rank()], this helper targets compact-support
#' finite-rank L-stat rules. It is not the internal engine used by the
#' user-level exact DNN wrappers, which instead update the full rank matrix
#' under deletion and reevaluate the DNN rank-weight rule.
#'
#' Works with:
#' - `phi` as a numeric **matrix** (n_ref x q): row `j` equals \eqn{\phi(j)}.
#' - `phi` as a **function** `function(j, t)` that returns a length-`q` numeric,
#'   allowing target-dependent contributions.
#'
#' @param ranks Integer matrix (T x (s+1)): for each target `t`, `ranks[t, k]`
#'   gives the 1-based index of the k-th neighbor, for `k = 1, ..., s+1`.
#'   (Must exclude `t` itself if targets are a subset of reference.)
#' @param a Numeric length-`s` vector of L-stat coefficients \eqn{(a_1,\dots,a_s)}.
#' @param phi Either a numeric matrix (n_ref x q) with row `j` equal to \eqn{\phi(j)},
#'   or a function `phi(j, t)` returning a length-`q` numeric.
#' @param normalize Logical; if `TRUE` (default) divide aggregates by `sum(a)`.
#'
#' @return A list with:
#' \describe{
#'   \item{theta0}{numeric length-`q` baseline aggregate}
#'   \item{theta_loo}{numeric matrix (n_ref x q) of delete-one aggregates}
#'   \item{info}{list with fields `s`, `T`, `q`, `normalize`}
#' }
#' @examples
#' set.seed(1)
#' n <- 20
#' p <- 3
#' s <- 4
#' X <- matrix(rnorm(n * p), n, p)
#' # toy ranks: nearest-neighbor indices precomputed elsewhere (T x (s+1))
#' # here we fake something deterministic just for the example
#' ranks <- cbind(matrix(rep(1:s, each = n), n, s), s + 1L)
#' y <- rnorm(n)
#' a <- rep(1, s)
#' out <- jackknife_slow_rank(ranks, a, cbind(y))
#' str(out)
#' @export
jackknife_slow_rank <- function(ranks, a, phi, normalize = TRUE) {
  ranks <- as.matrix(ranks)
  a <- as.numeric(a)

  Tn <- nrow(ranks)
  K <- ncol(ranks)
  s <- length(a)
  if (K < s + 1L) stop("ncol(ranks) must be at least s+1.")

  # Build per-(t,k) contribution array y_tk: (T x (s+1) x q)
  if (is.matrix(phi) && is.numeric(phi)) {
    PHI <- as.matrix(phi)
    n_ref <- nrow(PHI)
    q <- ncol(PHI)
    y_tk <- array(0.0, dim = c(Tn, s + 1L, q))
    for (k in seq_len(s + 1L)) {
      idx <- ranks[, k]
      if (any(idx < 1L | idx > n_ref)) stop("ranks indices out of PHI range.")
      y_tk[, k, ] <- PHI[idx, , drop = FALSE]
    }
  } else if (is.function(phi)) {
    n_ref <- max(ranks)
    # probe q from first call
    q_probe <- phi(ranks[1L, 1L], 1L)
    q <- length(q_probe)
    y_tk <- array(0.0, dim = c(Tn, s + 1L, q))
    for (t in seq_len(Tn)) {
      for (k in seq_len(s + 1L)) {
        j <- ranks[t, k]
        y_tk[t, k, ] <- phi(j, t)
      }
    }
  } else {
    stop("`phi` must be either a numeric matrix (n_ref x q) or a function(j, t).")
  }

  # Normalization constant
  Z <- if (isTRUE(normalize)) sum(a) else 1
  if (isTRUE(normalize) && Z == 0) stop("sum(a) == 0 with normalize=TRUE.")

  # Baseline per-target sums: S_t/Z where S_t = sum_{k<=s} a_k * y_{t,k}
  S_t_divZ <- matrix(0.0, Tn, q)
  for (k in seq_len(s)) {
    # y_tk[, k, ] drops to (Tn × q) if q>1, but to length-Tn vector if q==1.
    # Wrap in matrix(...) to guarantee (Tn × q).
    S_t_divZ <- S_t_divZ + a[k] * matrix(y_tk[, k, ], Tn, q)
  }
  S_t_divZ <- S_t_divZ / Z
  theta0 <- colMeans(S_t_divZ)

  # Initialize theta_loo to theta0
  theta_loo <- matrix(rep(theta0, each = n_ref), n_ref, q, byrow = FALSE)

  # Slow (reference) delete-one loop:
  # For each i, scan all targets, check position r where i appears, and
  # rebuild the shifted sum with neighbor s+1 pulled in.
  for (i in seq_len(n_ref)) {
    delta <- numeric(q)
    for (t in seq_len(Tn)) {
      # look only within first s+1 neighbors
      pos <- which(ranks[t, 1L:(s + 1L)] == i)
      if (length(pos) == 0L) next
      r <- pos[1L] # ensure scalar
      if (r > s) next # present only at s+1 -> no effect

      base <- S_t_divZ[t, ]
      term1 <- (a[r] / Z) * y_tk[t, r, ]

      # shift block sum (scale by 1/Z to match base scale)
      term2 <- if (s >= 2L && r <= s - 1L) {
        acc <- rep(0.0, q)
        for (k in r:(s - 1L)) acc <- acc - (a[k + 1L] - a[k]) * y_tk[t, k + 1L, ]
        acc / Z
      } else {
        rep(0.0, q)
      }

      term3 <- (a[s] / Z) * y_tk[t, s + 1L, ]
      Snew <- base - term1 + term2 + term3
      delta <- delta + (Snew - base)
    }
    theta_loo[i, ] <- theta0 + delta / Tn
  }

  list(
    theta0    = as.numeric(theta0),
    theta_loo = theta_loo,
    info      = list(s = s, T = Tn, q = q, normalize = normalize)
  )
}
