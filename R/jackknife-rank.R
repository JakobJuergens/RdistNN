#' Fast delete-one jackknife for rank-weighted functionals
#' (auto C++ when available, R fallback otherwise)
#'
#' Computes the baseline aggregate `theta0` and all delete-one aggregates
#' `theta_loo` for rank-weighted L-stat rules with finite support s.
#'
#' - If `phi` is a numeric **matrix** (n_ref x q) with rows `phi[j, ] = phi(j)` and
#'   the compiled symbol `rustat_jackknife_fast_rank()` is available, this calls
#'   the **C++ fast path**.
#' - Otherwise, it uses a general **R fallback** that also supports target-dependent
#'   per-neighbor contributions via a function `phi(j, t)`.
#'
#' @param ranks integer matrix (T x (s+1)) of neighbor indices (1-based).
#'              Each row t lists the first s+1 neighbors of target t.
#' @param a numeric length-s vector of L-stat coefficients a_1..a_s.
#' @param phi either:
#'   - numeric matrix (n_ref x q): row j is the q-vector phi(j), or
#'   - function `function(j, t)` returning a numeric length-q vector (allows
#'     target dependence).
#'   For a scalar \eqn{\phi}, use `cbind(y)` as the matrix case.
#' @param normalize logical; if TRUE (default) divide by sum(a).
#'
#' @return list with:
#'   - `theta0`: numeric length-q baseline aggregate,
#'   - `theta_loo`: numeric matrix (n_ref x q) of delete-one aggregates,
#'   - `info`: list(s, T, q, normalize).
#' @examples
#' set.seed(1)
#' n <- 20
#' p <- 2
#' X <- matrix(rnorm(n * p), n, p)
#' ranks <- cbind(matrix(sample.int(n, n * 3, replace = TRUE), n, 3), sample.int(n, n, replace = TRUE))
#' a <- c(1, 1, 1)
#' phi <- cbind(rnorm(n))
#' out <- jackknife_fast_rank(ranks, a, phi, normalize = TRUE)
#' str(out$theta0)
#' @export
jackknife_fast_rank <- function(ranks, a, phi, normalize = TRUE) {
  ranks <- as.matrix(ranks)
  a <- as.numeric(a)

  # Prefer compiled path when phi is matrix and the symbol is present
  if (is.matrix(phi) && is.numeric(phi) &&
    exists("rustat_jackknife_fast_rank", mode = "function")) {
    PHI <- as.matrix(phi)
    return(rustat_jackknife_fast_rank(ranks, a, PHI, normalize))
  }

  # Otherwise fall back to the general R implementation
  .jackknife_fast_rank_R(ranks = ranks, a = a, phi = phi, normalize = normalize)
}

# --- Internal R engine (supports phi as matrix or function(j, t)) ---------

.jackknife_fast_rank_R <- function(ranks, a, phi, normalize = TRUE) {
  stopifnot(is.matrix(ranks))
  Tn <- nrow(ranks)
  K <- ncol(ranks)
  s <- length(a)
  if (K < s + 1L) stop("ncol(ranks) must be at least s+1.")

  # Normalization constant
  Z <- if (isTRUE(normalize)) sum(a) else 1
  if (isTRUE(normalize) && Z == 0) stop("sum(a) == 0 with normalize=TRUE.")

  # Prepare y_tk: array (T x (s+1) x q) of per-(t,k) contributions
  if (is.matrix(phi) && is.numeric(phi)) {
    PHI <- as.matrix(phi)
    n_ref <- nrow(PHI)
    q <- ncol(PHI)
    y_tk <- array(0.0, dim = c(Tn, s + 1L, q))
    for (k in seq_len(s + 1L)) {
      idx <- ranks[, k]
      if (any(idx < 1L) || any(idx > n_ref)) stop("ranks indices out of PHI range.")
      y_tk[, k, ] <- PHI[idx, , drop = FALSE]
    }
  } else if (is.function(phi)) {
    n_ref <- max(ranks)
    # Probe q from the first call
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

  # Baseline per-target sums: S_t = sum_{k<=s} a_k * y_{t,k}
  S_t_divZ <- matrix(0.0, Tn, q)
  for (k in seq_len(s)) {
    S_t_divZ <- S_t_divZ + a[k] * matrix(y_tk[, k, ], Tn, q)
  }
  S_t_divZ <- S_t_divZ / Z

  # theta0 = mean_t S_t_divZ[t,]
  theta0 <- colMeans(S_t_divZ)

  # pre_t[t, r] = sum_{k=r}^{s-1} (a_{k+1}-a_k) * y_{t, k+1}, r=1..s-1
  if (s >= 2L) {
    pre_t <- array(0.0, dim = c(Tn, s - 1L, q))
    for (t in seq_len(Tn)) {
      acc <- rep(0.0, q)
      for (r in (s - 1L):1L) {
        acc <- acc + (a[r] - a[r + 1L]) * y_tk[t, r + 1L, ]
        pre_t[t, r, ] <- acc
      }
    }
  } else {
    pre_t <- array(0.0, dim = c(Tn, 0L, q)) # unused when s=1
  }

  # Build adjacency: for each i, list of (t, r) where i is rank r (1..s+1)
  n_ref_adj <- if (exists("PHI")) nrow(PHI) else max(ranks)
  affected <- vector("list", length = n_ref_adj)
  for (t in seq_len(Tn)) {
    for (k in seq_len(s + 1L)) {
      i <- ranks[t, k]
      if (is.null(affected[[i]])) {
        affected[[i]] <- cbind(integer(0L), integer(0L))
      }
      affected[[i]] <- rbind(affected[[i]], c(t, k))
    }
  }

  # Initialize theta_loo to theta0
  theta_loo <- matrix(rep(theta0, each = n_ref_adj), n_ref_adj, q, byrow = FALSE)

  # Per-i updates: only r <= s affect the sum (support)
  for (i in seq_len(n_ref_adj)) {
    lst <- affected[[i]]
    if (is.null(lst) || nrow(lst) == 0L) next
    delta <- rep(0.0, q)
    for (row in seq_len(nrow(lst))) {
      t <- lst[row, 1L]
      r <- lst[row, 2L] # already scalar 1..s+1 here
      if (r > s) next

      base <- S_t_divZ[t, ]
      term1 <- (a[r] / Z) * y_tk[t, r, ]
      term2 <- if (s >= 2L && r <= s - 1L) pre_t[t, r, ] / Z else rep(0.0, q)
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
