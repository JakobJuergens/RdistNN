#' Fast delete-one jackknife for rank-weighted functionals
#' (auto C++ when available, R fallback otherwise)
#'
#' Computes the baseline aggregate `theta0` and all delete-one aggregates
#' `theta_loo` for generic rank-weighted L-stat rules with finite support `s`.
#'
#' This utility is for rules that are fully characterized by:
#'
#' - the first `s + 1` ranked observations for each target, and
#' - a finite coefficient vector `a_1, ..., a_s`.
#'
#' It is therefore appropriate for compact-support L-statistics. The user-level
#' DNN wrappers (`np_dnn()`, `cate_dnn()`, `plqr_dnn()`) do **not** call this
#' helper directly for exact DNN weights. Those wrappers keep the full rank
#' matrix and reevaluate the DNN rank-weight rule after delete-one updates.
#'
#' - If `phi` is a numeric **matrix** (n_ref x q) with rows `phi[j, ] = phi(j)` and
#'   the compiled symbol `rustat_jackknife_fast_rank()` is available, this calls
#'   the **C++ fast path**.
#' - Otherwise, it uses a general **R fallback** that also supports target-dependent
#'   per-neighbor contributions via a function `phi(j, t)`.
#'
#' @param ranks Integer matrix (T x (s+1)) of positive neighbor indices
#'   (1-based). Each row t lists the first s+1 neighbors of target t.
#' @param a Finite numeric length-s vector of L-stat coefficients a_1..a_s.
#' @param phi either:
#'   - finite numeric matrix (n_ref x q): row j is the q-vector phi(j), or
#'   - function `function(j, t)` returning a finite numeric length-q vector
#'     (allows target dependence).
#'   For a scalar \eqn{\phi}, use `cbind(y)` as the matrix case.
#' @param normalize Scalar logical; if TRUE (default) divide by sum(a).
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
  args <- .jackknife_validate_finite_support_args(ranks, a, normalize)
  ranks <- args$ranks
  a <- args$a
  normalize <- args$normalize

  # Prefer compiled path when phi is matrix and the symbol is present
  if (is.matrix(phi) && is.numeric(phi)) {
    PHI <- .jackknife_validate_phi_matrix(phi, ranks, arg = "phi")
    if (exists("rustat_jackknife_fast_rank", mode = "function")) {
      return(rustat_jackknife_fast_rank(ranks, a, PHI, normalize))
    }
    return(.jackknife_fast_rank_R(ranks = ranks, a = a, phi = PHI, normalize = normalize))
  }

  # Otherwise fall back to the general R implementation
  .jackknife_fast_rank_R(ranks = ranks, a = a, phi = phi, normalize = normalize)
}

.jackknife_validate_finite_support_args <- function(ranks, a, normalize) {
  ranks <- .validate_positive_index_matrix(ranks, arg = "ranks")
  a <- .validate_finite_numeric_vector(a, arg = "a", min_length = 1L)
  normalize <- .validate_scalar_logical(normalize, "normalize")
  if (ncol(ranks) < length(a) + 1L) {
    .abort_bad_arg(
      "ncol(ranks)",
      "must be at least s+1",
      expected = sprintf("at least %d columns because length(a)=%d", length(a) + 1L, length(a)),
      actual = sprintf("%d columns", ncol(ranks)),
      quoted = FALSE
    )
  }
  if (normalize && sum(a) == 0) {
    .abort_bad_arg("a", "must have nonzero sum when `normalize = TRUE`", actual = "sum(a)=0")
  }
  list(ranks = ranks, a = a, normalize = normalize)
}

.jackknife_validate_phi_matrix <- function(phi, ranks, arg = "phi") {
  PHI <- .validate_matrix(
    phi,
    arg = arg,
    min_rows = 1L,
    min_cols = 1L,
    numeric = TRUE,
    finite = TRUE,
    strict = TRUE
  )
  max_rank <- max(ranks)
  if (max_rank > nrow(PHI)) {
    .abort_bad_arg(
      "ranks",
      "indices out of PHI range",
      expected = sprintf("indices <= nrow(%s)=%d", arg, nrow(PHI)),
      actual = sprintf("max(ranks)=%d", max_rank)
    )
  }
  PHI
}

.jackknife_eval_phi <- function(phi, j, t, q = NULL) {
  value <- phi(j, t)
  if (!is.numeric(value) || length(value) == 0L) {
    .abort_bad_arg("phi(j, t)", "must return a finite numeric vector", quoted = FALSE)
  }
  value <- as.numeric(value)
  if (!is.null(q) && length(value) != q) {
    .abort_bad_arg(
      "phi(j, t)",
      "must return a consistent length",
      expected = sprintf("length %d", q),
      actual = sprintf("length %d", length(value)),
      quoted = FALSE
    )
  }
  .abort_if_na_or_nonfinite(value, arg = "phi(j, t)", quoted = FALSE)
  value
}

# --- Internal R engine (supports phi as matrix or function(j, t)) ---------

.jackknife_fast_rank_R <- function(ranks, a, phi, normalize = TRUE) {
  args <- .jackknife_validate_finite_support_args(ranks, a, normalize)
  ranks <- args$ranks
  a <- args$a
  normalize <- args$normalize
  Tn <- nrow(ranks)
  s <- length(a)

  # Normalization constant
  Z <- if (normalize) sum(a) else 1

  # Prepare y_tk: array (T x (s+1) x q) of per-(t,k) contributions
  if (is.matrix(phi) && is.numeric(phi)) {
    PHI <- .jackknife_validate_phi_matrix(phi, ranks, arg = "phi")
    n_ref <- nrow(PHI)
    q <- ncol(PHI)
    y_tk <- array(0.0, dim = c(Tn, s + 1L, q))
    for (k in seq_len(s + 1L)) {
      idx <- ranks[, k]
      y_tk[, k, ] <- PHI[idx, , drop = FALSE]
    }
  } else if (is.function(phi)) {
    n_ref <- max(ranks)
    # Probe q from the first call
    q_probe <- .jackknife_eval_phi(phi, ranks[1L, 1L], 1L)
    q <- length(q_probe)
    y_tk <- array(0.0, dim = c(Tn, s + 1L, q))
    y_tk[1L, 1L, ] <- q_probe
    for (t in seq_len(Tn)) {
      for (k in seq_len(s + 1L)) {
        if (t == 1L && k == 1L) next
        j <- ranks[t, k]
        y_tk[t, k, ] <- .jackknife_eval_phi(phi, j, t, q = q)
      }
    }
  } else {
    stop("`phi` must be either a numeric matrix (n_ref x q) or a function(j, t).", call. = FALSE)
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
