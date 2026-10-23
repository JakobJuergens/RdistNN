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
#'   the package was compiled with C++ support, this calls the **C++ fast
#'   path**.
#' - Otherwise, it uses a general **R fallback** that also supports target-dependent
#'   per-neighbor contributions via a function `phi(j, t)`.
#' With matrix `phi`, `n_ref = nrow(phi)`. With function-valued `phi`, the
#' reference sample size is inferred as `max(ranks)`.
#'
#' @param ranks Integer matrix (T x (s+1)) of positive neighbor indices
#'   (1-based). Each row t lists the first s+1 neighbors of target t.
#'   Use distinct indices within each row for observation deletion. For repeated
#'   indices, the repeated-index convention adds each position's replacement change.
#' @param a Finite numeric length-s vector of L-stat coefficients a_1..a_s.
#' @param phi either:
#'   - finite numeric matrix (n_ref x q): row j is the q-vector phi(j), or
#'   - function `function(j, t)` returning a finite numeric length-q vector
#'     (allows target dependence).
#'   For a scalar \eqn{\phi}, use `cbind(y)` as the matrix case.
#' @param normalize Scalar logical; if TRUE (default) divide by sum(a).
#'
#' @return list with:
#'   - `theta0`: numeric length-q baseline aggregate, averaged across targets,
#'   - `theta_loo`: numeric matrix (n_ref x q) of delete-one aggregates,
#'   - `info`: list(s, T, q, normalize).
#' @examples
#' set.seed(1)
#' n <- 20
#' X <- matrix(rnorm(n * 2), n, 2)
#' wt <- weights_dnn(X, X0 = rbind(c(0, 0), c(1, 0)), k = 1L, s = 3L, return = "ord")
#' # A compact-support rule using the first three neighbors, plus a replacement.
#' ranks <- t(wt$ord[1:4, , drop = FALSE])
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
    cpp <- .rdistnn_compiled_fun("rustat_jackknife_fast_rank")
    if (!is.null(cpp)) {
      return(cpp(ranks, a, PHI, normalize))
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
  if (normalize && .compensated_sum(a) == 0) {
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

  .jackknife_finite_accumulation_R(ranks, a, y_tk, n_ref, normalize)
}
