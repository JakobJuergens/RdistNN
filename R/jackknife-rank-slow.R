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
#'   gives the positive 1-based index of the k-th neighbor, for `k = 1, ..., s+1`.
#'   (Must exclude `t` itself if targets are a subset of reference.)
#' @param a Finite numeric length-`s` vector of L-stat coefficients \eqn{(a_1,\dots,a_s)}.
#' @param phi Either a finite numeric matrix (n_ref x q) with row `j` equal to
#'   \eqn{\phi(j)}, or a function `phi(j, t)` returning a finite length-`q`
#'   numeric vector.
#' @param normalize Scalar logical; if `TRUE` (default) divide aggregates by `sum(a)`.
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
  args <- .jackknife_validate_finite_support_args(ranks, a, normalize)
  ranks <- args$ranks
  a <- args$a
  normalize <- args$normalize

  Tn <- nrow(ranks)
  s <- length(a)

  # Build per-(t,k) contribution array y_tk: (T x (s+1) x q)
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
    # probe q from first call
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

  Z <- if (normalize) .compensated_sum_R(a) else 1
  baseline <- matrix(0, Tn, q)
  for (t in seq_len(Tn)) for (component in seq_len(q)) {
    baseline[t, component] <- .compensated_sum_R(a * y_tk[t, seq_len(s), component]) / Z
  }
  theta0 <- vapply(seq_len(q), function(component) {
    .compensated_sum_R(baseline[, component]) / Tn
  }, numeric(1L))
  theta_loo <- matrix(0, n_ref, q)
  # Independent reference: rebuild each single-position neighbor deletion.
  # Repeated indices retain the fast API's additive single-position convention.
  for (i in seq_len(n_ref)) {
    terms <- lapply(seq_len(q), function(component) baseline[, component])
    for (t in seq_len(Tn)) {
      positions <- which(ranks[t, seq_len(s)] == i)
      for (position in positions) {
        keep <- seq_len(s + 1L) != position
        for (component in seq_len(q)) {
          replacement <- .compensated_sum_R(a * y_tk[t, keep, component]) / Z
          terms[[component]] <- c(terms[[component]], -baseline[t, component], replacement)
        }
      }
    }
    for (component in seq_len(q)) {
      theta_loo[i, component] <- .compensated_sum_R(terms[[component]]) / Tn
    }
  }
  list(theta0 = theta0, theta_loo = theta_loo,
    info = list(s = s, T = Tn, q = q, normalize = normalize))
}
