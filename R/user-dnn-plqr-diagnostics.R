# PLQR-specific diagnostics for plqr_result objects.

.plqr_gram_diagnostics <- function(W, d_resid_mat = NULL, rank_tol = sqrt(.Machine$double.eps)) {
  w_mat <- as.matrix(W)
  n <- nrow(w_mat)
  m <- ncol(w_mat)

  unavailable <- function() {
    data.frame(
      gram_available = rep(FALSE, m),
      gram_min_eig = rep(NA_real_, m),
      gram_max_eig = rep(NA_real_, m),
      gram_condition_number = rep(NA_real_, m),
      gram_rank_tol = rep(as.numeric(rank_tol), m),
      stringsAsFactors = FALSE
    )
  }

  if (is.null(d_resid_mat)) {
    return(unavailable())
  }

  d_resid <- as.matrix(d_resid_mat)
  if (!is.numeric(d_resid)) {
    .abort_bad_arg("object$.model$d_resid_mat", "must be numeric", quoted = FALSE)
  }
  if (nrow(d_resid) != n) {
    .abort_bad_arg(
      "object$.model$d_resid_mat",
      "must have row count equal to nrow(object$W)",
      expected = sprintf("nrow(object$W)=%d", n),
      actual = sprintf("nrow=%d", nrow(d_resid)),
      quoted = FALSE
    )
  }
  .abort_if_na_or_nonfinite(d_resid, arg = "object$.model$d_resid_mat", quoted = FALSE)

  gram_min_eig <- gram_max_eig <- gram_condition_number <- rep(NA_real_, m)
  for (j in seq_len(m)) {
    G <- crossprod(d_resid, d_resid * as.numeric(w_mat[, j]))
    G <- (G + t(G)) / 2
    eig <- eigen(G, symmetric = TRUE, only.values = TRUE)$values
    gram_min_eig[j] <- min(eig)
    gram_max_eig[j] <- max(eig)
    if (gram_min_eig[j] <= rank_tol) {
      gram_condition_number[j] <- Inf
    } else {
      gram_condition_number[j] <- gram_max_eig[j] / gram_min_eig[j]
    }
  }

  data.frame(
    gram_available = rep(TRUE, m),
    gram_min_eig = gram_min_eig,
    gram_max_eig = gram_max_eig,
    gram_condition_number = gram_condition_number,
    gram_rank_tol = rep(as.numeric(rank_tol), m),
    stringsAsFactors = FALSE
  )
}

.plqr_jacobian_diagnostics <- function(object, W, rank_tol = sqrt(.Machine$double.eps)) {
  w_mat <- as.matrix(W)
  m <- ncol(w_mat)

  jacobian_available <- rep(FALSE, m)
  jacobian_min_singular_value <- rep(NA_real_, m)
  jacobian_condition_number <- rep(NA_real_, m)
  jacobian_rank <- rep(NA_integer_, m)
  density_min_support <- rep(NA_real_, m)
  density_mean_weighted <- rep(NA_real_, m)
  density_median_weighted <- rep(NA_real_, m)
  density_max_support <- rep(NA_real_, m)
  theta_vcov_available <- rep(FALSE, m)

  jacobian <- object$jacobian
  if (is.list(jacobian) && length(jacobian) == m) {
    for (j in seq_len(m)) {
      J <- jacobian[[j]]
      if (!is.matrix(J) || !is.numeric(J) || any(!is.finite(J))) {
        next
      }
      jacobian_available[j] <- TRUE
      sv <- svd(J, nu = 0L, nv = 0L)$d
      jacobian_min_singular_value[j] <- min(sv)
      jacobian_rank[j] <- sum(sv > rank_tol * max(sv, 1))
      jacobian_condition_number[j] <- if (min(sv) <= rank_tol) {
        Inf
      } else {
        max(sv) / min(sv)
      }
    }
  }

  theta_vcov <- object$theta_vcov
  if (is.list(theta_vcov) && length(theta_vcov) == m) {
    theta_vcov_available <- vapply(
      theta_vcov,
      function(v) is.matrix(v) && is.numeric(v) && all(is.finite(v)),
      logical(1L)
    )
  }

  model <- object$.model
  density <- if (is.list(model)) model$density else NULL
  if (is.numeric(density) && length(density) == nrow(w_mat)) {
    support_tol <- sqrt(.Machine$double.eps)
    for (j in seq_len(m)) {
      w_j <- as.numeric(w_mat[, j])
      support <- w_j > support_tol
      if (any(support)) {
        # Support extrema: plain min/max of the density over observations
        # carrying weight, not weighted statistics.
        density_min_support[j] <- min(density[support])
        density_max_support[j] <- max(density[support])
        density_mean_weighted[j] <- sum(w_j * density)
        density_median_weighted[j] <- .dnn_weighted_quantile(
          density[support],
          w_j[support],
          prob = 0.5
        )
      }
    }
  }

  data.frame(
    jacobian_available = jacobian_available,
    jacobian_min_singular_value = jacobian_min_singular_value,
    jacobian_condition_number = jacobian_condition_number,
    jacobian_rank = jacobian_rank,
    density_min_support = density_min_support,
    density_mean_weighted = density_mean_weighted,
    density_median_weighted = density_median_weighted,
    density_max_support = density_max_support,
    theta_vcov_available = theta_vcov_available,
    stringsAsFactors = FALSE
  )
}

#' @name diagnostics
#' @export
diagnostics.plqr_result <- function(object, candidates = "default", ...) {
  if (!.plot_candidates_is_default(candidates)) {
    return(.dnn_candidate_diagnostics(object, candidates = candidates))
  }
  out <- .dnn_weight_diagnostics(
    W = .dnn_diag_require_W(object),
    meta = object$meta,
    dx = object$dx,
    tuning_payload = .dnn_diag_tuning_payload(object)
  )
  model <- object$.model
  d_resid_mat <- if (is.list(model)) model$d_resid_mat else NULL
  out <- .dnn_diag_append(out, .plqr_gram_diagnostics(object$W, d_resid_mat = d_resid_mat))
  .dnn_diag_append(out, .plqr_jacobian_diagnostics(object, W = object$W))
}
