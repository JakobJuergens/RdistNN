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

#' @name diagnostics
#' @export
diagnostics.plqr_result <- function(object, ...) {
  out <- .dnn_weight_diagnostics(
    W = .dnn_diag_require_W(object),
    meta = object$meta,
    dx = object$dx,
    tuning_payload = .dnn_diag_tuning_payload(object)
  )
  model <- object$.model
  d_resid_mat <- if (is.list(model)) model$d_resid_mat else NULL
  .dnn_diag_append(out, .plqr_gram_diagnostics(object$W, d_resid_mat = d_resid_mat))
}
