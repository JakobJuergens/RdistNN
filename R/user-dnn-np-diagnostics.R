# Nonparametric-regression-specific diagnostics for np_result objects.

.np_outcome_diagnostics <- function(W, y_vec = NULL) {
  w_mat <- as.matrix(W)
  n <- nrow(w_mat)
  m <- ncol(w_mat)

  unavailable <- function() {
    data.frame(
      outcome_available = rep(FALSE, m),
      weighted_mean_y = rep(NA_real_, m),
      weighted_sd_y = rep(NA_real_, m),
      weighted_y_q10 = rep(NA_real_, m),
      weighted_y_q50 = rep(NA_real_, m),
      weighted_y_q90 = rep(NA_real_, m),
      stringsAsFactors = FALSE
    )
  }

  if (is.null(y_vec)) {
    return(unavailable())
  }

  y_num <- as.numeric(y_vec)
  if (length(y_num) != n) {
    .abort_bad_arg(
      "object$.model$y_vec",
      "must have length equal to nrow(object$W)",
      expected = sprintf("nrow(object$W)=%d", n),
      actual = sprintf("length=%d", length(y_num)),
      quoted = FALSE
    )
  }
  .abort_if_na_or_nonfinite(y_num, arg = "object$.model$y_vec", quoted = FALSE)

  w_sum <- colSums(w_mat)
  good <- is.finite(w_sum) & w_sum > 0
  weighted_mean_y <- rep(NA_real_, m)
  weighted_sd_y <- rep(NA_real_, m)
  weighted_mean_y[good] <- as.numeric(crossprod(y_num, w_mat))[good] / w_sum[good]

  for (j in which(good)) {
    centered <- y_num - weighted_mean_y[j]
    weighted_sd_y[j] <- sqrt(sum(w_mat[, j] * centered^2) / w_sum[j])
  }

  weighted_y_q10 <- vapply(
    seq_len(m),
    function(j) .dnn_weighted_quantile(y_num, w_mat[, j], prob = 0.1),
    numeric(1L)
  )
  weighted_y_q50 <- vapply(
    seq_len(m),
    function(j) .dnn_weighted_quantile(y_num, w_mat[, j], prob = 0.5),
    numeric(1L)
  )
  weighted_y_q90 <- vapply(
    seq_len(m),
    function(j) .dnn_weighted_quantile(y_num, w_mat[, j], prob = 0.9),
    numeric(1L)
  )

  data.frame(
    outcome_available = rep(TRUE, m),
    weighted_mean_y = weighted_mean_y,
    weighted_sd_y = weighted_sd_y,
    weighted_y_q10 = weighted_y_q10,
    weighted_y_q50 = weighted_y_q50,
    weighted_y_q90 = weighted_y_q90,
    stringsAsFactors = FALSE
  )
}

#' @name diagnostics
#' @export
diagnostics.np_result <- function(object, ...) {
  out <- .dnn_weight_diagnostics(
    W = .dnn_diag_require_W(object),
    meta = object$meta,
    dx = object$dx,
    tuning_payload = .dnn_diag_tuning_payload(object)
  )
  model <- object$.model
  y_vec <- if (is.list(model)) model$y_vec else object$y
  .dnn_diag_append(out, .np_outcome_diagnostics(object$W, y_vec = y_vec))
}
