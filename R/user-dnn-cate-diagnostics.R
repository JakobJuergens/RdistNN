# CATE-specific diagnostics for cate_result objects.

.cate_overlap_diagnostics <- function(W, t_vec = NULL, e = NULL) {
  w_mat <- as.matrix(W)
  n <- nrow(w_mat)
  m <- ncol(w_mat)

  unavailable <- function(propensity_available = FALSE) {
    data.frame(
      overlap_available = rep(FALSE, m),
      treated_mass = rep(NA_real_, m),
      control_mass = rep(NA_real_, m),
      treated_share = rep(NA_real_, m),
      propensity_available = rep(propensity_available, m),
      weighted_mean_propensity = rep(NA_real_, m),
      weighted_propensity_q10 = rep(NA_real_, m),
      weighted_propensity_q90 = rep(NA_real_, m),
      stringsAsFactors = FALSE
    )
  }

  if (is.null(t_vec)) {
    return(unavailable(propensity_available = FALSE))
  }

  t_num <- as.numeric(t_vec)
  if (length(t_num) != n) {
    .abort_bad_arg(
      "object$.model$t_vec",
      "must have length equal to nrow(object$W)",
      expected = sprintf("nrow(object$W)=%d", n),
      actual = sprintf("length=%d", length(t_num)),
      quoted = FALSE
    )
  }
  .abort_if_na_or_nonfinite(t_num, arg = "object$.model$t_vec", quoted = FALSE)
  if (any(!(t_num %in% c(0, 1)))) {
    .abort_bad_arg("object$.model$t_vec", "must be binary in {0,1}", quoted = FALSE)
  }

  treated_mass <- as.numeric(crossprod(t_num, w_mat))
  control_mass <- as.numeric(crossprod(1 - t_num, w_mat))
  denom <- treated_mass + control_mass
  treated_share <- rep(NA_real_, m)
  good <- is.finite(denom) & denom > 0
  treated_share[good] <- treated_mass[good] / denom[good]

  propensity_available <- !is.null(e)
  weighted_mean_propensity <- rep(NA_real_, m)
  weighted_propensity_q10 <- rep(NA_real_, m)
  weighted_propensity_q90 <- rep(NA_real_, m)
  if (isTRUE(propensity_available)) {
    e_num <- as.numeric(e)
    if (length(e_num) != n) {
      .abort_bad_arg(
        "object$.model$e",
        "must have length equal to nrow(object$W)",
        expected = sprintf("nrow(object$W)=%d", n),
        actual = sprintf("length=%d", length(e_num)),
        quoted = FALSE
      )
    }
    .abort_if_na_or_nonfinite(e_num, arg = "object$.model$e", quoted = FALSE)
    w_sum <- colSums(w_mat)
    good_w <- is.finite(w_sum) & w_sum > 0
    weighted_mean_propensity[good_w] <- as.numeric(crossprod(e_num, w_mat))[good_w] / w_sum[good_w]
    weighted_propensity_q10 <- vapply(
      seq_len(m),
      function(j) .dnn_weighted_quantile(e_num, w_mat[, j], prob = 0.1),
      numeric(1L)
    )
    weighted_propensity_q90 <- vapply(
      seq_len(m),
      function(j) .dnn_weighted_quantile(e_num, w_mat[, j], prob = 0.9),
      numeric(1L)
    )
  }

  data.frame(
    overlap_available = rep(TRUE, m),
    treated_mass = treated_mass,
    control_mass = control_mass,
    treated_share = treated_share,
    propensity_available = rep(propensity_available, m),
    weighted_mean_propensity = weighted_mean_propensity,
    weighted_propensity_q10 = weighted_propensity_q10,
    weighted_propensity_q90 = weighted_propensity_q90,
    stringsAsFactors = FALSE
  )
}

#' @name diagnostics
#' @export
diagnostics.cate_result <- function(object, ...) {
  out <- .dnn_weight_diagnostics(
    W = .dnn_diag_require_W(object),
    meta = object$meta,
    dx = object$dx,
    tuning_payload = .dnn_diag_tuning_payload(object)
  )
  model <- object$.model
  t_vec <- if (is.list(model)) model$t_vec else NULL
  e <- if (is.list(model)) model$e else NULL
  .dnn_diag_append(out, .cate_overlap_diagnostics(object$W, t_vec = t_vec, e = e))
}
