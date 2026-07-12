# CATE-specific diagnostics for cate_result objects.

.cate_paper_diagnostics <- function(
  W, x_mat, x0_mat, signal, tau_hat, t_vec = NULL, radius_mass = 0.90
) {
  w_mat <- as.matrix(W)
  x_mat <- as.matrix(x_mat)
  x0_mat <- as.matrix(x0_mat)
  signal <- as.numeric(signal)
  tau_hat <- as.numeric(tau_hat)
  n <- nrow(w_mat)
  m <- ncol(w_mat)

  if (nrow(x_mat) != n || ncol(x_mat) != ncol(x0_mat)) {
    .abort_bad_state("CATE paper diagnostics require conformable X, X0, and W.")
  }
  if (nrow(x0_mat) == 1L && m > 1L) {
    x0_mat <- x0_mat[rep(1L, m), , drop = FALSE]
  }
  if (nrow(x0_mat) != m) {
    .abort_bad_state("CATE paper diagnostics require one target per W column.")
  }
  score_available <- length(signal) == n && length(tau_hat) == m
  col_sums <- colSums(w_mat)
  if (any(!is.finite(w_mat)) || any(!is.finite(col_sums)) || any(col_sums <= 0)) {
    .abort_bad_state("CATE paper diagnostic weights must be finite with positive column sums.")
  }
  w_mat <- sweep(w_mat, 2L, col_sums, FUN = "/")

  physical_radius_90 <- vapply(seq_len(m), function(j) {
    distance <- sqrt(rowSums(sweep(x_mat, 2L, x0_mat[j, ], FUN = "-")^2))
    ord <- order(distance)
    idx <- which(cumsum(w_mat[ord, j]) >= radius_mass)[1L]
    distance[ord[[idx]]]
  }, numeric(1L))

  arm_ess <- function(arm) {
    if (is.null(t_vec)) return(rep(NA_real_, m))
    idx <- as.numeric(t_vec) == arm
    vapply(seq_len(m), function(j) {
      w <- w_mat[idx, j]
      if (length(w) == 0L || sum(w^2) <= 0) NA_real_ else sum(w)^2 / sum(w^2)
    }, numeric(1L))
  }

  top5_score_share <- if (score_available) {
    vapply(seq_len(m), function(j) {
      contribution <- abs(w_mat[, j] * (signal - tau_hat[[j]]))
      total <- sum(contribution)
      if (!is.finite(total) || total <= 0) return(NA_real_)
      sum(sort(contribution, decreasing = TRUE)[seq_len(min(5L, n))]) / total
    }, numeric(1L))
  } else {
    rep(NA_real_, m)
  }

  data.frame(
    physical_radius_90 = physical_radius_90,
    arm_ess_available = rep(!is.null(t_vec), m),
    treated_ess = arm_ess(1L),
    control_ess = arm_ess(0L),
    score_contribution_available = rep(score_available, m),
    top5_score_share = top5_score_share,
    stringsAsFactors = FALSE
  )
}

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
diagnostics.cate_result <- function(object, candidates = "default", ...) {
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
  t_vec <- if (is.list(model)) model$t_vec else NULL
  e <- if (is.list(model)) model$e else NULL
  out <- .dnn_diag_append(out, .cate_overlap_diagnostics(object$W, t_vec = t_vec, e = e))
  .dnn_diag_append(
    out,
    .cate_paper_diagnostics(
      W = object$W,
      x_mat = model$x_mat,
      x0_mat = object$x0,
      signal = model$psi,
      tau_hat = object$tau,
      t_vec = t_vec
    )
  )
}
