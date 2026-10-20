# Diagnostics methods for user-level DNN estimator results.

.dnn_diag_top_share <- function(w, k) {
  k_use <- max(1L, min(length(w), as.integer(k)))
  sum(sort(as.numeric(w), decreasing = TRUE)[seq_len(k_use)])
}

.dnn_weighted_quantile <- function(x, w, prob = 0.5) {
  if (!is.numeric(x) || !is.numeric(w) || length(x) != length(w)) {
    return(NA_real_)
  }
  if (!is.numeric(prob) || length(prob) != 1L || is.na(prob) ||
      !is.finite(prob) || prob < 0 || prob > 1) {
    return(NA_real_)
  }
  if (.has_na_only(x) || .has_nonfinite_non_na(x)) {
    return(NA_real_)
  }
  if (.has_na_only(w) || .has_nonfinite_non_na(w)) {
    return(NA_real_)
  }
  if (any(w < 0)) {
    return(NA_real_)
  }

  ord <- order(x, seq_along(x))
  x_ord <- as.numeric(x[ord])
  w_ord <- as.numeric(w[ord])
  idx <- .weighted_cdf_index(w_ord, as.numeric(prob))
  if (is.na(idx)) return(NA_real_)
  x_ord[idx]
}

# Original-coordinate physical radius, independent of fitted localization.
# Callers provide canonical training/target matrices and nonnegative weights.
.dnn_physical_radius <- function(W, X, X0, prob = .9) {
  vapply(seq_len(ncol(W)), function(j) {
    .dnn_weighted_quantile(dist_l2(X, X0[j, ]), W[, j], prob)
  }, numeric(1L))
}

.dnn_tuning_diag_fields <- function(tuning_payload, m) {
  if (is.null(tuning_payload) || !is.list(tuning_payload)) {
    return(list(
      tuning_source = rep(NA_character_, m),
      tuning_name = rep(NA_character_, m),
      tuning_default_s = rep(NA_integer_, m),
      tuning_grid_available = rep(NA, m),
      tuning_grid_n_candidates = rep(NA_integer_, m),
      tuning_lambda_default = rep(NA_real_, m)
    ))
  }

  source <- if (!is.null(tuning_payload$source)) {
    as.character(tuning_payload$source)
  } else {
    NA_character_
  }
  name <- if (!is.null(tuning_payload$default$meta$name)) {
    as.character(tuning_payload$default$meta$name)
  } else if (!is.null(source)) {
    source
  } else {
    NA_character_
  }
  default_s <- if (!is.null(tuning_payload$default$s)) {
    as.integer(tuning_payload$default$s)
  } else {
    NA_integer_
  }
  has_grid <- !is.null(tuning_payload$grid$candidates) &&
    is.data.frame(tuning_payload$grid$candidates)
  grid_n <- if (isTRUE(has_grid)) nrow(tuning_payload$grid$candidates) else NA_integer_
  lambda <- if (!is.null(tuning_payload$default$meta$lambda)) {
    as.numeric(tuning_payload$default$meta$lambda)
  } else {
    NA_real_
  }

  list(
    tuning_source = rep(source, m),
    tuning_name = rep(name, m),
    tuning_default_s = rep(default_s, m),
    tuning_grid_available = rep(has_grid, m),
    tuning_grid_n_candidates = rep(as.integer(grid_n), m),
    tuning_lambda_default = rep(lambda, m)
  )
}

.dnn_weight_diagnostics <- function(W, meta = NULL, dx = NULL, tuning_payload = NULL) {
  w_mat <- as.matrix(W)
  if (!is.numeric(w_mat)) {
    .abort_bad_arg("W", "must be a numeric matrix", quoted = FALSE)
  }
  if (nrow(w_mat) < 1L || ncol(w_mat) < 1L) {
    .abort_bad_arg("W", "must have positive nrow and ncol", quoted = FALSE)
  }
  .abort_if_na_or_nonfinite(w_mat, arg = "W", quoted = FALSE)

  n <- nrow(w_mat)
  m <- ncol(w_mat)

  meta_n <- if (!is.null(meta$n)) as.integer(meta$n) else as.integer(n)
  meta_p <- if (!is.null(meta$p)) as.integer(meta$p) else NA_integer_
  meta_k <- if (!is.null(meta$k)) as.integer(meta$k) else NA_integer_
  meta_s <- if (!is.null(meta$s)) as.integer(meta$s) else NA_integer_
  meta_method <- if (!is.null(meta$method)) as.character(meta$method) else NA_character_
  meta_distance <- if (!is.null(meta$distance)) as.character(meta$distance) else NA_character_

  n_eff <- 1 / colSums(w_mat^2)
  max_weight <- apply(w_mat, 2L, max)
  top1 <- apply(w_mat, 2L, .dnn_diag_top_share, k = 1L)
  top5 <- apply(w_mat, 2L, .dnn_diag_top_share, k = 5L)
  top10 <- apply(w_mat, 2L, .dnn_diag_top_share, k = 10L)
  n_nonzero <- colSums(abs(w_mat) > sqrt(.Machine$double.eps))

  radius_rank_k <- rep(NA_real_, m)
  radius_rank_s <- rep(NA_real_, m)
  weighted_mean_distance <- rep(NA_real_, m)
  weighted_median_distance <- rep(NA_real_, m)
  distance_available <- rep(FALSE, m)

  if (!is.null(dx)) {
    d_mat <- as.matrix(dx)
    if (!is.numeric(d_mat)) {
      .abort_bad_arg("dx", "must be a numeric matrix", quoted = FALSE)
    }
    if (nrow(d_mat) != nrow(w_mat) || ncol(d_mat) != ncol(w_mat)) {
      .abort_bad_arg(
        "dx",
        "must have the same dimensions as `W`",
        expected = sprintf("nrow=%d, ncol=%d", nrow(w_mat), ncol(w_mat)),
        actual = sprintf("nrow=%d, ncol=%d", nrow(d_mat), ncol(d_mat)),
        quoted = FALSE
      )
    }
    .abort_if_na_or_nonfinite(d_mat, arg = "dx", quoted = FALSE)

    if (!is.na(meta_k) && meta_k >= 1L) {
      k_use <- max(1L, min(n, meta_k))
      radius_rank_k <- apply(d_mat, 2L, function(d) sort(d, partial = k_use)[k_use])
    }
    if (!is.na(meta_s) && meta_s >= 1L) {
      s_use <- max(1L, min(n, meta_s))
      radius_rank_s <- apply(d_mat, 2L, function(d) sort(d, partial = s_use)[s_use])
    }

    denom <- colSums(w_mat)
    good <- is.finite(denom) & denom > 0
    weighted_mean_distance[good] <- colSums(d_mat * w_mat)[good] / denom[good]

    weighted_median_distance <- vapply(
      seq_len(m),
      function(j) .dnn_weighted_quantile(d_mat[, j], w_mat[, j], prob = 0.5),
      numeric(1L)
    )
    distance_available <- rep(TRUE, m)
  }

  tuning <- .dnn_tuning_diag_fields(tuning_payload = tuning_payload, m = m)
  out <- data.frame(
    target_id = paste0("target_", seq_len(m)),
    n = rep(meta_n, m),
    p = rep(meta_p, m),
    k = rep(meta_k, m),
    s = rep(meta_s, m),
    method = rep(meta_method, m),
    distance = rep(meta_distance, m),
    n_eff = as.numeric(n_eff),
    max_weight = as.numeric(max_weight),
    top1_weight_share = as.numeric(top1),
    top5_weight_share = as.numeric(top5),
    top10_weight_share = as.numeric(top10),
    n_nonzero_weight = as.integer(n_nonzero),
    distance_diagnostics_available = as.logical(distance_available),
    radius_rank_k = as.numeric(radius_rank_k),
    radius_rank_s = as.numeric(radius_rank_s),
    weighted_mean_distance = as.numeric(weighted_mean_distance),
    weighted_median_distance = as.numeric(weighted_median_distance),
    tuning_source = tuning$tuning_source,
    tuning_name = tuning$tuning_name,
    tuning_default_s = tuning$tuning_default_s,
    tuning_grid_available = tuning$tuning_grid_available,
    tuning_grid_n_candidates = tuning$tuning_grid_n_candidates,
    tuning_lambda_default = tuning$tuning_lambda_default,
    stringsAsFactors = FALSE
  )
  rownames(out) <- out$target_id
  if (!is.null(tuning_payload$grid$candidates) && is.data.frame(tuning_payload$grid$candidates)) {
    attr(out, "tuning_grid") <- tuning_payload$grid$candidates
  }
  class(out) <- c("dnn_diagnostics", "data.frame")
  out
}

.dnn_diag_require_W <- function(object) {
  if (is.null(object$W)) {
    .abort_bad_arg(
      "object$W",
      "is required for diagnostics; refit with `return` including \"W\"",
      quoted = FALSE
    )
  }
  object$W
}

.dnn_diag_tuning_payload <- function(object) {
  if (!is.null(object$meta$tuning)) {
    object$meta$tuning
  } else {
    object$tuning
  }
}

.dnn_diag_append <- function(out, extra) {
  for (nm in names(extra)) {
    out[[nm]] <- extra[[nm]]
  }
  class(out) <- c("dnn_diagnostics", "data.frame")
  out
}

.dnn_candidate_diag_frame <- function(grid, selected) {
  diag <- grid$diagnostics
  required <- c(
    "candidate_id", "computed_key", "status", "reason", "s", "lambda",
    "ess", "max_weight", "top5_weight_share", "degenerate"
  )
  missing <- setdiff(required, names(diag))
  if (length(missing) > 0L) {
    .abort_bad_state(sprintf(
      "`object$tuning_grid$diagnostics` is missing required column(s): %s.",
      paste(missing, collapse = ", ")
    ))
  }
  diag$candidate_id <- as.character(diag$candidate_id)
  if (anyNA(diag$candidate_id) || any(!nzchar(diag$candidate_id))) {
    .abort_bad_state("`object$tuning_grid$diagnostics$candidate_id` must contain non-empty strings.")
  }
  duplicated_id <- unique(diag$candidate_id[duplicated(diag$candidate_id)])
  if (length(duplicated_id) > 0L) {
    .abort_bad_state(sprintf(
      "`object$tuning_grid$diagnostics` contains duplicate candidate row(s): %s.",
      paste(duplicated_id, collapse = ", ")
    ))
  }

  selected_id <- as.character(selected$candidate_id)
  missing_row <- setdiff(selected_id, diag$candidate_id)
  if (length(missing_row) > 0L) {
    .abort_bad_state(sprintf(
      "`object$tuning_grid$diagnostics` is missing row(s) for selected candidate(s): %s.",
      paste(missing_row, collapse = ", ")
    ))
  }

  diag <- diag[diag$candidate_id %in% selected_id, , drop = FALSE]
  diag$.candidate_order <- match(diag$candidate_id, selected_id)
  diag[order(diag$.candidate_order), , drop = FALSE]
}

.dnn_candidate_diagnostics <- function(object, candidates = "all") {
  grid <- .plot_require_grid(object, fn = "diagnostics")
  selected <- .plot_resolve_candidates(object, candidates = candidates, fn = "diagnostics")
  diag <- .dnn_candidate_diag_frame(grid = grid, selected = selected)
  diag$.candidate_order <- NULL
  idx <- match(diag$candidate_id, selected$candidate_id)
  top5 <- as.numeric(diag$top5_weight_share)

  meta <- object$meta
  out <- data.frame(
    candidate_id = diag$candidate_id,
    candidate_label = selected$candidate_label[idx],
    computed_key = diag$computed_key,
    is_default = as.logical(selected$is_default[idx]),
    status = diag$status,
    reason = diag$reason,
    n = rep(if (!is.null(meta$n)) as.integer(meta$n) else NA_integer_, nrow(diag)),
    p = rep(if (!is.null(meta$p)) as.integer(meta$p) else NA_integer_, nrow(diag)),
    k = rep(if (!is.null(meta$k)) as.integer(meta$k) else NA_integer_, nrow(diag)),
    s = as.integer(diag$s),
    lambda = suppressWarnings(as.numeric(diag$lambda)),
    method = rep(if (!is.null(meta$method)) as.character(meta$method) else NA_character_, nrow(diag)),
    distance = rep(if (!is.null(meta$distance)) as.character(meta$distance) else NA_character_, nrow(diag)),
    n_eff = as.numeric(diag$ess),
    max_weight = as.numeric(diag$max_weight),
    top1_weight_share = rep(NA_real_, nrow(diag)),
    top5_weight_share = top5,
    top10_weight_share = rep(NA_real_, nrow(diag)),
    n_nonzero_weight = rep(NA_integer_, nrow(diag)),
    distance_diagnostics_available = rep(FALSE, nrow(diag)),
    radius_rank_k = rep(NA_real_, nrow(diag)),
    radius_rank_s = rep(NA_real_, nrow(diag)),
    weighted_mean_distance = rep(NA_real_, nrow(diag)),
    weighted_median_distance = rep(NA_real_, nrow(diag)),
    degenerate = as.logical(diag$degenerate),
    stringsAsFactors = FALSE
  )
  rownames(out) <- out$candidate_id
  class(out) <- c("dnn_diagnostics", "data.frame")
  out
}

#' Compute diagnostics for fitted localized estimators
#'
#' Returns a compact diagnostics table for supported fit objects. Per-target
#' DNN diagnostics report fixed top-1, top-5, and top-10 weight concentration
#' shares. Distance diagnostics such as `radius_rank_k` and `radius_rank_s`
#' are ordinary nearest-rank scale summaries; `radius_rank_s` is not the DNN
#' support radius. For exact DNN fits the finite-population support edge is
#' rank `n - s + k`.
#'
#' PLQR fits additionally report residualized-design Gram conditioning
#' (`gram_min_eig`, `gram_condition_number`) and, when direct-sandwich
#' coefficient inference is available, target-level Jacobian conditioning:
#' `jacobian_available`, `jacobian_min_singular_value`,
#' `jacobian_condition_number`, `jacobian_rank`, density summaries over the
#' weight support (`density_min_support`, `density_max_support`, plus the
#' weighted `density_mean_weighted` and `density_median_weighted`), and
#' `theta_vcov_available`. A target can have an
#' acceptable residualized Gram but an ill-conditioned density-weighted
#' Jacobian, so coefficient-inference reliability should be judged from the
#' Jacobian block.
#'
#' @section Experimental methods:
#' Methods for `plqr_result` objects are experimental, as is all PLQR
#' functionality in this release; see [plqr_dnn()].
#'
#' @param object A fitted estimator object, currently an object returned by
#'   [np_dnn()], [cate_dnn()], [plqr_dnn()], [affine_ratio_dnn()] or [clate_dnn()].
#' @param candidates Candidate rows to report for grid fits. The default
#'   `"default"` preserves the per-target diagnostics table for the selected
#'   fit. Use `"all"` or a character vector of candidate ids to report
#'   candidate-level grid diagnostics from `object$tuning_grid`. Candidate
#'   grid diagnostics include `top5_weight_share`; the cutoff is recorded as
#'   `top_n = 5` and is independent of the estimator's `k`.
#' @param ... Not used.
#' @return An object of class `dnn_diagnostics`.
#' @name diagnostics
#' @examples
#' set.seed(1)
#' x_mat <- matrix(rnorm(40 * 2), 40, 2)
#' y_vec <- x_mat[, 1] + rnorm(40, sd = 0.2)
#' fit_np <- np_dnn(X = x_mat, x0 = c(0, 0), Y = y_vec, se = FALSE)
#' diagnostics(fit_np)
#'
#' demo_cate <- cate_demo_data(n = 50, p = 2, seed = 2)
#' fit_cate <- cate_dnn(
#'   data = demo_cate$data,
#'   x_cols = c("X1", "X2"),
#'   y_col = "Y",
#'   t_col = "D",
#'   m0_col = "m0_oracle", m1_col = "m1_oracle",
#'   e_col = "e_oracle",
#'   x0 = demo_cate$x0,
#'   se = FALSE
#' )
#' diagnostics(fit_cate)
#'
#' demo_plqr <- plqr_demo_data(n = 60, p = 2, p_d = 1, seed = 3)
#' fit_plqr <- plqr_dnn(
#'   data = demo_plqr$data,
#'   x_cols = c("X1", "X2"),
#'   D_cols = "D1",
#'   g_col = "g_oracle",
#'   mD_cols = "mD1",
#'   y_col = "Y",
#'   theta_eval = demo_plqr$truth$theta_true,
#'   x0 = demo_plqr$x0,
#'   se = FALSE
#' )
#' diagnostics(fit_plqr)
#' @export
diagnostics <- function(object, candidates = "default", ...) {
  UseMethod("diagnostics")
}

#' @rdname diagnostics
#' @export
diagnostics.default <- function(object, candidates = "default", ...) {
  .abort_bad_arg(
    "object",
    sprintf("diagnostics() is not implemented for class '%s'", class(object)[1L]),
    quoted = FALSE
  )
}

#' @rdname diagnostics
#' @param x A `dnn_diagnostics` object.
#' @param ... Not used.
#' @return `x`, invisibly.
#' @export
print.dnn_diagnostics <- function(x, ...) {
  cat("<dnn_diagnostics>\n", sep = "")
  if ("candidate_id" %in% names(x)) {
    ok <- as.character(x$status) == "ok"
    cat("  candidates: ", nrow(x), " (available: ", sum(ok), ")\n", sep = "")
  } else {
    cat("  targets: ", nrow(x), "\n", sep = "")
  }
  finite_neff <- x$n_eff[is.finite(x$n_eff)]
  if (length(finite_neff) > 0L) {
    cat(
      "  n_eff  : min=", round(min(finite_neff), 4),
      ", median=", round(stats::median(finite_neff), 4),
      ", max=", round(max(finite_neff), 4), "\n",
      sep = ""
    )
  }
  finite_top5 <- x$top5_weight_share[is.finite(x$top5_weight_share)]
  finite_max <- x$max_weight[is.finite(x$max_weight)]
  if (length(finite_top5) > 0L && length(finite_max) > 0L) {
    cat(
      "  weights: max(top5 share)=", round(max(finite_top5), 4),
      ", max(weight)=", round(max(finite_max), 4), "\n",
      sep = ""
    )
  }
  cat(
    "  dist   : available for ", sum(x$distance_diagnostics_available), "/", nrow(x),
    " target(s)\n",
    sep = ""
  )
  if ("outcome_available" %in% names(x)) {
    cat(
      "  outcome: available for ", sum(x$outcome_available), "/", nrow(x),
      " target(s)\n",
      sep = ""
    )
  }
  if ("overlap_available" %in% names(x)) {
    cat(
      "  overlap: available for ", sum(x$overlap_available), "/", nrow(x),
      " target(s)\n",
      sep = ""
    )
  }
  fields <- c(treated_share = "treated share", treated_ess = "treated ESS",
    control_ess = "control ESS", weighted_propensity_q10 = "propensity q10",
    weighted_propensity_q90 = "propensity q90",
    top5_score_share = "top-five absolute score-contribution share",
    instrument_1_mass = "instrument 1 mass", instrument_0_mass = "instrument 0 mass",
    propensity_clipped_mass = "clipped propensity mass",
    denominator_stat = "first-stage statistic",
    jk_top5_sq_share = "top-five squared jackknife-contribution share")
  for (field in names(fields)) .dnn_print_diagnostic_range(x, field, fields[[field]])
  if ("gram_available" %in% names(x)) {
    cat(
      "  gram   : available for ", sum(x$gram_available), "/", nrow(x),
      " target(s)\n",
      sep = ""
    )
  }
  invisible(x)
}
