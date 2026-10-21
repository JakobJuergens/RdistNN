# PLQR moment-grid callbacks; coefficient sandwich policy stays in the core.

.plqr_dnn_rank_weight_spec <- function(method, k, s, normalize) {
  if (identical(method, "exact")) {
    rank_weight_dnn_exact(k = as.integer(k), s = as.integer(s))
  } else {
    rank_weight_dnn_asymp(k = as.integer(k), s = as.integer(s), normalize = normalize)
  }
}

.plqr_dnn_grid_estimates <- function(moment, se) {
  component <- rownames(moment)
  if (is.null(component)) {
    component <- paste0("comp_", seq_len(nrow(moment)))
  }
  out <- data.frame(
    target_id = rep(seq_len(ncol(moment)), each = nrow(moment)),
    component = rep(component, times = ncol(moment)),
    estimate = as.numeric(moment),
    stringsAsFactors = FALSE
  )
  out$se <- if (is.null(se)) rep(NA_real_, nrow(out)) else as.numeric(se)
  out
}

.plqr_dnn_grid_default_precomputed <- function(
  candidates, moment_hat, se_hat, n, k, method, normalize
) {
  default_row <- candidates[candidates$is_default, , drop = FALSE]
  spec <- .plqr_dnn_rank_weight_spec(
    method = method,
    k = k,
    s = default_row$s[[1L]],
    normalize = normalize
  )
  tab <- .rank_weight_table(spec, n)
  result <- list(
    estimates = .plqr_dnn_grid_estimates(moment = moment_hat, se = se_hat),
    diagnostics = .dnn_grid_diagnostics_from_table(tab, n = n, ctx = list(k = k))
  )
  stats::setNames(list(result), as.character(default_row$computed_key[[1L]]))
}

.plqr_dnn_grid_se_fn <- function(
  ranks, score_mat, x_mat, x0_mat, distance, method, k, normalize, ctx, jackknife
) {
  force(ranks)
  force(score_mat)
  force(x_mat)
  force(x0_mat)
  force(distance)
  force(method)
  force(k)
  force(normalize)
  force(ctx)
  force(jackknife)

  function(candidate_row, tab_jk) {
    s <- as.integer(candidate_row$s[[1L]])
    se <- if (identical(jackknife, "fast")) {
      spec <- .plqr_dnn_rank_weight_spec(method = method, k = k, s = s, normalize = normalize)
      loo <- if (is.null(tab_jk)) {
        .dnn_loo_loop(score_mat, ranks, spec, center = TRUE)
      } else {
        .dnn_loo_centered(
          score_mat,
          ranks,
          tab_jk,
          tab_full = attr(tab_jk, "tab_full", exact = TRUE)
        )$loo
      }
      .dnn_jackknife_se_from_loo(
        loo = loo,
        q = ncol(score_mat),
        n = nrow(ranks),
        m = ncol(ranks)
      )
    } else {
      .plqr_dnn_jackknife_se_loo(
        x_mat = x_mat,
        x0_mat = x0_mat,
        score_mat = score_mat,
        distance = distance,
        method = method,
        k = k,
        s = s,
        normalize = normalize,
        ctx = ctx,
        jackknife = "loo",
        precomputed_ranks = NULL
      )
    }

    out <- .plqr_dnn_grid_estimates(moment = se, se = NULL)
    data.frame(
      target_id = out$target_id,
      component = out$component,
      se = out$estimate,
      stringsAsFactors = FALSE
    )
  }
}

.plqr_dnn_warn_infeasible_grid <- function(candidates) {
  .dnn_warn_infeasible_grid(candidates)
}

.plqr_dnn_grid_meta <- function(
  tuning_payload, engine_meta, n, p, p_d, k, method, distance, se, jackknife
) {
  .dnn_estimator_grid_meta(tuning_payload, engine_meta, n, p, k, method, distance, se, jackknife, p_d = p_d)
}

.plqr_dnn_tuning_grid <- function(
  candidates,
  ranks,
  moment_hat,
  se_hat,
  x_mat,
  x0_mat,
  score_mat,
  distance,
  method,
  k,
  normalize,
  ctx,
  se,
  jackknife,
  tuning_payload,
  reuse_default = TRUE
) {
  .plqr_dnn_warn_infeasible_grid(candidates)
  se_fn <- if (isTRUE(se) && !identical(jackknife, "none")) {
    .plqr_dnn_grid_se_fn(
      ranks = ranks,
      score_mat = score_mat,
      x_mat = x_mat,
      x0_mat = x0_mat,
      distance = distance,
      method = method,
      k = k,
      normalize = normalize,
      ctx = ctx,
      jackknife = jackknife
    )
  } else {
    NULL
  }

  grid <- .dnn_grid_eval(
    candidates = candidates,
    ranks = ranks,
    weight_spec_fn = function(candidate_row) {
      .plqr_dnn_rank_weight_spec(
        method = method,
        k = k,
        s = candidate_row$s[[1L]],
        normalize = normalize
      )
    },
    estimate_fn = function(W, candidate_row) {
      moment <- .dnn_estimate_from_signal(score_mat, W, n = nrow(x_mat), arg = "score")
      rownames(moment) <- paste0("comp_", seq_len(nrow(moment)))
      colnames(moment) <- paste0("target_", seq_len(ncol(moment)))
      .plqr_dnn_grid_estimates(moment = moment, se = NULL)
    },
    se_fn = se_fn,
    ctx = list(k = k),
    precomputed = if (reuse_default) .plqr_dnn_grid_default_precomputed(
      candidates = candidates,
      moment_hat = moment_hat,
      se_hat = se_hat,
      n = nrow(x_mat),
      k = k,
      method = method,
      normalize = normalize
    ) else NULL
  )
  grid$meta <- .plqr_dnn_grid_meta(
    tuning_payload = tuning_payload,
    engine_meta = grid$meta,
    n = nrow(x_mat),
    p = ncol(x_mat),
    p_d = ncol(score_mat),
    k = k,
    method = method,
    distance = distance,
    se = se,
    jackknife = jackknife
  )
  grid
}
