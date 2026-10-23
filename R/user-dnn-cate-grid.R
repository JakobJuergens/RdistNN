# CATE tuning-grid callbacks over frozen signals and fitted geometry.

.cate_dnn_rank_weight_spec <- function(method, k, s, normalize) {
  .dnn_rank_weight_spec(method, k, s, normalize)
}

.cate_dnn_grid_estimates <- function(tau, se) {
  out <- data.frame(
    target_id = seq_along(tau),
    component = NA_character_,
    estimate = as.numeric(tau),
    stringsAsFactors = FALSE
  )
  out$se <- if (is.null(se)) rep(NA_real_, length(tau)) else as.numeric(se)
  out
}

.cate_dnn_grid_default_precomputed <- function(
  candidates, tau_hat, se_hat, n, k, method, normalize
) {
  default_row <- candidates[candidates$is_default, , drop = FALSE]
  spec <- .cate_dnn_rank_weight_spec(
    method = method,
    k = k,
    s = default_row$s[[1L]],
    normalize = normalize
  )
  tab <- .rank_weight_table(spec, n)
  result <- list(
    estimates = .cate_dnn_grid_estimates(tau = tau_hat, se = se_hat),
    diagnostics = .dnn_grid_diagnostics_from_table(tab, n = n, ctx = list(k = k))
  )
  stats::setNames(list(result), as.character(default_row$computed_key[[1L]]))
}

.cate_dnn_grid_se_fn <- function(
  ranks, psi, x_mat, x0_mat, distance, method, k, normalize, ctx, jackknife, groups = NULL
) {
  force(ranks)
  force(psi)
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
      spec <- .cate_dnn_rank_weight_spec(method = method, k = k, s = s, normalize = normalize)
      loo <- if (is.null(tab_jk)) {
        .dnn_loo_loop(psi, ranks, spec, center = TRUE, groups = groups)
      } else {
        .dnn_loo_centered(
          psi,
          ranks,
          tab_jk,
          tab_full = attr(tab_jk, "tab_full", exact = TRUE),
          groups = groups
        )$loo
      }
      as.numeric(.dnn_jackknife_se_from_loo(
        loo = loo,
        q = 1L,
        n = nrow(ranks),
        m = ncol(ranks)
      )[1L, ])
    } else {
      .cate_dnn_jackknife_se_loo(
        x_mat = x_mat,
        x0_mat = x0_mat,
        psi = psi,
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

    data.frame(
      target_id = seq_along(se),
      component = NA_character_,
      se = as.numeric(se),
      stringsAsFactors = FALSE
    )
  }
}

.cate_dnn_warn_infeasible_grid <- function(candidates) {
  .dnn_warn_infeasible_grid(candidates)
}

.cate_dnn_grid_meta <- function(tuning_payload, engine_meta, n, p, k, method, distance, se, jackknife) {
  .dnn_estimator_grid_meta(tuning_payload, engine_meta, n, p, k, method, distance, se, jackknife)
}

.cate_dnn_tuning_grid <- function(
  candidates,
  ranks,
  tau_hat,
  se_hat,
  x_mat,
  x0_mat,
  psi,
  distance,
  method,
  k,
  normalize,
  ctx,
  se,
  jackknife,
  tuning_payload,
  reuse_default = TRUE, groups = NULL
) {
  .cate_dnn_warn_infeasible_grid(candidates)
  se_fn <- if (isTRUE(se) && !identical(jackknife, "none")) {
    .cate_dnn_grid_se_fn(
      ranks = ranks,
      psi = psi,
      x_mat = x_mat,
      x0_mat = x0_mat,
      distance = distance,
      method = method,
      k = k,
      normalize = normalize,
      ctx = ctx,
      jackknife = jackknife, groups = groups
    )
  } else {
    NULL
  }

  grid <- .dnn_grid_eval(
    candidates = candidates,
    ranks = ranks,
    weight_spec_fn = function(candidate_row) {
      .cate_dnn_rank_weight_spec(
        method = method,
        k = k,
        s = candidate_row$s[[1L]],
        normalize = normalize
      )
    },
    estimate_fn = function(W, candidate_row) {
      tau <- as.numeric(.dnn_grid_mean_estimate(
        psi, W, n = nrow(x_mat), arg = "pseudo_outcome", candidate_row = candidate_row)[1L, ])
      data.frame(
        target_id = seq_along(tau),
        component = NA_character_,
        estimate = tau,
        stringsAsFactors = FALSE
      )
    },
    se_fn = se_fn,
    ctx = list(k = k, groups = groups),
    precomputed = if (reuse_default && !.dnn_has_ties(groups)) .cate_dnn_grid_default_precomputed(
      candidates = candidates,
      tau_hat = tau_hat,
      se_hat = se_hat,
      n = nrow(x_mat),
      k = k,
      method = method,
      normalize = normalize
    ) else NULL
  )
  grid$meta <- .cate_dnn_grid_meta(
    tuning_payload = tuning_payload,
    engine_meta = grid$meta,
    n = nrow(x_mat),
    p = ncol(x_mat),
    k = k,
    method = method,
    distance = distance,
    se = se,
    jackknife = jackknife
  )
  grid
}
