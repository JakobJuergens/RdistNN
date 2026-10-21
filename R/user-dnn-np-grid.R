# NP tuning-grid callbacks over retained outcomes and fitted geometry.

.np_dnn_rank_weight_spec <- function(method, k, s, normalize) {
  .dnn_rank_weight_spec(method, k, s, normalize)
}

.np_dnn_grid_estimates <- function(mu, se) {
  out <- data.frame(
    target_id = seq_along(mu),
    component = NA_character_,
    estimate = as.numeric(mu),
    stringsAsFactors = FALSE
  )
  out$se <- if (is.null(se)) rep(NA_real_, length(mu)) else as.numeric(se)
  out
}

.np_dnn_grid_default_precomputed <- function(
  candidates, mu_hat, se_hat, n, k, method, normalize
) {
  default_row <- candidates[candidates$is_default, , drop = FALSE]
  spec <- .np_dnn_rank_weight_spec(
    method = method,
    k = k,
    s = default_row$s[[1L]],
    normalize = normalize
  )
  tab <- .rank_weight_table(spec, n)
  result <- list(
    estimates = .np_dnn_grid_estimates(mu = mu_hat, se = se_hat),
    diagnostics = .dnn_grid_diagnostics_from_table(tab, n = n, ctx = list(k = k))
  )
  stats::setNames(list(result), as.character(default_row$computed_key[[1L]]))
}

.np_dnn_grid_se_fn <- function(
  ranks, y_vec, x_mat, x0_mat, distance, method, k, normalize, ctx, jackknife
) {
  force(ranks)
  force(y_vec)
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
      spec <- .np_dnn_rank_weight_spec(method = method, k = k, s = s, normalize = normalize)
      loo <- if (is.null(tab_jk)) {
        .dnn_loo_loop(y_vec, ranks, spec, center = TRUE)
      } else {
        .dnn_loo_centered(
          y_vec,
          ranks,
          tab_jk,
          tab_full = attr(tab_jk, "tab_full", exact = TRUE)
        )$loo
      }
      as.numeric(.dnn_jackknife_se_from_loo(
        loo = loo,
        q = 1L,
        n = nrow(ranks),
        m = ncol(ranks)
      )[1L, ])
    } else {
      .np_dnn_jackknife_se_loo(
        x_mat = x_mat,
        x0_mat = x0_mat,
        y = y_vec,
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

.np_dnn_warn_infeasible_grid <- function(candidates) {
  .dnn_warn_infeasible_grid(candidates)
}

.np_dnn_grid_meta <- function(tuning_payload, engine_meta, n, p, k, method, distance, se, jackknife) {
  .dnn_estimator_grid_meta(tuning_payload, engine_meta, n, p, k, method, distance, se, jackknife)
}

.np_dnn_tuning_grid <- function(
  candidates,
  ranks,
  mu_hat,
  se_hat,
  x_mat,
  x0_mat,
  y_vec,
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
  .np_dnn_warn_infeasible_grid(candidates)
  se_fn <- if (isTRUE(se) && !identical(jackknife, "none")) {
    .np_dnn_grid_se_fn(
      ranks = ranks,
      y_vec = y_vec,
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
      .np_dnn_rank_weight_spec(
        method = method,
        k = k,
        s = candidate_row$s[[1L]],
        normalize = normalize
      )
    },
    estimate_fn = function(W, candidate_row) {
      mu <- as.numeric(.dnn_estimate_from_signal(y_vec, W, n = nrow(x_mat), arg = "Y")[1L, ])
      data.frame(
        target_id = seq_along(mu),
        component = NA_character_,
        estimate = mu,
        stringsAsFactors = FALSE
      )
    },
    se_fn = se_fn,
    ctx = list(k = k),
    precomputed = if (reuse_default) .np_dnn_grid_default_precomputed(
      candidates = candidates,
      mu_hat = mu_hat,
      se_hat = se_hat,
      n = nrow(x_mat),
      k = k,
      method = method,
      normalize = normalize
    ) else NULL
  )
  grid$meta <- .np_dnn_grid_meta(
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
