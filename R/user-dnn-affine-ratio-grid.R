.affine_grid_estimates <- function(fit) {
  data.frame(target_id = seq_along(fit$theta), component = NA_character_,
             estimate = fit$theta,
             se = if (is.null(fit$se)) rep(NA_real_, length(fit$theta)) else fit$se)
}

.affine_ratio_tuning_grid <- function(default, signal, ranks, model, active, jackknife, tol, tuning,
                                     reuse_default = TRUE) {
  if (!is.null(default$.model$x_mat)) {
    model$groups <- default$.model$groups
  }
  candidates <- tuning$grid$candidates
  .dnn_warn_infeasible_grid(candidates)
  cache <- new.env(parent = emptyenv())
  default_row <- candidates[candidates$is_default, , drop = FALSE]
  key <- default_row$computed_key[[1L]]
  cache[[key]] <- default
  spec <- function(row) {
    candidate_model <- model
    candidate_model$s <- row$s[[1L]]
    .affine_rank_weight(candidate_model)
  }
  precomputed <- if (reuse_default && !.dnn_has_ties(model$groups)) stats::setNames(list(list(
    estimates = .affine_grid_estimates(default),
    diagnostics = .dnn_grid_diagnostics_from_table(.rank_weight_table(spec(default_row), nrow(signal)),
                                                  nrow(signal), list(k = model$k))
  )), key) else NULL
  grid <- .dnn_grid_eval(candidates, ranks, weight_spec_fn = spec,
    estimate_fn = function(W, row) {
      fit <- .affine_ratio_evaluate(signal, W, ranks, spec(row), active, jackknife, tol, groups = model$groups)
      cache[[row$computed_key[[1L]]]] <- fit
      .affine_grid_estimates(fit)
    },
    se_fn = if (active) function(row, tab_jk) .affine_grid_estimates(cache[[row$computed_key[[1L]]]]) else NULL,
    ctx = list(k = model$k, groups = model$groups), precomputed = precomputed)
  # The common grid status describes candidate evaluation. Preserve the
  # separate per-target numerical status needed for ratio inference.
  grid$estimates$ratio_status <- grid$estimates$covariance_status <- NA_character_
  grid$estimates$wald_available <- FALSE
  for (i in seq_len(nrow(grid$estimates))) {
    row <- grid$estimates[i, ]
    if (row$status != "ok") next
    fit <- cache[[row$computed_key]]
    j <- row$target_id
    grid$estimates$ratio_status[i] <- fit$status$ratio_status[j]
    grid$estimates$covariance_status[i] <- fit$status$covariance_status[j]
    grid$estimates$wald_available[i] <- fit$status$wald_available[j]
  }
  grid$meta <- c(list(name = tuning$default$meta$name, kx = default$meta$p,
                     k = model$k, method = model$method, normalize = model$normalize,
                     distance = default$meta$distance, jackknife = default$meta$jackknife), grid$meta)
  grid
}
