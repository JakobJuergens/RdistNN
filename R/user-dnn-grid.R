# Shared candidate evaluation, cache reuse, failure rows and diagnostics.

.dnn_grid_candidate_s <- function(candidate_row) {
  if (!("s" %in% names(candidate_row))) {
    return(NA_integer_)
  }
  as.integer(candidate_row$s[[1L]])
}

.dnn_grid_candidate_lambda <- function(candidate_row) {
  if (!("lambda" %in% names(candidate_row))) {
    return(NA_real_)
  }
  suppressWarnings(as.numeric(candidate_row$lambda[[1L]]))
}

.dnn_grid_diag_top_n <- function(default = 5L) {
  .validate_scalar_int(default, "default")
}

.dnn_grid_check_candidates <- function(candidates) {
  if (!is.data.frame(candidates)) {
    .abort_bad_arg("candidates", "must be a data.frame")
  }
  required <- c(
    "candidate_id", "computed_key", "is_default", "feasible",
    "reason", "status"
  )
  missing <- setdiff(required, names(candidates))
  if (length(missing) > 0L) {
    .abort_bad_state(sprintf(
      "`candidates` is missing required column(s): %s.",
      paste(missing, collapse = ", ")
    ))
  }
  if (nrow(candidates) < 1L) {
    .abort_bad_state("`candidates` must contain at least one row.")
  }

  candidates <- as.data.frame(candidates, stringsAsFactors = FALSE)
  candidates$candidate_id <- as.character(candidates$candidate_id)
  candidates$computed_key <- as.character(candidates$computed_key)
  candidates$is_default <- as.logical(candidates$is_default)
  candidates$feasible <- as.logical(candidates$feasible)
  candidates$reason <- as.character(candidates$reason)
  candidates$status <- as.character(candidates$status)

  if (anyNA(candidates$candidate_id) || any(!nzchar(candidates$candidate_id))) {
    .abort_bad_state("`candidates$candidate_id` must contain non-empty strings.")
  }
  if (anyNA(candidates$computed_key) || any(!nzchar(candidates$computed_key))) {
    .abort_bad_state("`candidates$computed_key` must contain non-empty strings.")
  }
  if (anyNA(candidates$is_default) || sum(candidates$is_default) != 1L) {
    .abort_bad_state("`candidates$is_default` must contain exactly one TRUE value.")
  }
  if (anyNA(candidates$feasible)) {
    .abort_bad_state("`candidates$feasible` cannot contain missing values.")
  }
  if (anyNA(candidates$status) || any(!nzchar(candidates$status))) {
    .abort_bad_state("`candidates$status` must contain non-empty strings.")
  }
  candidates
}

.dnn_grid_check_estimates <- function(estimates) {
  if (!is.data.frame(estimates)) {
    .dnn_abort("callback_invalid", "`estimate_fn()` must return a data.frame.")
  }
  required <- c("target_id", "component", "estimate")
  missing <- setdiff(required, names(estimates))
  if (length(missing) > 0L) {
    .dnn_abort("callback_invalid", sprintf(
      "`estimate_fn()` output is missing required column(s): %s.",
      paste(missing, collapse = ", ")
    ))
  }
  if (nrow(estimates) < 1L) {
    .dnn_abort("callback_invalid", "`estimate_fn()` output must contain at least one row.")
  }
  if (!is.numeric(estimates$estimate)) {
    .dnn_abort("callback_invalid", "`estimate_fn()` column `estimate` must be numeric.")
  }
  estimates <- as.data.frame(estimates, stringsAsFactors = FALSE)
  estimates$estimate <- as.numeric(estimates$estimate)
  estimates
}

.dnn_grid_check_se <- function(se) {
  if (!is.data.frame(se)) {
    .dnn_abort("callback_invalid", "`se_fn()` must return a data.frame.")
  }
  required <- c("target_id", "component", "se")
  missing <- setdiff(required, names(se))
  if (length(missing) > 0L) {
    .dnn_abort("callback_invalid", sprintf(
      "`se_fn()` output is missing required column(s): %s.",
      paste(missing, collapse = ", ")
    ))
  }
  if (!is.numeric(se$se)) {
    .dnn_abort("callback_invalid", "`se_fn()` column `se` must be numeric.")
  }
  se <- as.data.frame(se, stringsAsFactors = FALSE)
  se$se <- as.numeric(se$se)
  se
}

.dnn_grid_target_key <- function(x) {
  component <- as.character(x$component)
  component[is.na(x$component)] <- "<NA>"
  paste(as.character(x$target_id), component, sep = "\r")
}

.dnn_grid_merge_estimate_se <- function(estimates, se) {
  if (is.null(se)) {
    estimates$se <- NA_real_
    return(estimates)
  }
  idx <- match(.dnn_grid_target_key(estimates), .dnn_grid_target_key(se))
  estimates$se <- se$se[idx]
  estimates
}

.dnn_grid_diagnostics_from_table <- function(tab, n, ctx) {
  top_n <- .dnn_grid_diag_top_n()
  if (is.null(tab)) {
    return(data.frame(
      ess = NA_real_,
      max_weight = NA_real_,
      top_n = as.integer(top_n),
      top5_weight_share = NA_real_,
      degenerate = NA,
      stringsAsFactors = FALSE
    ))
  }

  abs_tab <- abs(tab)
  denom <- sum(abs_tab)
  ess <- 1 / sum(tab^2)
  top5_share <- NA_real_
  if (!is.na(top_n) && denom > 0) {
    top_n_use <- min(length(tab), as.integer(top_n))
    top5_share <- sum(abs_tab[seq_len(top_n_use)]) / denom
  }
  data.frame(
    ess = ess,
    max_weight = max(abs_tab),
    top_n = as.integer(top_n),
    top5_weight_share = top5_share,
    degenerate = ess >= 0.9 * n,
    stringsAsFactors = FALSE
  )
}

.dnn_grid_compute_key <- function(
  candidate_row, ranks, weight_spec_fn, estimate_fn, se_fn, ctx, group_facts = NULL
) {
  n <- nrow(ranks)
  spec <- weight_spec_fn(candidate_row)
  if (!is_rank_weight(spec)) {
    .dnn_abort("callback_invalid", "`weight_spec_fn()` must return a rank_weight.")
  }

  tab <- .rank_weight_table(spec, n)
  W <- if (is.null(tab)) {
    rank_weight_eval_batch(spec, ranks, ctx = ctx)
  } else {
    matrix(tab[ranks], nrow = n, ncol = ncol(ranks))
  }

  groups <- if (is.list(ctx)) ctx$groups else NULL
  W <- .dnn_average_distance_weights(W, groups, facts = group_facts)
  estimates <- .dnn_grid_check_estimates(estimate_fn(W, candidate_row))
  se <- NULL
  if (!is.null(se_fn)) {
    tab_jk <- .rank_weight_table(spec, n - 1L)
    if (!is.null(tab_jk) && !is.null(tab)) {
      attr(tab_jk, "tab_full") <- tab
    }
    se <- .dnn_grid_check_se(se_fn(candidate_row, tab_jk))
  }

  list(
    estimates = .dnn_grid_merge_estimate_se(estimates, se),
    diagnostics = if (!.dnn_has_ties(groups, facts = group_facts)) {
      .dnn_grid_diagnostics_from_table(tab, n = n, ctx = ctx)
    } else {
      .dnn_grid_diagnostics_from_weights(W)
    }
  )
}

.dnn_grid_estimate_template <- function(results, m) {
  for (result in results) {
    if (!is.null(result$estimates)) {
      return(result$estimates[, c("target_id", "component"), drop = FALSE])
    }
  }
  data.frame(
    target_id = seq_len(m),
    component = NA_character_,
    stringsAsFactors = FALSE
  )
}

.dnn_grid_na_estimates <- function(template) {
  out <- template
  out$estimate <- NA_real_
  out$se <- NA_real_
  out
}

.dnn_grid_attach_estimate_candidate <- function(estimates, candidate_row) {
  data.frame(
    candidate_id = as.character(candidate_row$candidate_id[[1L]]),
    computed_key = as.character(candidate_row$computed_key[[1L]]),
    lambda = .dnn_grid_candidate_lambda(candidate_row),
    s = .dnn_grid_candidate_s(candidate_row),
    target_id = estimates$target_id,
    component = estimates$component,
    estimate = estimates$estimate,
    se = estimates$se,
    status = as.character(candidate_row$status[[1L]]),
    reason = as.character(candidate_row$reason[[1L]]),
    stringsAsFactors = FALSE
  )
}

.dnn_grid_empty_diagnostics <- function(ctx = NULL) {
  top_n <- .dnn_grid_diag_top_n()
  data.frame(
    ess = NA_real_,
    max_weight = NA_real_,
    top_n = as.integer(top_n),
    top5_weight_share = NA_real_,
    degenerate = NA,
    stringsAsFactors = FALSE
  )
}

.dnn_grid_diag_field <- function(diagnostic, name, default) {
  if (name %in% names(diagnostic)) {
    return(diagnostic[[name]][[1L]])
  }
  default
}

.dnn_grid_attach_diagnostic_candidate <- function(diagnostic, candidate_row) {
  top5 <- .dnn_grid_diag_field(diagnostic, "top5_weight_share", NA_real_)
  data.frame(
    candidate_id = as.character(candidate_row$candidate_id[[1L]]),
    computed_key = as.character(candidate_row$computed_key[[1L]]),
    lambda = .dnn_grid_candidate_lambda(candidate_row),
    s = .dnn_grid_candidate_s(candidate_row),
    ess = diagnostic$ess[[1L]],
    max_weight = diagnostic$max_weight[[1L]],
    top_n = .dnn_grid_diag_field(diagnostic, "top_n", .dnn_grid_diag_top_n()),
    top5_weight_share = top5,
    degenerate = diagnostic$degenerate[[1L]],
    status = as.character(candidate_row$status[[1L]]),
    reason = as.character(candidate_row$reason[[1L]]),
    stringsAsFactors = FALSE
  )
}

.dnn_grid_meta <- function(candidates, ranks, unique_keys) {
  list(
    n = as.integer(nrow(ranks)),
    m = as.integer(ncol(ranks)),
    n_candidates = as.integer(nrow(candidates)),
    n_feasible = as.integer(sum(candidates$feasible)),
    n_skipped = as.integer(sum(!candidates$feasible)),
    n_failed = as.integer(sum(candidates$status == "failed")),
    n_unique = as.integer(length(unique_keys))
  )
}

# NP/CATE own this declaration: a validated weighted-mean range failure is
# an expected optional-candidate outcome. The main order keeps its original
# numerical error, and generic callback/kernel failures remain fatal.
.dnn_grid_mean_estimate <- function(signal, W, n, arg, candidate_row) {
  tryCatch(.dnn_estimate_from_signal(signal, W, n = n, arg = arg),
    rdistnn_numerical_error = function(err) {
      if (isTRUE(candidate_row$is_default[[1L]])) stop(err)
      .dnn_abort("candidate_failure", conditionMessage(err),
        argument = err$argument, parent = err)
    })
}

.dnn_grid_eval <- function(
  candidates,
  ranks,
  weight_spec_fn,
  estimate_fn,
  se_fn = NULL,
  ctx = NULL,
  precomputed = NULL
) {
  candidates <- .dnn_grid_check_candidates(candidates)
  ranks <- .rank_check_mat(ranks)
  group_facts <- .dnn_group_facts(if (is.list(ctx)) ctx$groups else NULL)
  if (!is.function(weight_spec_fn)) {
    .abort_bad_arg("weight_spec_fn", "must be a function", quoted = FALSE)
  }
  if (!is.function(estimate_fn)) {
    .abort_bad_arg("estimate_fn", "must be a function", quoted = FALSE)
  }
  if (!is.null(se_fn) && !is.function(se_fn)) {
    .abort_bad_arg("se_fn", "must be NULL or a function", quoted = FALSE)
  }

  feasible_keys <- unique(candidates$computed_key[candidates$feasible])
  if (length(feasible_keys) < 1L) {
    .abort_bad_arg("candidates", "must contain at least one feasible candidate")
  }
  if (nrow(candidates) > 1L && length(unique(candidates$computed_key)) == 1L) {
    .dnn_warn("grid_single_value",
      "All feasible tuning candidates share the same computed key; the grid is effectively a single value."
    )
  }

  results <- setNames(vector("list", length(feasible_keys)), feasible_keys)
  if (!is.null(precomputed)) {
    if (!is.list(precomputed) || is.null(names(precomputed))) {
      .abort_bad_arg("precomputed", "must be a named list", quoted = FALSE)
    }
    for (key in intersect(names(precomputed), feasible_keys)) {
      results[[key]] <- precomputed[[key]]
    }
  }
  for (key in feasible_keys) {
    if (!is.null(results[[key]])) {
      next
    }
    idx <- which(candidates$feasible & candidates$computed_key == key)
    candidate_row <- candidates[idx[[1L]], , drop = FALSE]
    key_is_default <- any(candidates$is_default[idx])
    result <- tryCatch(
      .dnn_grid_compute_key(
        candidate_row = candidate_row,
        ranks = ranks,
        weight_spec_fn = weight_spec_fn,
        estimate_fn = estimate_fn,
        se_fn = se_fn,
        ctx = ctx,
        group_facts = group_facts
      ),
      rdistnn_candidate_failure = function(err) {
        if (key_is_default) {
          stop(err)
        }
        list(error = err, reason = conditionMessage(err))
      }
    )
    if (!is.null(result$error)) {
      failed_idx <- candidates$feasible & candidates$computed_key == key
      candidates$status[failed_idx] <- "failed"
      candidates$reason[failed_idx] <- result$reason
    }
    results[[key]] <- result
  }

  failed <- candidates$status == "failed"
  if (any(failed)) {
    .dnn_warn("partial_result",
      sprintf(paste0("Evaluation failed for %d optional tuning candidate(s). ",
        "Inspect `tuning_grid$candidates` for all candidate statuses and reasons."), sum(failed)),
      affected = sum(failed), examples = candidates$candidate_id[failed])
  }

  template <- .dnn_grid_estimate_template(results, m = ncol(ranks))
  estimates <- vector("list", nrow(candidates))
  diagnostics <- vector("list", nrow(candidates))
  for (i in seq_len(nrow(candidates))) {
    candidate_row <- candidates[i, , drop = FALSE]
    key <- as.character(candidate_row$computed_key[[1L]])
    result <- results[[key]]
    ok <- isTRUE(candidate_row$feasible[[1L]]) &&
      identical(as.character(candidate_row$status[[1L]]), "ok") &&
      !is.null(result) &&
      is.null(result$error)

    base_est <- if (ok) result$estimates else .dnn_grid_na_estimates(template)
    base_diag <- if (ok) result$diagnostics else .dnn_grid_empty_diagnostics(ctx = ctx)
    estimates[[i]] <- .dnn_grid_attach_estimate_candidate(base_est, candidate_row)
    diagnostics[[i]] <- .dnn_grid_attach_diagnostic_candidate(base_diag, candidate_row)
  }

  list(
    candidates = candidates,
    estimates = do.call(rbind, estimates),
    diagnostics = do.call(rbind, diagnostics),
    meta = .dnn_grid_meta(candidates, ranks, feasible_keys)
  )
}

# Shared descriptive grid policy; estimator adapters own inference and statuses.
.dnn_warn_infeasible_grid <- function(candidates) {
  skipped <- candidates[!candidates$feasible, c("candidate_id", "reason"), drop = FALSE]
  if (nrow(skipped) < 1L) {
    return(invisible(NULL))
  }
  .dnn_warn("infeasible_candidates",
    sprintf(paste0("Skipped infeasible tuning candidate(s): %d. ",
      "Inspect `tuning_grid$candidates` for all candidate statuses and reasons."), nrow(skipped)),
    affected = nrow(skipped), examples = skipped$candidate_id)
  invisible(NULL)
}

.dnn_estimator_grid_meta <- function(tuning_payload, engine_meta, n, p, k, method,
                                     distance, se, jackknife, p_d = NULL) {
  grid_meta <- tuning_payload$grid$meta
  name <- if (!is.null(grid_meta$name)) {
    as.character(grid_meta$name)
  } else if (!is.null(tuning_payload$default$meta$name)) {
    as.character(tuning_payload$default$meta$name)
  } else {
    as.character(tuning_payload$source)
  }
  out <- list(
    name = name,
    n = as.integer(n),
    kx = as.integer(p),
    k = as.integer(k),
    method = method,
    distance = .dnn_distance_name(distance),
    jackknife = if (isTRUE(se)) jackknife else "none",
    n_candidates = engine_meta$n_candidates,
    n_feasible = engine_meta$n_feasible,
    n_skipped = engine_meta$n_skipped,
    n_failed = engine_meta$n_failed,
    n_unique = engine_meta$n_unique
  )
  if (!is.null(p_d)) {
    out <- append(out, list(p_d = as.integer(p_d)), after = match("kx", names(out)))
  }
  out
}

# Candidate-level summaries use actual target weights after block averaging.
.dnn_grid_diagnostics_from_weights <- function(W) {
  rows <- lapply(seq_len(ncol(W)), function(j) {
    .dnn_grid_diagnostics_from_table(sort(abs(W[, j]), decreasing = TRUE), nrow(W), NULL)
  })
  out <- rows[[1L]]
  out$ess <- min(vapply(rows, function(row) row$ess, numeric(1L)))
  for (nm in c("max_weight", "top5_weight_share")) {
    out[[nm]] <- max(vapply(rows, function(row) row[[nm]], numeric(1L)))
  }
  out$degenerate <- any(vapply(rows, function(row) row$degenerate, logical(1L)))
  out
}
