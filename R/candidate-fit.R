# Per-candidate views over grid fits.

.candidate_grid_payload <- function(object) {
  grid <- object[["tuning_grid", exact = TRUE]]
  if (is.null(grid)) {
    .dnn_abort("extra_unavailable",
      paste0("This fit does not contain a `tuning_grid` payload; ",
        "refit with `tuning_mode = \"grid\"` before calling candidate_fit()."))
  }
  scope <- .dnn_reporting_scope()
  cached <- if (!is.null(scope)) scope$validated_grid else NULL
  if (!is.null(cached) && identical(grid, cached$grid) &&
      identical(object$meta$s, cached$fitted_s)) return(grid)
  if (!is.list(grid) || !all(c("candidates", "estimates", "diagnostics", "meta") %in% names(grid)) ||
      !is.data.frame(grid$candidates) ||
      !is.data.frame(grid$estimates) || !is.data.frame(grid$diagnostics) || !is.list(grid$meta)) {
    .dnn_abort("state_invalid", "`object$tuning_grid` is malformed.",
      hint = "Refit with the current estimator before reusing this fit.")
  }
  .dnn_grid_check_candidates(grid$candidates, retained = TRUE)
  jackknife <- grid$meta[["jackknife", exact = TRUE]]
  if (!is.character(jackknife) || !is.null(dim(jackknife)) || length(jackknife) != 1L ||
      is.na(jackknife) || !jackknife %in% c("none", "fast", "loo")) {
    .dnn_abort("state_invalid", "`object$tuning_grid$meta$jackknife` is missing or invalid.",
      hint = "Refit with the current estimator before reusing this fit.")
  }
  default_s <- grid$candidates$s[grid$candidates$is_default]
  if (!is.numeric(default_s) || !is.null(dim(default_s)) || length(default_s) != 1L ||
      !is.finite(default_s) || default_s < 1 || default_s != trunc(default_s) ||
      !is.numeric(object$meta$s) || length(object$meta$s) != 1L || default_s != object$meta$s) {
    .dnn_abort("state_invalid", paste0(
      "`object$tuning_grid$candidates$s` for the default candidate must be a positive integer ",
      "matching the fitted order. Refit with the current estimator."))
  }
  .candidate_validate_estimate_rows(object, grid$estimates)
  .dnn_candidate_diag_frame(grid)
  if (!is.null(scope)) scope$validated_grid <- list(grid = grid, fitted_s = object$meta$s)
  grid
}

.candidate_resolve_row <- function(object, candidate) {
  grid <- .candidate_grid_payload(object)
  candidates <- grid$candidates
  if (missing(candidate) || length(candidate) != 1L) {
    .abort_bad_arg(
      "candidate",
      "must be a single candidate id string or integer row index",
      quoted = FALSE
    )
  }

  if (is.character(candidate)) {
    idx <- match(candidate, candidates$candidate_id)
    if (is.na(idx)) {
      .abort_bad_arg(
        "candidate",
        sprintf(
          "unknown candidate id '%s'; use one of: %s",
          .dnn_brief(candidate),
          .dnn_message_examples(candidates$candidate_id)$text
        ),
        quoted = FALSE
      )
    }
  } else {
    idx <- .validate_scalar_int(candidate, "candidate")
    if (idx < 1L || idx > nrow(candidates)) {
      .abort_bad_arg(
        "candidate",
        sprintf("row index must be in 1..%d", nrow(candidates)),
        quoted = FALSE
      )
    }
  }

  row <- candidates[idx, , drop = FALSE]
  status <- as.character(row$status[[1L]])
  if (!identical(status, "ok")) {
    reason <- as.character(row$reason[[1L]])
    .abort_bad_arg(
      "candidate",
      sprintf(
        "candidate '%s' is %s: %s",
        as.character(row$candidate_id[[1L]]),
        status,
        reason
      ),
      quoted = FALSE
    )
  }
  s <- row$s
  if (!is.numeric(s) || is.complex(s) || length(s) != 1L || !is.finite(s) ||
      s < 1 || s > .Machine$integer.max || s != trunc(s)) {
    .dnn_abort("state_invalid",
      "`object$tuning_grid$candidates$s` must be a positive integer for the selected candidate.")
  }
  .candidate_validate_retained_controls(object, row)
  row
}

.candidate_meta_int <- function(object, field) {
  value <- object$meta[[field, exact = TRUE]]
  if (is.null(value)) {
    .dnn_abort("state_invalid", sprintf("`object$meta$%s` is missing. Refit the current estimator.", field))
  }
  if (!is.numeric(value) || is.complex(value) || length(value) != 1L || !is.finite(value) ||
      value < -.Machine$integer.max || value > .Machine$integer.max || value != trunc(value)) {
    .dnn_abort("state_invalid", sprintf("`object$meta$%s` must be an integer scalar.", field))
  }
  as.integer(value)
}

.candidate_expected_m <- function(object) {
  meta_m <- .candidate_meta_int(object, "m")
  x0_m <- NULL
  if (!is.null(object$x0)) {
    x0_m <- nrow(as.matrix(object$x0))
  }
  if (!is.null(meta_m) && !is.null(x0_m) && x0_m != meta_m) {
    .dnn_abort("state_invalid", sprintf(
      "`object$x0` row count must match `object$meta$m` (got nrow(object$x0)=%d, m=%d).",
      x0_m,
      meta_m
    ))
  }
  meta_m
}

# Candidate controls are required even when no optional weights are rebuilt.
# Check the two retained control copies and stored computation identities before
# attaching a new order to previously calculated estimates.
.candidate_validate_retained_controls <- function(object, row) {
  .dnn_validate_retained({
    n <- .candidate_meta_int(object, "n")
    .candidate_expected_m(object)
    model <- object$.model
    required <- c("x_mat", "method", "k", "s", "normalize")
    if (!is.list(model)) {
      .dnn_abort("state_invalid", "This object lacks retained candidate controls. Refit with the current estimator.")
    }
    missing <- required[vapply(required, function(nm) is.null(model[[nm]]), logical(1L))]
    if (length(missing)) {
      .dnn_abort("state_invalid",
        sprintf("`object$.model$%s` is missing. Refit with the current estimator.", missing[[1L]]))
    }
    X <- .validate_matrix(model$x_mat, "object$.model$x_mat",
      min_rows = 1L, min_cols = 1L, numeric = TRUE, finite = TRUE)
    if (nrow(X) != n || ncol(X) != .candidate_meta_int(object, "p")) {
      .dnn_abort("state_invalid",
        "Retained covariates are inconsistent with `object$meta$n/p`. Refit the current estimator.")
    }
    method <- object$meta$method
    if (!is.character(method) || length(method) != 1L || is.na(method) ||
        !method %in% c("exact", "asymp")) {
      .dnn_abort("state_invalid", "`object$meta$method` must be 'exact' or 'asymp'. Refit the current estimator.")
    }
    if (!identical(model$method, method)) {
      .dnn_abort("state_invalid",
        "`object$.model$method` is inconsistent with `object$meta$method`. Refit the current estimator.")
    }
    .validate_scalar_logical(model$normalize, "object$.model$normalize")
    k <- .candidate_meta_int(object, "k")
    parent_s <- .candidate_meta_int(object, "s")
    model_k <- .validate_scalar_int(model$k, "object$.model$k")
    model_s <- .validate_scalar_int(model$s, "object$.model$s")
    if (k != model_k || parent_s != model_s) {
      .dnn_abort("state_invalid",
        "Retained `k/s` controls are inconsistent between metadata and model state. Refit the current estimator."
      )
    }
    s <- as.integer(row$s[[1L]])
    if (k < 1L || k > s || s > n || parent_s < k || parent_s > n) {
      .dnn_abort("state_invalid",
        "Retained candidate controls must satisfy 1 <= k <= s <= n. Refit the current estimator.")
    }
    # A scalar tuning display can request a new uncertainty policy. Its grid
    # owns candidate feasibility; the original fit's policy stays unchanged.
    jackknife <- object$tuning_grid$meta$jackknife
    if (jackknife %in% c("fast", "loo") && s >= n) {
      .dnn_abort("state_invalid",
        "Retained candidate jackknife controls must satisfy s < n. Refit the current estimator.")
    }
    key <- row$computed_key
    if (!is.character(key) || length(key) != 1L || is.na(key) ||
        !identical(key, sprintf("s_%d", s))) {
      .dnn_abort("state_invalid",
        "The selected candidate's `computed_key` is inconsistent with its order. Refit the current estimator."
      )
    }
    estimates <- object$tuning_grid$estimates
    required_rows <- c("candidate_id", "computed_key", "s")
    if (!all(required_rows %in% names(estimates))) {
      .dnn_abort("state_invalid",
        "Stored candidate estimates lack computation identities. Refit the current estimator.")
    }
    selected <- estimates[!is.na(estimates$candidate_id) &
      estimates$candidate_id == row$candidate_id[[1L]], , drop = FALSE]
    if (!nrow(selected)) {
      .dnn_abort("state_invalid", sprintf("No estimates recorded for candidate '%s'.", row$candidate_id[[1L]]))
    }
    .candidate_validate_order_identity(selected, row, "estimates")
    invisible(NULL)
  })
}

.candidate_validate_order_identity <- function(rows, candidate_row, payload) {
  s <- candidate_row$s[[1L]]
  key <- candidate_row$computed_key[[1L]]
  if (!is.numeric(rows$s) || is.complex(rows$s) ||
      anyNA(rows$s) || any(rows$s != s) ||
      !is.character(rows$computed_key) || anyNA(rows$computed_key) ||
      any(rows$computed_key != key)) {
    .dnn_abort("state_invalid",
      sprintf("Stored candidate %s are inconsistent with the selected order. Refit the current estimator.", payload)
    )
  }
  invisible(NULL)
}

.candidate_validate_estimate_rows <- function(object, rows, id = NULL) {
  required <- c("candidate_id", "computed_key", "lambda", "s", "target_id", "component",
    "estimate", "se", "status", "reason")
  missing <- setdiff(required, names(rows))
  if (length(missing) > 0L) {
    .dnn_abort("state_invalid", sprintf(
      "`object$tuning_grid$estimates` is missing required column(s): %s.",
      paste(missing, collapse = ", ")
    ), hint = "Refit with the current estimator before reusing this fit.")
  }
  invalid <- function(field, reason) {
    .dnn_abort("state_invalid", sprintf(
      "Stored candidate estimates `%s` %s. Refit with the current estimator.", field, reason))
  }
  for (field in c("candidate_id", "computed_key", "status", "reason")) {
    value <- rows[[field]]
    if (!is.character(value) || !is.null(dim(value)) || anyNA(value)) {
      invalid(field, "must retain non-missing character values")
    }
  }
  for (field in c("lambda", "s", "estimate", "se", "target_id")) {
    value <- rows[[field]]
    if (!is.numeric(value) || is.complex(value) || !is.null(dim(value))) {
      invalid(field, "must retain real numeric values")
    }
  }
  if (!is.character(rows$component) || !is.null(dim(rows$component))) {
    invalid("component", "must retain character values")
  }
  candidates <- object$tuning_grid$candidates
  index <- match(rows$candidate_id, candidates$candidate_id)
  if (anyNA(index)) invalid("candidate_id", "contains unknown candidate identities")
  for (field in c("status", "reason")) {
    if (!identical(rows[[field]], candidates[[field]][index])) {
      invalid(field, "is inconsistent with the candidate table")
    }
  }
  # The full retained table is checked before selection. Scientific availability
  # below belongs to a selected successful candidate, including partial ratios.
  if (is.null(id)) return(invisible(TRUE))

  paired <- inherits(object, "rdistnn_affine_ratio_result")
  for (field in c("estimate", "se")) {
    value <- rows[[field]]
    optional <- paired || (field == "se" && identical(object$tuning_grid$meta$jackknife, "none"))
    if (!is.numeric(value) || !is.null(dim(value)) ||
        any(!is.finite(value) & !(optional & is.na(value) & !is.nan(value))) ||
        (field == "se" && (any(value < 0, na.rm = TRUE) ||
          (!paired && optional && any(!is.na(value)))))) {
      .dnn_abort("state_invalid", sprintf(
        paste0("Stored candidate `%s` values have invalid type, range or availability. ",
          "Refit with the current estimator."), field))
    }
  }
  component <- rows$component
  scalar <- inherits(object, c("rdistnn_np_result", "rdistnn_cate_result"))
  if (!is.character(component) || !is.null(dim(component)) ||
      ((scalar || paired) && any(!is.na(component)))) {
    .dnn_abort("state_invalid", paste0(
      "Stored candidate `component` values are inconsistent with the fitted estimator. ",
      "Refit with the current estimator."))
  }
  if (paired) {
    invalid <- function(field, reason) {
      .dnn_abort("state_invalid",
        sprintf("Stored candidate `%s` %s. Refit with the current estimator.", field, reason))
    }
    .affine_validate_retained_status(rows, rows$estimate, rows$se,
      !identical(object$tuning_grid$meta$jackknife, "none"), invalid, "tuning_grid$estimates")
  }

  expected_m <- .candidate_expected_m(object)
  target_id <- rows$target_id
  if (!is.numeric(target_id) || is.complex(target_id) || !is.null(dim(target_id)) ||
      length(target_id) != nrow(rows) || any(!is.finite(target_id)) ||
      any(target_id < -.Machine$integer.max | target_id > .Machine$integer.max) ||
      any(target_id != trunc(target_id))) {
    .dnn_abort("state_invalid", sprintf("Estimate rows for candidate '%s' must have integer `target_id` values.", id))
  }
  target_id <- as.integer(target_id)
  if (!is.null(expected_m)) {
    invalid <- target_id < 1L | target_id > expected_m
    if (any(invalid)) {
      .dnn_abort("state_invalid", sprintf(
        "Estimate rows for candidate '%s' have `target_id` outside 1..%d.",
        id,
        expected_m
      ))
    }
    missing_target <- setdiff(seq_len(expected_m), unique(target_id))
    if (length(missing_target) > 0L) {
      .dnn_abort("state_invalid", sprintf(
        "Estimate rows for candidate '%s' are missing target_id value(s): %s.",
        id,
        paste(missing_target, collapse = ", ")
      ))
    }
  }

  if (inherits(object, c("rdistnn_np_result", "rdistnn_cate_result", "rdistnn_affine_ratio_result"))) {
    duplicated_target <- unique(target_id[duplicated(target_id)])
    if (length(duplicated_target) > 0L) {
      .dnn_abort("state_invalid", sprintf(
        "Estimate rows for candidate '%s' contain duplicate target_id value(s): %s.",
        id,
        paste(duplicated_target, collapse = ", ")
      ))
    }
  }

  if (inherits(object, "rdistnn_plqr_result") && !is.null(object$moment)) {
    expected_components <- rownames(as.matrix(object$moment))
    if (!is.null(expected_components)) {
      component <- as.character(rows$component)
      if (anyNA(component) || any(!nzchar(component))) {
        .dnn_abort("state_invalid",
          sprintf("Estimate rows for candidate '%s' must have non-empty `component` values.", id))
      }
      unknown_component <- setdiff(unique(component), expected_components)
      if (length(unknown_component) > 0L) {
        .dnn_abort("state_invalid", sprintf(
          "Estimate rows for candidate '%s' contain unknown component value(s): %s.",
          id,
          paste(unknown_component, collapse = ", ")
        ))
      }
      pair <- paste(target_id, component, sep = "\r")
      duplicated_pair <- unique(pair[duplicated(pair)])
      if (length(duplicated_pair) > 0L) {
        .dnn_abort("state_invalid",
          sprintf("Estimate rows for candidate '%s' contain duplicate target/component rows.", id))
      }
      if (!is.null(expected_m)) {
        expected_pair <- as.vector(outer(seq_len(expected_m), expected_components, paste, sep = "\r"))
        if (length(setdiff(expected_pair, pair)) > 0L) {
          .dnn_abort("state_invalid",
            sprintf("Estimate rows for candidate '%s' are missing target/component rows.", id))
        }
      }
    }
  }

  invisible(TRUE)
}

.candidate_estimates <- function(object, candidate_row) {
  grid <- .candidate_grid_payload(object)
  id <- as.character(candidate_row$candidate_id[[1L]])
  rows <- grid$estimates[grid$estimates$candidate_id == id, , drop = FALSE]
  if (nrow(rows) < 1L) {
    .dnn_abort("state_invalid", sprintf("No estimates recorded for candidate '%s'.", id))
  }
  .candidate_validate_estimate_rows(object, rows, id)
  rows[order(rows$target_id, rows$component), , drop = FALSE]
}

.candidate_rank_weight_spec <- function(object, candidate_row) {
  method <- object$meta$method
  if (!is.character(method) || length(method) != 1L || is.na(method) ||
      !method %in% c("exact", "asymp")) {
    .dnn_abort("state_invalid", "`object$meta$method` must be 'exact' or 'asymp'. Refit the current estimator.")
  }
  normalize <- .validate_scalar_logical(object$.model$normalize, "object$.model$normalize")
  k <- .candidate_meta_int(object, "k")
  s <- as.integer(candidate_row$s[[1L]])
  .dnn_rank_weight_spec(method, k, s, normalize)
}

.candidate_validate_cached_ranks <- function(object, ranks) {
  expected_n <- .candidate_meta_int(object, "n")
  if (!is.null(expected_n) && nrow(ranks) != expected_n) {
    .dnn_abort("state_invalid", sprintf(
      "Cached ranks row count must match `object$meta$n` (got nrow=%d, n=%d).",
      nrow(ranks),
      expected_n
    ))
  }
  expected_m <- .candidate_expected_m(object)
  if (!is.null(expected_m) && ncol(ranks) != expected_m) {
    .dnn_abort("state_invalid", sprintf(
      "Cached ranks column count must match the number of target points (got ncol=%d, m=%d).",
      ncol(ranks),
      expected_m
    ))
  }
  ranks
}

.candidate_cached_ranks <- function(object) {
  ranks <- .dnn_retained_ranking(object)$ranks
  .candidate_validate_cached_ranks(object, ranks)
}

.candidate_weights_from_ranks <- function(object, candidate_row) {
  ranks <- .candidate_cached_ranks(object)
  spec <- .candidate_rank_weight_spec(object, candidate_row)
  tab <- .rank_weight_table(spec, nrow(ranks))
  ctx <- object$.model$ctx
  W <- if (is.null(tab)) {
    rank_weight_eval_batch(spec, ranks, ctx = ctx)
  } else {
    matrix(tab[ranks], nrow = nrow(ranks), ncol = ncol(ranks))
  }
  .dnn_average_distance_weights(W, object$.model$groups)
}

.candidate_tuning_payload <- function(object, candidate_row) {
  # The fitted-state owner validates required tuning once in this operation.
  # Fresh scalar sensitivity grids do not transform a fitted candidate view.
  .dnn_check_report_state(object)
  payload <- object$meta[["tuning", exact = TRUE]]
  payload$default$s <- as.integer(candidate_row$s[[1L]])
  lambda <- candidate_row[["lambda", exact = TRUE]]
  if (!is.numeric(lambda) || !is.null(dim(lambda)) || length(lambda) != 1L) {
    .dnn_abort("state_invalid", paste0(
      "The selected candidate must retain its numeric `lambda` label. ",
      "Refit with the current estimator."))
  }
  payload$default$meta$lambda <- as.numeric(lambda)
  payload$default$meta$candidate_id <- as.character(candidate_row$candidate_id[[1L]])
  payload$default$feasibility <- list(
    k_le_s = as.integer(object$meta$k) <= as.integer(candidate_row$s[[1L]]),
    s_lt_n = identical(object$meta$jackknife, "none") ||
      as.integer(candidate_row$s[[1L]]) < as.integer(object$meta$n)
  )
  payload
}

.candidate_update_model <- function(object, candidate_row) {
  model <- object$.model
  model$s <- as.integer(candidate_row$s[[1L]])
  model$ranks <- .candidate_cached_ranks(object)
  model
}

.candidate_attach_common_view_fields <- function(out, object, candidate_row, W = NULL) {
  id <- as.character(candidate_row$candidate_id[[1L]])
  out$tuning_grid <- NULL
  out$meta$s <- as.integer(candidate_row$s[[1L]])
  out$meta$candidate_id <- id
  tuning_payload <- .candidate_tuning_payload(object, candidate_row)
  out$meta$tuning <- tuning_payload
  if ("tuning" %in% names(out)) {
    out$tuning <- tuning_payload
  }
  out$.model <- .candidate_update_model(object, candidate_row)

  keep_weights <- any(c("W", "ess", "neighbors") %in% names(object))
  facts <- .dnn_group_facts(object$.model$groups)
  # Untied table-based candidates need only a scalar table norm. Tied targets
  # require their actual averaged weights; reuse W if the adapter already has it.
  concentration <- NULL
  if (is.null(W) && !keep_weights && !any(facts$tied_targets)) {
    tab <- .rank_weight_table(.candidate_rank_weight_spec(object, candidate_row), object$meta$n)
    if (!is.null(tab)) concentration <- rep(1 / sum(tab^2), object$meta$m)
  }
  if (is.null(concentration)) {
    if (is.null(W)) W <- .candidate_weights_from_ranks(object, candidate_row)
    concentration <- 1 / colSums(W^2)
  }
  out$meta <- .dnn_attach_notice_facts(out$meta, object$.model$groups, concentration, facts)
  if (keep_weights) {
    if ("W" %in% names(object)) {
      out$W <- W
    }
    if ("ess" %in% names(object)) {
      out$ess <- concentration
    }
    if ("neighbors" %in% names(object)) {
      k_eff <- max(1L, min(nrow(W), as.integer(out$meta$k)))
      ranks <- .candidate_cached_ranks(object)
      out$neighbors <- apply(ranks, 2L, function(r) order(r)[seq_len(k_eff)])
      if (is.null(dim(out$neighbors))) {
        out$neighbors <- matrix(out$neighbors, nrow = k_eff)
      }
    }
  }
  out
}

.candidate_se_vector <- function(rows) {
  se <- as.numeric(rows$se)
  if (all(is.na(se))) {
    NULL
  } else {
    se
  }
}

.candidate_value_matrix <- function(rows, value_col, template) {
  component <- rownames(template)
  if (is.null(component)) {
    component <- unique(as.character(rows$component))
  }
  target_id <- seq_len(ncol(template))
  out <- matrix(NA_real_, nrow = length(component), ncol = length(target_id))
  rownames(out) <- component
  colnames(out) <- colnames(template)

  for (i in seq_len(nrow(rows))) {
    comp_idx <- match(as.character(rows$component[[i]]), component)
    target_idx <- match(as.integer(rows$target_id[[i]]), target_id)
    if (!is.na(comp_idx) && !is.na(target_idx)) {
      out[comp_idx, target_idx] <- rows[[value_col]][[i]]
    }
  }
  out
}

.candidate_np_view <- function(object, candidate) {
  row <- .candidate_resolve_row(object, candidate)
  rows <- .candidate_estimates(object, row)
  out <- object
  out$mu <- as.numeric(rows$estimate)
  out["se"] <- list(.candidate_se_vector(rows))
  .candidate_attach_common_view_fields(out, object, row)
}

.candidate_cate_view <- function(object, candidate) {
  row <- .candidate_resolve_row(object, candidate)
  rows <- .candidate_estimates(object, row)
  out <- object
  out$meta <- .cate_result_meta(out$meta)
  out$tau <- as.numeric(rows$estimate)
  out["se"] <- list(.candidate_se_vector(rows))
  .candidate_attach_common_view_fields(out, object, row)
}

.candidate_affine_view <- function(object, candidate) {
  row <- .candidate_resolve_row(object, candidate)
  model <- .affine_prediction_state(object)
  ranks <- .candidate_cached_ranks(object)
  W <- .candidate_weights_from_ranks(object, row)
  fit <- .affine_ratio_evaluate(model$signal, W, ranks,
    .candidate_rank_weight_spec(object, row), !is.null(object$se),
    object$meta$jackknife, model$denominator_tol, groups = model$groups)
  out <- object
  for (field in setdiff(names(fit), "deletions")) out[field] <- fit[field]
  if ("deletions" %in% names(object)) out$deletions <- fit$deletions
  if ("pseudo_values" %in% names(object)) {
    pv <- fit$deletions
    for (j in seq_along(fit$theta)) {
      full <- matrix(fit$components[j, ], nrow(W), 2L, byrow = TRUE)
      pv[, , j] <- full + (nrow(W) - 1) * (full - pv[, , j])
    }
    out$pseudo_values <- pv
  }
  .candidate_attach_common_view_fields(out, object, row, W = W)
}

.candidate_plqr_view <- function(object, candidate) {
  row <- .candidate_resolve_row(object, candidate)
  rows <- .candidate_estimates(object, row)
  out <- object
  out$moment <- .candidate_value_matrix(rows, "estimate", object$moment)
  # Grid mode carries moment-level jackknife SEs only; direct coefficient
  # inference is unavailable with tuning_mode = "grid".
  out["moment_se"] <- list(if (all(is.na(rows$se))) {
    NULL
  } else {
    .candidate_value_matrix(rows, "se", object$moment)
  })
  .candidate_attach_common_view_fields(out, object, row)
}

.candidate_ids_for_list <- function(object, candidates) {
  grid <- .candidate_grid_payload(object)
  ids <- grid$candidates$candidate_id
  if (missing(candidates) || identical(candidates, "all")) {
    return(ids)
  }
  if (is.character(candidates)) {
    if (length(candidates) < 1L || anyNA(candidates) || any(!nzchar(candidates))) {
      .abort_bad_arg(
        "candidates",
        "must contain one or more non-missing candidate id strings, or be \"all\"",
        quoted = FALSE
      )
    }
    unknown <- setdiff(candidates, ids)
    if (length(unknown) > 0L) {
      .abort_bad_arg(
        "candidates",
        sprintf(
          "unknown candidate id(s): %s; use one of: %s",
          .dnn_message_examples(unknown)$text,
          .dnn_message_examples(ids)$text
        ),
        quoted = FALSE
      )
    }
    return(candidates)
  }
  if (is.numeric(candidates) || is.integer(candidates)) {
    if (length(candidates) < 1L) {
      .abort_bad_arg(
        "candidates",
        "must contain one or more integer row indices, or be \"all\"",
        quoted = FALSE
      )
    }
    return(vapply(
      candidates,
      function(i) as.character(.candidate_resolve_row(object, i)$candidate_id[[1L]]),
      character(1L)
    ))
  }
  .abort_bad_arg(
    "candidates",
    "must be \"all\", candidate id strings, or integer row indices",
    quoted = FALSE
  )
}

#' Extract one candidate view from a grid fit
#'
#' `candidate_fit()` turns one row of a grid-mode fit into a scalar-fit view.
#' The returned object has the same class as `object`, with primary estimates
#' and standard errors taken from `object$tuning_grid`, candidate-dependent
#' extras such as `W`, `ess`, and `neighbors` rebuilt from cached ranks when
#' the parent carried them, and candidate-independent state carried over. The
#' rebuilt `neighbors` extra is the first `k` nearest observations, not the
#' full DNN weight support.
#'
#' Reconstructed weights use cached ranks and target-specific distance groups,
#' averaging within exact-distance blocks as in [weights_matrix()].
#' Candidate uncertainty redistributes each deleted member's block mass over
#' the survivors using the `n-1` rank table. Calls with ties warn once: the
#' accompanying theory assumes a continuous distance distribution and does not
#' establish inference under persistent ties. During experimental development,
#' refit saved models after package updates.
#'
#' Candidate IDs identify grid rows, not unique values of `s`. Repeated resolved
#' orders remain separate candidates, and plots distinguish their IDs in labels.
#' Extraction validates retained controls, target IDs and cached ranks; an
#' incomplete or inconsistent saved object must be refitted.
#'
#' @section Experimental methods:
#' Methods for `rdistnn_plqr_result` objects are experimental, as is all PLQR
#' functionality in this release; see [plqr_dnn()].
#'
#' @param object A grid fit returned by [np_dnn()], [cate_dnn()], or
#'   [clate_dnn()], [affine_ratio_dnn()] or [plqr_dnn()] with
#'   `tuning_mode = "grid"`.
#' @param candidate A single `candidate_id` string or integer row index into
#'   `object$tuning_grid$candidates`.
#' @param ... Reserved for future methods.
#'
#' @return A scalar-style view with the same class as the parent. Affine-ratio
#'   views rebuild paired components and covariance from retained signals and
#'   frozen ranks; nuisance models and geometry are never refitted.
#'
#' @examples
#' set.seed(1)
#' x_mat <- matrix(rnorm(40 * 2), 40, 2)
#' y_vec <- x_mat[, 1] + rnorm(40, sd = 0.2)
#' fit <- np_dnn(
#'   X = x_mat, x0 = c(0, 0), Y = y_vec,
#'   tuning_mode = "grid", se = FALSE
#' )
#' candidate_fit(fit, fit$tuning_grid$candidates$candidate_id[[1L]])
#' @export
candidate_fit <- function(object, candidate, ...) {
  .dnn_reuse_scope(object)
  UseMethod("candidate_fit")
}

#' @rdname candidate_fit
#' @export
candidate_fit.default <- function(object, candidate, ...) {
  .candidate_grid_payload(object)
  .abort_bad_arg(
    "object",
    sprintf("candidate_fit() is not implemented for class '%s'", class(object)[1L]),
    quoted = FALSE
  )
}

#' @rdname candidate_fit
#' @export
candidate_fit.rdistnn_np_result <- function(object, candidate, ...) {
  .dnn_reuse_scope(object)
  .candidate_np_view(object, candidate)
}

#' @rdname candidate_fit
#' @export
candidate_fit.rdistnn_cate_result <- function(object, candidate, ...) {
  .dnn_reuse_scope(object)
  .candidate_cate_view(object, candidate)
}

#' @rdname candidate_fit
#' @export
candidate_fit.rdistnn_affine_ratio_result <- function(object, candidate, ...) {
  .dnn_reuse_scope(object)
  .candidate_affine_view(object, candidate)
}

#' @rdname candidate_fit
#' @export
candidate_fit.rdistnn_plqr_result <- function(object, candidate, ...) {
  .dnn_reuse_scope(object)
  .candidate_plqr_view(object, candidate)
}

#' Extract multiple candidate views from a grid fit
#'
#' @param candidates `"all"` (default), a character vector of candidate ids, or
#'   an integer vector of candidate row indices.
#'
#' @return A named list of scalar-style candidate views.
#' @rdname candidate_fit
#' @export
candidate_fits <- function(object, candidates = "all", ...) {
  .dnn_reuse_scope(object)
  UseMethod("candidate_fits")
}

.candidate_fits_impl <- function(object, candidates = "all", ...) {
  ids <- .candidate_ids_for_list(object, candidates)
  out <- lapply(ids, function(id) candidate_fit(object, id, ...))
  stats::setNames(out, ids)
}

#' @rdname candidate_fit
#' @export
candidate_fits.default <- function(object, candidates = "all", ...) {
  .candidate_grid_payload(object)
  .abort_bad_arg(
    "object",
    sprintf("candidate_fits() is not implemented for class '%s'", class(object)[1L]),
    quoted = FALSE
  )
}

#' @rdname candidate_fit
#' @export
candidate_fits.rdistnn_np_result <- function(object, candidates = "all", ...) {
  .dnn_reuse_scope(object)
  .candidate_fits_impl(object, candidates = candidates, ...)
}

#' @rdname candidate_fit
#' @export
candidate_fits.rdistnn_cate_result <- function(object, candidates = "all", ...) {
  .dnn_reuse_scope(object)
  .candidate_fits_impl(object, candidates = candidates, ...)
}

#' @rdname candidate_fit
#' @export
candidate_fits.rdistnn_affine_ratio_result <- function(object, candidates = "all", ...) {
  .dnn_reuse_scope(object)
  .candidate_fits_impl(object, candidates = candidates, ...)
}

#' @rdname candidate_fit
#' @export
candidate_fits.rdistnn_plqr_result <- function(object, candidates = "all", ...) {
  .dnn_reuse_scope(object)
  .candidate_fits_impl(object, candidates = candidates, ...)
}
