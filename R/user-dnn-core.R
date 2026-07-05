# Shared internal helpers for user-level DNN estimators.

.dnn_resolve_data_cols <- function(cols, data, arg, allow_multiple = TRUE) {
  if (is.null(cols)) {
    .abort_bad_arg(arg, "must be provided")
  }

  if (is.character(cols)) {
    if (length(cols) == 0L || anyNA(cols) || any(!nzchar(cols))) {
      .abort_bad_arg(arg, "must be a non-empty character vector of column names")
    }
    missing <- setdiff(cols, names(data))
    if (length(missing) > 0L) {
      stop(sprintf(
        "Column(s) %s (specified by `%s`) not found in `data`.\nAvailable columns: %s",
        paste(sprintf("'%s'", missing), collapse = ", "),
        arg,
        paste(names(data), collapse = ", ")
      ), call. = FALSE)
    }
    resolved <- cols
  } else if (is.numeric(cols)) {
    if (length(cols) == 0L) {
      .abort_bad_arg(arg, "must be a non-empty integer vector of column indices")
    }
    if (.has_na_only(cols)) {
      .abort_bad_arg(arg, "cannot contain `NA` indices")
    }
    if (.has_nonfinite_non_na(cols)) {
      .abort_bad_arg(arg, "cannot contain non-finite indices (`Inf`, `-Inf`, or `NaN`)")
    }
    if (any(cols != as.integer(cols))) {
      .abort_bad_arg(arg, "must be a non-empty integer vector of column indices")
    }
    idx <- as.integer(cols)
    p <- ncol(data)
    if (any(idx < 1L | idx > p)) {
      .abort_bad_arg(arg, sprintf("contains out-of-range indices; valid range is 1..%d", p))
    }
    resolved <- names(data)[idx]
  } else {
    .abort_bad_arg(arg, "must be a character vector of column names or an integer vector of indices")
  }

  if (!isTRUE(allow_multiple) && length(resolved) != 1L) {
    .abort_bad_arg(arg, "must reference exactly one column")
  }
  resolved
}

# Extract and validate a single data-frame column by name or index.
.dnn_col_vec <- function(col, data, n, arg, ref = "data") {
  col_name <- .dnn_resolve_data_cols(col, data = data, arg = arg, allow_multiple = FALSE)
  .validate_vector_nrow(data[[col_name]], n, arg = arg, ref = ref, strict = TRUE)
}

# Extract and validate multiple data-frame columns by names/indices.
.dnn_col_mat <- function(cols, data, arg, out_arg = "X") {
  col_names <- .dnn_resolve_data_cols(cols, data = data, arg = arg, allow_multiple = TRUE)
  .validate_matrix_nonempty(data[, col_names, drop = FALSE], arg = out_arg)
}

.dnn_resolve_targets <- function(x_mat, X0 = NULL, x0 = NULL) {
  if (!is.null(X0) && !is.null(x0)) {
    stop("Provide either `X0` or `x0`, not both.", call. = FALSE)
  }
  if (is.null(X0)) {
    if (is.null(x0)) {
      stop("Provide `X0` (m x p) or a single `x0` (length p).", call. = FALSE)
    }
    x0_vec <- .validate_target_vector(x0, ncol(x_mat), x0_arg = "x0", x_arg = "X")
    return(matrix(x0_vec, nrow = 1L))
  }
  .validate_ncol_relation(x_mat, X0, y_arg = "X0", x_arg = "X", relation = "equal")
}

.dnn_x0_with_names <- function(x0_mat, x_mat) {
  # Result objects store the evaluation points so downstream consumers
  # (plotting, reporting) can reconstruct target profiles. Propagate
  # covariate names from the training matrix when the targets have none.
  if (is.null(colnames(x0_mat)) && !is.null(colnames(x_mat)) &&
      ncol(x0_mat) == ncol(x_mat)) {
    colnames(x0_mat) <- colnames(x_mat)
  }
  x0_mat
}

.dnn_candidate_format_kappa <- function(kappa) {
  label <- vapply(
    kappa,
    function(x) format(x, trim = TRUE, scientific = FALSE),
    character(1L)
  )
  label <- gsub("\\.", "p", label)
  gsub("-", "m", label)
}

.dnn_candidate_unique_ids <- function(base_ids) {
  seen <- new.env(parent = emptyenv())
  ids <- character(length(base_ids))
  for (i in seq_along(base_ids)) {
    base <- base_ids[[i]]
    count <- if (exists(base, envir = seen, inherits = FALSE)) {
      get(base, envir = seen, inherits = FALSE) + 1L
    } else {
      1L
    }
    assign(base, count, envir = seen)
    ids[[i]] <- if (count == 1L) base else sprintf("%s__%d", base, count)
  }
  ids
}

.dnn_candidate_default_index <- function(candidates, default) {
  idx_s <- which(as.integer(candidates$s) == as.integer(default$s))
  if (length(idx_s) == 0L) {
    return(integer(0L))
  }

  default_kappa <- default$meta[["kappa", exact = TRUE]]
  default_kappa_num <- suppressWarnings(as.numeric(default_kappa))
  if (!is.null(default_kappa) && length(default_kappa_num) == 1L &&
      !is.na(default_kappa_num) && is.finite(default_kappa_num)) {
    kappa <- suppressWarnings(as.numeric(as.character(candidates$kappa)))
    tol <- sqrt(.Machine$double.eps) * pmax(1, abs(default_kappa_num))
    idx_kappa <- idx_s[!is.na(kappa[idx_s]) & abs(kappa[idx_s] - default_kappa_num) <= tol]
    if (length(idx_kappa) > 0L) {
      return(idx_kappa[1L])
    }
  }

  idx_s[1L]
}

.dnn_prepend_synthetic_default_candidate <- function(candidates, default) {
  synthetic <- candidates[NA_integer_, , drop = FALSE]
  rownames(synthetic) <- NULL
  synthetic$s <- as.integer(default$s)
  synthetic$kappa <- NA_real_
  synthetic$feasible <- TRUE
  synthetic$reason <- "ok"
  rbind(synthetic, candidates)
}

.dnn_resolve_candidate_schema <- function(grid, default) {
  if (!is.list(grid) || !is.data.frame(grid$candidates)) {
    .abort_bad_state("`tuning$grid$candidates` must be a data.frame.")
  }
  candidates <- as.data.frame(grid$candidates, stringsAsFactors = FALSE)
  if (nrow(candidates) < 1L) {
    .abort_bad_state("`tuning$grid$candidates` must contain at least one row.")
  }
  required <- c("s", "feasible", "reason")
  missing <- setdiff(required, names(candidates))
  if (length(missing) > 0L) {
    .abort_bad_state(sprintf(
      "`tuning$grid$candidates` is missing required column(s): %s.",
      paste(missing, collapse = ", ")
    ))
  }

  candidates$s <- as.integer(candidates$s)
  candidates$feasible <- as.logical(candidates$feasible)
  candidates$reason <- as.character(candidates$reason)
  if (!("kappa" %in% names(candidates))) {
    candidates$kappa <- NA_real_
  }
  candidates$kappa <- suppressWarnings(as.numeric(as.character(candidates$kappa)))

  default_idx <- .dnn_candidate_default_index(candidates, default)
  if (length(default_idx) == 0L) {
    candidates <- .dnn_prepend_synthetic_default_candidate(candidates, default)
    default_idx <- 1L
  }

  candidates$is_default <- FALSE
  candidates$is_default[default_idx[1L]] <- TRUE
  candidates$status <- ifelse(candidates$feasible, "ok", "infeasible")
  candidates$computed_key <- sprintf("s_%d", as.integer(candidates$s))

  base_ids <- sprintf("s_%d", as.integer(candidates$s))
  has_kappa <- !is.na(candidates$kappa)
  base_ids[has_kappa] <- paste0(
    "kappa_",
    .dnn_candidate_format_kappa(candidates$kappa[has_kappa])
  )
  candidates$candidate_id <- .dnn_candidate_unique_ids(base_ids)

  schema_cols <- c(
    "candidate_id", "computed_key", "is_default", "feasible",
    "reason", "status", "kappa", "s"
  )
  candidates <- candidates[, c(schema_cols, setdiff(names(candidates), schema_cols)), drop = FALSE]
  rownames(candidates) <- NULL
  grid$candidates <- candidates
  grid
}

.dnn_resolve_s <- function(s, tuning, n, p, k, se, jackknife, ctx = NULL) {
  needs_jk <- isTRUE(se) && !identical(jackknife, "none")

  tuning_source <- "tuning"
  if (is.null(tuning) && is.null(s)) {
    # Neither `s` nor `tuning` supplied: default to the paper undersmoothing
    # rule s = floor(n^((kx + 1)/(kx + 2))) via tuning_dnn_paper(). A fixed
    # constant default (for example s = 1) yields uniform weights and hence
    # no localization at all, so it must never be the silent default.
    tuning <- tuning_dnn_paper()
    tuning_source <- "paper_default"
  }

  if (is.null(tuning)) {
    s_out <- .validate_scalar_int(s, "s")
    k_ok <- as.integer(k) <= as.integer(s_out)
    s_ok <- if (needs_jk) as.integer(s_out) < as.integer(n) else TRUE
    if (!k_ok) {
      .abort_bad_arg("k", sprintf("must satisfy k <= s (got k=%d, s=%d)", as.integer(k), as.integer(s_out)))
    }
    if (!s_ok) {
      .abort_bad_arg(
        "s",
        "must satisfy s < n when delete-1 jackknife is requested",
        expected = sprintf("s < n=%d", as.integer(n)),
        actual = sprintf("s=%d (with `se`=%s, `jackknife`='%s')", as.integer(s_out), as.character(se), jackknife)
      )
    }
    default <- list(
      s = as.integer(s_out),
      meta = list(name = "manual", n = as.integer(n), kx = as.integer(p)),
      feasibility = list(k_le_s = k_ok, s_lt_n = s_ok)
    )
    grid <- .dnn_resolve_candidate_schema(
      list(
        candidates = data.frame(
          kappa = NA_real_,
          s = as.integer(s_out),
          feasible = TRUE,
          reason = "ok",
          stringsAsFactors = FALSE
        ),
        meta = list(name = "manual", n = as.integer(n), kx = as.integer(p))
      ),
      default = default
    )
    payload <- list(
      source = "manual",
      default = default,
      grid = grid
    )
    return(list(s = as.integer(s_out), tuning = payload, tuning_used = FALSE))
  }

  if (!is.null(s)) {
    .abort_bad_arg(
      "s/tuning",
      "cannot provide both explicit `s` and `tuning`; omit `s` when `tuning` is used",
      quoted = FALSE
    )
  }

  default <- tuning_eval(
    tuning = tuning,
    n = n,
    kx = p,
    k = k,
    se = se,
    jackknife = jackknife,
    ctx = ctx,
    mode = "default"
  )
  grid <- tuning_eval(
    tuning = tuning,
    n = n,
    kx = p,
    k = k,
    se = se,
    jackknife = jackknife,
    ctx = ctx,
    mode = "grid"
  )

  if (!is.list(default$feasibility)) {
    .abort_bad_state("`tuning_eval(..., mode='default')` must return `feasibility`.")
  }

  s_out <- .validate_scalar_int(default$s, "s")
  if (!isTRUE(default$feasibility$k_le_s)) {
    .abort_bad_arg("k", sprintf("must satisfy k <= s (got k=%d, s=%d)", as.integer(k), as.integer(s_out)))
  }
  if (needs_jk && !isTRUE(default$feasibility$s_lt_n)) {
    .abort_bad_arg(
      "s",
      "must satisfy s < n when delete-1 jackknife is requested",
      expected = sprintf("s < n=%d", as.integer(n)),
      actual = sprintf("s=%d (with `se`=%s, `jackknife`='%s')", as.integer(s_out), as.character(se), jackknife)
    )
  }

  payload <- list(
    source = tuning_source,
    default = default,
    grid = .dnn_resolve_candidate_schema(grid, default = default)
  )
  list(s = as.integer(s_out), tuning = payload, tuning_used = TRUE)
}

.dnn_grid_candidate_s <- function(candidate_row) {
  if (!("s" %in% names(candidate_row))) {
    return(NA_integer_)
  }
  as.integer(candidate_row$s[[1L]])
}

.dnn_grid_candidate_kappa <- function(candidate_row) {
  if (!("kappa" %in% names(candidate_row))) {
    return(NA_real_)
  }
  suppressWarnings(as.numeric(candidate_row$kappa[[1L]]))
}

.dnn_grid_top_k <- function(ctx) {
  if (is.list(ctx) && !is.null(ctx$k)) {
    return(.validate_scalar_int(ctx$k, "ctx$k"))
  }
  if (is.list(ctx) && is.list(ctx$meta) && !is.null(ctx$meta$k)) {
    return(.validate_scalar_int(ctx$meta$k, "ctx$meta$k"))
  }
  NA_integer_
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
    .abort_bad_state("`estimate_fn()` must return a data.frame.")
  }
  required <- c("target_id", "component", "estimate")
  missing <- setdiff(required, names(estimates))
  if (length(missing) > 0L) {
    .abort_bad_state(sprintf(
      "`estimate_fn()` output is missing required column(s): %s.",
      paste(missing, collapse = ", ")
    ))
  }
  if (nrow(estimates) < 1L) {
    .abort_bad_state("`estimate_fn()` output must contain at least one row.")
  }
  if (!is.numeric(estimates$estimate)) {
    .abort_bad_state("`estimate_fn()` column `estimate` must be numeric.")
  }
  estimates <- as.data.frame(estimates, stringsAsFactors = FALSE)
  estimates$estimate <- as.numeric(estimates$estimate)
  estimates
}

.dnn_grid_check_se <- function(se) {
  if (!is.data.frame(se)) {
    .abort_bad_state("`se_fn()` must return a data.frame.")
  }
  required <- c("target_id", "component", "se")
  missing <- setdiff(required, names(se))
  if (length(missing) > 0L) {
    .abort_bad_state(sprintf(
      "`se_fn()` output is missing required column(s): %s.",
      paste(missing, collapse = ", ")
    ))
  }
  if (!is.numeric(se$se)) {
    .abort_bad_state("`se_fn()` column `se` must be numeric.")
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
  if (is.null(tab)) {
    return(data.frame(
      ess = NA_real_,
      max_weight = NA_real_,
      top_share = NA_real_,
      degenerate = NA,
      stringsAsFactors = FALSE
    ))
  }

  abs_tab <- abs(tab)
  denom <- sum(abs_tab)
  ess <- 1 / sum(tab^2)
  k_top <- .dnn_grid_top_k(ctx)
  top_share <- NA_real_
  if (!is.na(k_top) && denom > 0) {
    top_n <- min(length(tab), as.integer(k_top))
    top_share <- sum(abs_tab[seq_len(top_n)]) / denom
  }
  data.frame(
    ess = ess,
    max_weight = max(abs_tab),
    top_share = top_share,
    degenerate = ess >= 0.9 * n,
    stringsAsFactors = FALSE
  )
}

.dnn_grid_compute_key <- function(
  candidate_row, ranks, weight_spec_fn, estimate_fn, se_fn, ctx
) {
  n <- nrow(ranks)
  spec <- weight_spec_fn(candidate_row)
  if (!is_rank_weight(spec)) {
    .abort_bad_state("`weight_spec_fn()` must return a rank_weight.")
  }

  tab <- .rank_weight_table(spec, n)
  W <- if (is.null(tab)) {
    rank_weight_eval_batch(spec, ranks, ctx = ctx)
  } else {
    matrix(tab[ranks], nrow = n, ncol = ncol(ranks))
  }

  estimates <- .dnn_grid_check_estimates(estimate_fn(W, candidate_row))
  se <- NULL
  if (!is.null(se_fn)) {
    tab_jk <- .rank_weight_table(spec, n - 1L)
    se <- .dnn_grid_check_se(se_fn(candidate_row, tab_jk))
  }

  list(
    estimates = .dnn_grid_merge_estimate_se(estimates, se),
    diagnostics = .dnn_grid_diagnostics_from_table(tab, n = n, ctx = ctx)
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
    kappa = .dnn_grid_candidate_kappa(candidate_row),
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

.dnn_grid_empty_diagnostics <- function() {
  data.frame(
    ess = NA_real_,
    max_weight = NA_real_,
    top_share = NA_real_,
    degenerate = NA,
    stringsAsFactors = FALSE
  )
}

.dnn_grid_attach_diagnostic_candidate <- function(diagnostic, candidate_row) {
  data.frame(
    candidate_id = as.character(candidate_row$candidate_id[[1L]]),
    computed_key = as.character(candidate_row$computed_key[[1L]]),
    kappa = .dnn_grid_candidate_kappa(candidate_row),
    s = .dnn_grid_candidate_s(candidate_row),
    ess = diagnostic$ess[[1L]],
    max_weight = diagnostic$max_weight[[1L]],
    top_share = diagnostic$top_share[[1L]],
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

.dnn_grid_eval <- function(
  candidates,
  ranks,
  weight_spec_fn,
  estimate_fn,
  se_fn = NULL,
  ctx = NULL
) {
  candidates <- .dnn_grid_check_candidates(candidates)
  ranks <- .rank_check_mat(ranks)
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
  if (sum(candidates$feasible) > 1L && length(feasible_keys) == 1L) {
    warning(
      "All feasible tuning candidates share the same computed key; the grid is effectively a single value.",
      call. = FALSE
    )
  }

  results <- setNames(vector("list", length(feasible_keys)), feasible_keys)
  for (key in feasible_keys) {
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
        ctx = ctx
      ),
      error = function(err) {
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
    base_diag <- if (ok) result$diagnostics else .dnn_grid_empty_diagnostics()
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

.dnn_distance_name <- function(distance) {
  if (is.character(distance)) {
    return(distance)
  }
  if (is.list(distance) && !is.null(distance$name)) {
    return(distance$name)
  }
  "user_distance"
}

.dnn_warn_degenerate_weights <- function(W) {
  n <- nrow(W)
  if (is.null(n) || n < 2L) {
    return(invisible(NULL))
  }
  ess_max <- max(1 / colSums(W^2))
  if (ess_max >= 0.9 * n) {
    warning(
      sprintf(
        paste0(
          "DNN weights are (near-)uniform at one or more target points ",
          "(max effective sample size %.1f of n = %d), so the fit provides ",
          "little or no localization. This typically means `s` is too small ",
          "relative to `k` (for example s = 1, or k = s). Increase `s`, or ",
          "omit both `s` and `tuning` to use the default paper tuning rule."
        ),
        ess_max, n
      ),
      call. = FALSE
    )
  }
  invisible(NULL)
}

.dnn_weights_with_extras <- function(
  x_mat, x0_mat, distance, method, k, s, normalize, ctx, return
) {
  keep_extras <- intersect(return, c("ranks", "ord", "dx"))
  wt <- weights_dnn(
    X = x_mat,
    X0 = x0_mat,
    distance = distance,
    k = as.integer(k),
    s = as.integer(s),
    method = method,
    normalize = normalize,
    ctx = ctx,
    return = unique(c("W", keep_extras))
  )
  .dnn_warn_degenerate_weights(wt$W)
  list(
    wt = wt,
    W = wt$W,
    m = ncol(wt$W)
  )
}

.dnn_signal_matrix <- function(signal, n, arg = "signal") {
  if (is.null(dim(signal))) {
    signal <- .validate_vector_nrow(signal, n, arg = arg, ref = "X", finite = TRUE)
    return(matrix(as.numeric(signal), ncol = 1L))
  }

  sig <- as.matrix(signal)
  if (nrow(sig) != n) {
    .abort_bad_state(sprintf("nrow(%s) must equal nrow(X).", arg))
  }
  if (!is.numeric(sig)) {
    .abort_bad_arg(arg, "must be numeric")
  }
  .abort_if_na_or_nonfinite(sig, arg = arg)
  sig
}

.dnn_estimate_from_signal <- function(signal, W, n, arg = "signal") {
  sig <- .dnn_signal_matrix(signal, n = n, arg = arg)
  crossprod(sig, W)
}

.dnn_attach_common_extras <- function(out, return, wt, W, k, tuning_payload = NULL) {
  if ("W" %in% return) out$W <- W
  if ("ranks" %in% return && !is.null(wt$ranks)) out$ranks <- wt$ranks
  if ("ord" %in% return && !is.null(wt$ord)) out$ord <- wt$ord
  if ("dx" %in% return && !is.null(wt$dx)) out$dx <- wt$dx
  if ("tuning" %in% return && !is.null(tuning_payload)) out$tuning <- tuning_payload
  if ("ess" %in% return) out$ess <- 1 / colSums(W^2)
  if ("neighbors" %in% return) {
    k_eff <- max(1L, min(nrow(W), as.integer(k)))
    if (!is.null(wt$ord)) {
      out$neighbors <- wt$ord[seq_len(k_eff), , drop = FALSE]
    } else if (!is.null(wt$ranks)) {
      out$neighbors <- apply(wt$ranks, 2L, function(r) order(r)[seq_len(k_eff)])
    }
  }
  out
}

.dnn_jackknife_se_from_loo <- function(loo, q, n, m) {
  if (q == 1L) {
    loo_bar <- colMeans(loo)
    centered <- sweep(loo, 2L, loo_bar, FUN = "-")
    var_jk <- ((n - 1) / n) * colSums(centered^2)
    return(matrix(sqrt(pmax(var_jk, 0)), nrow = 1L, ncol = m))
  }

  se <- matrix(NA_real_, nrow = q, ncol = m)
  for (j in seq_len(q)) {
    mj <- matrix(loo[, j, ], nrow = n, ncol = m)
    mj_bar <- colMeans(mj)
    centered <- sweep(mj, 2L, mj_bar, FUN = "-")
    var_jk <- ((n - 1) / n) * colSums(centered^2)
    se[j, ] <- sqrt(pmax(var_jk, 0))
  }
  se
}

.dnn_loo_prefix <- function(sig, ranks, tab_jk, tab_full = NULL) {
  ranks <- .rank_check_mat(ranks)
  n <- nrow(ranks)
  m <- ncol(ranks)
  if (n < 2L) {
    .abort_bad_arg("ranks", "must have at least two rows for delete-1 jackknife")
  }

  sig <- .dnn_signal_matrix(sig, n = n, arg = "sig")
  q <- ncol(sig)

  if (!is.numeric(tab_jk)) {
    .abort_bad_arg("tab_jk", "must be numeric", quoted = FALSE)
  }
  if (length(tab_jk) != n - 1L) {
    .abort_bad_arg(
      "tab_jk",
      "length must equal nrow(ranks) - 1",
      expected = sprintf("length %d", n - 1L),
      actual = sprintf("length %d", length(tab_jk)),
      quoted = FALSE
    )
  }
  .abort_if_na_or_nonfinite(tab_jk, arg = "tab_jk", quoted = FALSE)
  tab_jk <- as.numeric(tab_jk)

  if (!is.null(tab_full)) {
    if (!is.numeric(tab_full)) {
      .abort_bad_arg("tab_full", "must be numeric", quoted = FALSE)
    }
    if (length(tab_full) != n) {
      .abort_bad_arg(
        "tab_full",
        "length must equal nrow(ranks)",
        expected = sprintf("length %d", n),
        actual = sprintf("length %d", length(tab_full)),
        quoted = FALSE
      )
    }
    .abort_if_na_or_nonfinite(tab_full, arg = "tab_full", quoted = FALSE)

    res <- jackknife_rank_full(ranks, sig, as.numeric(tab_full), tab_jk)
    loo <- res$theta_loo
    if (q == 1L) {
      loo <- matrix(loo, nrow = n, ncol = m)
    }
    return(loo)
  }

  if (q == 1L) {
    loo <- matrix(NA_real_, nrow = n, ncol = m)
  } else {
    loo <- array(NA_real_, dim = c(n, q, m))
  }

  tab_left <- c(tab_jk, 0)
  tab_right <- c(0, tab_jk)

  for (j in seq_len(m)) {
    r <- ranks[, j]
    ord <- order(r)
    for (component in seq_len(q)) {
      sig_by_rank <- sig[ord, component]
      left_sum <- cumsum(sig_by_rank * tab_left)
      right_sum <- rev(cumsum(rev(sig_by_rank * tab_right)))
      loo_j <- c(0, left_sum)[r] + c(right_sum, 0)[r + 1L]

      if (q == 1L) {
        loo[, j] <- loo_j
      } else {
        loo[, component, j] <- loo_j
      }
    }
  }

  loo
}

.dnn_loo_loop <- function(sig, ranks, rw) {
  if (!is_rank_weight(rw)) {
    .abort_bad_arg("rw", "must be a rank_weight", quoted = FALSE)
  }
  ranks <- .rank_check_mat(ranks)
  n <- nrow(ranks)
  m <- ncol(ranks)
  if (n < 2L) {
    .abort_bad_arg("ranks", "must have at least two rows for delete-1 jackknife")
  }

  sig <- .dnn_signal_matrix(sig, n = n, arg = "sig")
  q <- ncol(sig)

  if (q == 1L) {
    loo <- matrix(NA_real_, nrow = n, ncol = m)
  } else {
    loo <- array(NA_real_, dim = c(n, q, m))
  }

  idx_all <- seq_len(n)
  for (i in idx_all) {
    keep <- (idx_all != i)

    r_i <- matrix(ranks[i, ], nrow = n - 1L, ncol = m, byrow = TRUE)
    rank_matrix_i <- ranks[keep, , drop = FALSE]
    rank_matrix_i <- rank_matrix_i - (rank_matrix_i > r_i)
    if (!is.integer(rank_matrix_i)) storage.mode(rank_matrix_i) <- "integer"

    W_i <- rank_weight_eval_batch(rw, rank_matrix_i)
    est_i <- crossprod(sig[keep, , drop = FALSE], W_i)
    if (q == 1L) {
      loo[i, ] <- as.numeric(est_i[1L, ])
    } else {
      loo[i, , ] <- est_i
    }
  }

  loo
}

.dnn_jackknife_se_loo_matrix <- function(
  x_mat, x0_mat, signal, distance, method, k, s, normalize, ctx,
  signal_arg = "signal",
  jackknife = c("fast", "loo"),
  precomputed_ranks = NULL
) {
  jackknife <- match.arg(jackknife)
  n <- nrow(x_mat)
  if (as.integer(s) >= n) {
    stop(sprintf(
      "Delete-1 jackknife requires s < n (got s=%d, n=%d). Use `se = FALSE` or reduce `s`.",
      as.integer(s), n
    ), call. = FALSE)
  }

  sig <- .dnn_signal_matrix(signal, n = n, arg = signal_arg)
  q <- ncol(sig)
  m <- nrow(x0_mat)

  if (identical(jackknife, "fast")) {
    ranks_full <- precomputed_ranks
    if (is.null(ranks_full)) {
      ranks_full <- weights_dnn(
        X = x_mat,
        X0 = x0_mat,
        distance = distance,
        k = k,
        s = s,
        method = method,
        normalize = normalize,
        ctx = ctx,
        return = "ranks"
      )$ranks
    }

    ranks_full <- as.matrix(ranks_full)
    if (nrow(ranks_full) != n || ncol(ranks_full) != m) {
      .abort_bad_state("`precomputed_ranks` must be an n x m ranks matrix.")
    }
    if (!is.integer(ranks_full)) storage.mode(ranks_full) <- "integer"

    rw <- if (identical(method, "exact")) {
      rank_weight_dnn_exact(k = as.integer(k), s = as.integer(s))
    } else {
      rank_weight_dnn_asymp(k = as.integer(k), s = as.integer(s), normalize = normalize)
    }

    tab_jk <- .rank_weight_table(rw, n - 1L)
    loo <- if (is.null(tab_jk)) {
      .dnn_loo_loop(sig, ranks_full, rw)
    } else {
      tab_full <- .rank_weight_table(rw, n)
      .dnn_loo_prefix(sig, ranks_full, tab_jk, tab_full = tab_full)
    }
    return(.dnn_jackknife_se_from_loo(loo = loo, q = q, n = n, m = m))
  }

  if (q == 1L) {
    loo <- matrix(NA_real_, nrow = n, ncol = m)
    for (i in seq_len(n)) {
      keep <- (seq_len(n) != i)
      wt_i <- weights_dnn(
        X = x_mat[keep, , drop = FALSE],
        X0 = x0_mat,
        distance = distance,
        k = k,
        s = s,
        method = method,
        normalize = normalize,
        ctx = ctx,
        return = "W"
      )
      loo[i, ] <- as.numeric(crossprod(sig[keep, 1L], wt_i$W))
    }
    return(.dnn_jackknife_se_from_loo(loo = loo, q = q, n = n, m = m))
  }

  loo <- array(NA_real_, dim = c(n, q, m))
  for (i in seq_len(n)) {
    keep <- (seq_len(n) != i)
    wt_i <- weights_dnn(
      X = x_mat[keep, , drop = FALSE],
      X0 = x0_mat,
      distance = distance,
      k = k,
      s = s,
      method = method,
      normalize = normalize,
      ctx = ctx,
      return = "W"
    )
    loo[i, , ] <- crossprod(sig[keep, , drop = FALSE], wt_i$W)
  }
  .dnn_jackknife_se_from_loo(loo = loo, q = q, n = n, m = m)
}
