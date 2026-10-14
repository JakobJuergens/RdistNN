# Per-candidate views over grid fits.

.candidate_grid_payload <- function(object) {
  grid <- object$tuning_grid
  if (is.null(grid)) {
    .abort_bad_arg(
      "object",
      "does not contain a `tuning_grid` payload; refit with `tuning_mode = \"grid\"` before calling candidate_fit()",
      quoted = FALSE
    )
  }
  if (!is.list(grid) || !is.data.frame(grid$candidates) ||
      !is.data.frame(grid$estimates) || !is.data.frame(grid$diagnostics)) {
    .abort_bad_state("`object$tuning_grid` is malformed.")
  }
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
          candidate,
          paste(candidates$candidate_id, collapse = ", ")
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
  row
}

.candidate_meta_int <- function(object, field) {
  value <- object$meta[[field, exact = TRUE]]
  if (is.null(value)) {
    return(NULL)
  }
  if (!is.numeric(value) || length(value) != 1L || !is.finite(value) || value != as.integer(value)) {
    .abort_bad_state(sprintf("`object$meta$%s` must be an integer scalar.", field))
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
    .abort_bad_state(sprintf(
      "`object$x0` row count must match `object$meta$m` (got nrow(object$x0)=%d, m=%d).",
      x0_m,
      meta_m
    ))
  }
  if (!is.null(meta_m)) {
    return(meta_m)
  }
  x0_m
}

.candidate_validate_estimate_rows <- function(object, rows, id) {
  required <- c("target_id", "component", "estimate", "se")
  missing <- setdiff(required, names(rows))
  if (length(missing) > 0L) {
    .abort_bad_state(sprintf(
      "`object$tuning_grid$estimates` is missing required column(s): %s.",
      paste(missing, collapse = ", ")
    ))
  }

  expected_m <- .candidate_expected_m(object)
  target_id <- suppressWarnings(as.integer(rows$target_id))
  if (length(target_id) != nrow(rows) || anyNA(target_id)) {
    .abort_bad_state(sprintf("Estimate rows for candidate '%s' must have integer `target_id` values.", id))
  }
  if (!is.null(expected_m)) {
    invalid <- target_id < 1L | target_id > expected_m
    if (any(invalid)) {
      .abort_bad_state(sprintf(
        "Estimate rows for candidate '%s' have `target_id` outside 1..%d.",
        id,
        expected_m
      ))
    }
    missing_target <- setdiff(seq_len(expected_m), unique(target_id))
    if (length(missing_target) > 0L) {
      .abort_bad_state(sprintf(
        "Estimate rows for candidate '%s' are missing target_id value(s): %s.",
        id,
        paste(missing_target, collapse = ", ")
      ))
    }
  }

  if (inherits(object, "np_result") || inherits(object, "cate_result")) {
    duplicated_target <- unique(target_id[duplicated(target_id)])
    if (length(duplicated_target) > 0L) {
      .abort_bad_state(sprintf(
        "Estimate rows for candidate '%s' contain duplicate target_id value(s): %s.",
        id,
        paste(duplicated_target, collapse = ", ")
      ))
    }
  }

  if (inherits(object, "plqr_result") && !is.null(object$moment)) {
    expected_components <- rownames(as.matrix(object$moment))
    if (!is.null(expected_components)) {
      component <- as.character(rows$component)
      if (anyNA(component) || any(!nzchar(component))) {
        .abort_bad_state(sprintf("Estimate rows for candidate '%s' must have non-empty `component` values.", id))
      }
      unknown_component <- setdiff(unique(component), expected_components)
      if (length(unknown_component) > 0L) {
        .abort_bad_state(sprintf(
          "Estimate rows for candidate '%s' contain unknown component value(s): %s.",
          id,
          paste(unknown_component, collapse = ", ")
        ))
      }
      pair <- paste(target_id, component, sep = "\r")
      duplicated_pair <- unique(pair[duplicated(pair)])
      if (length(duplicated_pair) > 0L) {
        .abort_bad_state(sprintf("Estimate rows for candidate '%s' contain duplicate target/component rows.", id))
      }
      if (!is.null(expected_m)) {
        expected_pair <- as.vector(outer(seq_len(expected_m), expected_components, paste, sep = "\r"))
        if (length(setdiff(expected_pair, pair)) > 0L) {
          .abort_bad_state(sprintf("Estimate rows for candidate '%s' are missing target/component rows.", id))
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
  rows <- rows[order(rows$target_id, rows$component), , drop = FALSE]
  if (nrow(rows) < 1L) {
    .abort_bad_state(sprintf("No estimates recorded for candidate '%s'.", id))
  }
  .candidate_validate_estimate_rows(object, rows, id)
  rows
}

.candidate_rank_weight_spec <- function(object, candidate_row) {
  method <- if (!is.null(object$meta$method)) as.character(object$meta$method) else "exact"
  normalize <- if (is.list(object$.model) && !is.null(object$.model$normalize)) {
    isTRUE(object$.model$normalize)
  } else {
    TRUE
  }
  k <- as.integer(object$meta$k)
  s <- as.integer(candidate_row$s[[1L]])
  if (identical(method, "exact")) {
    rank_weight_dnn_exact(k = k, s = s)
  } else {
    rank_weight_dnn_asymp(k = k, s = s, normalize = normalize)
  }
}

.candidate_validate_cached_ranks <- function(object, ranks) {
  expected_n <- .candidate_meta_int(object, "n")
  if (!is.null(expected_n) && nrow(ranks) != expected_n) {
    .abort_bad_state(sprintf(
      "Cached ranks row count must match `object$meta$n` (got nrow=%d, n=%d).",
      nrow(ranks),
      expected_n
    ))
  }
  expected_m <- .candidate_expected_m(object)
  if (!is.null(expected_m) && ncol(ranks) != expected_m) {
    .abort_bad_state(sprintf(
      "Cached ranks column count must match the number of target points (got ncol=%d, m=%d).",
      ncol(ranks),
      expected_m
    ))
  }
  ranks
}

.candidate_cached_ranks <- function(object) {
  ranks <- NULL
  if (is.list(object$.model) && !is.null(object$.model$ranks)) {
    ranks <- object$.model$ranks
  } else if (!is.null(object$ranks)) {
    ranks <- object$ranks
  }
  if (is.null(ranks)) {
    .abort_bad_arg(
      "object",
      "does not store cached ranks for candidate reconstruction; refit with `tuning_mode = \"grid\"`",
      quoted = FALSE
    )
  }
  .candidate_validate_cached_ranks(object, .rank_check_mat(ranks))
}

.candidate_weights_from_ranks <- function(object, candidate_row) {
  ranks <- .candidate_cached_ranks(object)
  spec <- .candidate_rank_weight_spec(object, candidate_row)
  tab <- .rank_weight_table(spec, nrow(ranks))
  ctx <- if (is.list(object$.model)) object$.model$ctx else NULL
  if (is.null(tab)) {
    rank_weight_eval_batch(spec, ranks, ctx = ctx)
  } else {
    matrix(tab[ranks], nrow = nrow(ranks), ncol = ncol(ranks))
  }
}

.candidate_tuning_payload <- function(object, candidate_row) {
  payload <- if (!is.null(object$meta$tuning)) {
    object$meta$tuning
  } else {
    object$tuning
  }
  if (is.null(payload) || !is.list(payload)) {
    return(NULL)
  }
  payload$default$s <- as.integer(candidate_row$s[[1L]])
  if (is.null(payload$default$meta) || !is.list(payload$default$meta)) {
    payload$default$meta <- list()
  }
  payload$default$meta$lambda <- if ("lambda" %in% names(candidate_row)) {
    suppressWarnings(as.numeric(candidate_row$lambda[[1L]]))
  } else {
    NA_real_
  }
  payload$default$meta$candidate_id <- as.character(candidate_row$candidate_id[[1L]])
  payload$default$feasibility <- list(
    k_le_s = as.integer(object$meta$k) <= as.integer(candidate_row$s[[1L]]),
    s_lt_n = as.integer(candidate_row$s[[1L]]) < as.integer(object$meta$n)
  )
  payload
}

.candidate_update_model <- function(object, candidate_row) {
  model <- object$.model
  if (!is.list(model)) {
    return(model)
  }
  model$s <- as.integer(candidate_row$s[[1L]])
  if (is.null(model$ranks)) {
    model$ranks <- .candidate_cached_ranks(object)
  }
  model
}

.candidate_attach_common_view_fields <- function(out, object, candidate_row) {
  id <- as.character(candidate_row$candidate_id[[1L]])
  out$tuning_grid <- NULL
  out$meta$s <- as.integer(candidate_row$s[[1L]])
  out$meta$candidate_id <- id
  tuning_payload <- .candidate_tuning_payload(object, candidate_row)
  if (!is.null(tuning_payload)) {
    out$meta$tuning <- tuning_payload
  }
  if ("tuning" %in% names(out) && !is.null(tuning_payload)) {
    out$tuning <- tuning_payload
  }
  out$.model <- .candidate_update_model(object, candidate_row)

  if (any(c("W", "ess", "neighbors") %in% names(object))) {
    W <- .candidate_weights_from_ranks(object, candidate_row)
    if ("W" %in% names(object)) {
      out$W <- W
    }
    if ("ess" %in% names(object)) {
      out$ess <- 1 / colSums(W^2)
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
  out$se <- .candidate_se_vector(rows)
  .candidate_attach_common_view_fields(out, object, row)
}

.candidate_cate_view <- function(object, candidate) {
  row <- .candidate_resolve_row(object, candidate)
  rows <- .candidate_estimates(object, row)
  out <- object
  out$tau <- as.numeric(rows$estimate)
  out$se <- .candidate_se_vector(rows)
  .candidate_attach_common_view_fields(out, object, row)
}

.candidate_plqr_view <- function(object, candidate) {
  row <- .candidate_resolve_row(object, candidate)
  rows <- .candidate_estimates(object, row)
  out <- object
  out$moment <- .candidate_value_matrix(rows, "estimate", object$moment)
  # Grid mode carries moment-level jackknife SEs only; direct coefficient
  # inference is unavailable with tuning_mode = "grid".
  out$moment_se <- if (all(is.na(rows$se))) {
    NULL
  } else {
    .candidate_value_matrix(rows, "se", object$moment)
  }
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
          paste(unknown, collapse = ", "),
          paste(ids, collapse = ", ")
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
#' @section Experimental methods:
#' Methods for `plqr_result` objects are experimental, as is all PLQR
#' functionality in this release; see [plqr_dnn()].
#'
#' @param object A grid fit returned by [np_dnn()], [cate_dnn()], or
#'   [plqr_dnn()] with
#'   `tuning_mode = "grid"`.
#' @param candidate A single `candidate_id` string or integer row index into
#'   `object$tuning_grid$candidates`.
#' @param ... Reserved for future methods.
#'
#' @return A scalar-style `np_result`, `cate_result`, or `plqr_result` view.
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
candidate_fit.np_result <- function(object, candidate, ...) {
  .candidate_np_view(object, candidate)
}

#' @rdname candidate_fit
#' @export
candidate_fit.cate_result <- function(object, candidate, ...) {
  .candidate_cate_view(object, candidate)
}

#' @rdname candidate_fit
#' @export
candidate_fit.plqr_result <- function(object, candidate, ...) {
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
candidate_fits.np_result <- function(object, candidates = "all", ...) {
  .candidate_fits_impl(object, candidates = candidates, ...)
}

#' @rdname candidate_fit
#' @export
candidate_fits.cate_result <- function(object, candidates = "all", ...) {
  .candidate_fits_impl(object, candidates = candidates, ...)
}

#' @rdname candidate_fit
#' @export
candidate_fits.plqr_result <- function(object, candidates = "all", ...) {
  .candidate_fits_impl(object, candidates = candidates, ...)
}
