# Retained-state and numerical adapters for DNN/TDNN regression.
# Public result methods consume these helpers without changing estimand classes.

.np_prediction_state <- function(object) {
  .dnn_validate_retained({
    model <- object$.model
    .dnn_check_tie_model(model)
    if (!is.list(model) || is.null(model$x_mat) || is.null(model$y_vec) ||
        is.null(model$k) || is.null(model$kind)) {
      .dnn_abort("state_invalid", paste0(
        "This object does not contain fitted model state for prediction. ",
        "Refit with the current np_dnn() or np_tdnn() implementation."))
    }
    model$x_mat <- .validate_matrix(model$x_mat, "object$.model$x_mat",
      min_rows = 1L, min_cols = 1L, numeric = TRUE, finite = TRUE)
    n <- nrow(model$x_mat)
    model$y_vec <- .validate_vector_nrow(model$y_vec, n, "object$.model$y_vec",
      ref = "object$.model$x_mat", strict = TRUE, finite = TRUE)
    model$k <- .validate_scalar_int(model$k, "object$.model$k")
    model$normalize <- .validate_scalar_logical(model$normalize, "object$.model$normalize")
    if (!is.character(model$method) || length(model$method) != 1L ||
        is.na(model$method) || !model$method %in% c("exact", "asymp")) {
      .abort_bad_arg("object$.model$method", "must be one of 'exact' or 'asymp'")
    }
    if (identical(model$kind, "tdnn")) {
      model$s1 <- .validate_scalar_int(model$s1, "object$.model$s1")
      model$s2 <- .validate_scalar_int(model$s2, "object$.model$s2")
      model$dim <- .validate_scalar_int(model$dim, "object$.model$dim")
      if (model$s1 >= model$s2 || model$s2 > n) {
        .abort_bad_arg("object$.model$s1/s2", sprintf("must satisfy s1 < s2 <= n=%d", n))
      }
    } else if (identical(model$kind, "dnn")) {
      model$s <- .validate_scalar_int(model$s, "object$.model$s")
      if (model$k > model$s || model$s > n) {
        .abort_bad_arg("object$.model$k/s", sprintf("must satisfy k <= s <= n=%d", n))
      }
    } else {
      .abort_bad_arg("object$.model$kind", "must be 'dnn' or 'tdnn'")
    }
    model$distance <- .resolve_distance(model$distance)
    model
  })
}

.np_rank_weight <- function(model) {
  if (identical(model$kind, "dnn")) {
    return(.dnn_rank_weight_spec(model$method, model$k, model$s, model$normalize))
  }
  if (identical(model$method, "exact")) {
    rank_weight_tdnn_exact(model$s1, model$s2, model$dim)
  } else {
    rank_weight_tdnn_asymp(model$s1, model$s2, model$dim, normalize = model$normalize)
  }
}

.np_prediction_weights <- function(model, x0_mat, return = "W", inference = FALSE) {
  distance <- if (inference) model$inference_distance else model$distance
  ctx <- if (inference) model$inference_ctx else model$ctx
  if (identical(model$kind, "tdnn")) {
    weights_tdnn(X = model$x_mat, X0 = x0_mat, distance = distance,
      method = model$method, s1 = model$s1, s2 = model$s2, dim = model$dim,
      normalize = model$normalize, ctx = ctx, return = return)
  } else {
    weights_dnn(X = model$x_mat, X0 = x0_mat, distance = distance,
      method = model$method, k = model$k, s = model$s,
      normalize = model$normalize, ctx = ctx, return = return)
  }
}

# Target columns share deletion row identities. Subtracting a target-specific
# common offset preserves their cross covariance and small local variation.
.np_covariance_deletions <- function(model, x0_mat, ranks = NULL, jackknife = "fast", groups = NULL) {
  n <- nrow(model$x_mat)
  order <- if (identical(model$kind, "tdnn")) model$s2 else model$s
  if (order >= n) {
    .abort_bad_arg("s", paste0(
      "must satisfy s < n (s2 < n for TDNN) for delete-1 covariance; ",
      "reduce the order or request point predictions only"))
  }
  if (identical(jackknife, "loo")) {
    loo <- .jackknife_recomputed_centered(matrix(model$y_vec, ncol = 1L), function(i) {
      deleted <- model
      deleted$x_mat <- model$x_mat[-i, , drop = FALSE]
      .np_prediction_weights(deleted, x0_mat, inference = TRUE)$W
    })
  } else {
    if (is.null(ranks) || is.null(groups)) {
      wt <- .np_prediction_weights(model, x0_mat, "ranks")
      if (is.null(ranks)) ranks <- wt$ranks
      groups <- wt$groups
    }
    rw <- .np_rank_weight(model)
    loo <- .dnn_loo_centered(model$y_vec, ranks, .rank_weight_table(rw, n - 1L),
      tab_full = .rank_weight_table(rw, n), groups = groups)$loo
  }
  loo
}

.np_marginal_variance <- function(model, x0_mat, ranks = NULL, jackknife = "fast", groups = NULL) {
  deletions <- .np_covariance_deletions(model, x0_mat, ranks, jackknife, groups)
  n <- nrow(deletions)
  .covariance_matrix(deletions, (n - 1) / n, diagonal = TRUE)
}

.np_joint_covariance <- function(model, x0_mat, ranks = NULL, jackknife = "fast", groups = NULL) {
  deletions <- .np_covariance_deletions(model, x0_mat, ranks, jackknife, groups)
  n <- nrow(deletions)
  V <- .covariance_matrix(deletions, (n - 1) / n)
  dimnames(V) <- rep(list(paste0("target_", seq_len(nrow(x0_mat)))), 2L)
  V
}

.np_dnn_jackknife_se_loo <- function(
  x_mat, x0_mat, y, distance, method, k, s, normalize, ctx,
  jackknife = c("fast", "loo"),
  precomputed_ranks = NULL, precomputed_groups = NULL
) {
  jackknife <- match.arg(jackknife)
  as.numeric(
    .dnn_jackknife_se_loo_matrix(
      x_mat = x_mat,
      x0_mat = x0_mat,
      signal = y,
      distance = distance,
      method = method,
      k = k,
      s = s,
      normalize = normalize,
      ctx = ctx,
      signal_arg = "Y",
      jackknife = jackknife,
      precomputed_ranks = precomputed_ranks,
      precomputed_groups = precomputed_groups
    )[1L, ]
  )
}
