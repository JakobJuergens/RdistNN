# Retained-state and numerical adapters for DNN/TDNN regression.
# Public result methods consume these helpers without changing estimand classes.

# Resolve the named metric and its setup once. Keep the unfitted specification
# separately for the established explicit-deletion (loo) reference path.
.np_fitted_geometry <- function(distance, x_mat, ctx) {
  spec <- .resolve_distance(distance)
  fitted_ctx <- ctx
  if (is.null(fitted_ctx) && is.function(spec$setup)) fitted_ctx <- spec$setup(x_mat)
  frozen <- spec
  frozen$setup <- NULL
  list(distance = frozen, ctx = fitted_ctx, inference_distance = spec, inference_ctx = ctx)
}

.np_prediction_state <- function(object) {
  model <- object$.model
  if (!is.list(model) || is.null(model$x_mat) || is.null(model$y_vec) ||
      is.null(model$k) || is.null(model$kind)) {
    .abort_bad_state(paste0(
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
.np_centered_deletions <- function(model, x0_mat, ranks = NULL, jackknife = "fast") {
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
    if (is.null(ranks)) ranks <- .np_prediction_weights(model, x0_mat, "ranks")$ranks
    rw <- .np_rank_weight(model)
    loo <- .dnn_loo_centered(model$y_vec, ranks, .rank_weight_table(rw, n - 1L),
      tab_full = .rank_weight_table(rw, n))$loo
  }
  sweep(loo, 2L, colMeans(loo), "-")
}

.np_marginal_variance <- function(model, x0_mat, ranks = NULL, jackknife = "fast") {
  centered <- .np_centered_deletions(model, x0_mat, ranks, jackknife)
  n <- nrow(centered)
  (n - 1) / n * colSums(centered^2)
}

.np_joint_covariance <- function(model, x0_mat, ranks = NULL, jackknife = "fast") {
  centered <- .np_centered_deletions(model, x0_mat, ranks, jackknife)
  n <- nrow(centered)
  V <- (n - 1) / n * crossprod(centered)
  dimnames(V) <- rep(list(paste0("target_", seq_len(nrow(x0_mat)))), 2L)
  V
}
