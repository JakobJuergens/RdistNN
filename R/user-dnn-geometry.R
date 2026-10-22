# Shared retained geometry for core DNN/TDNN estimators.

# Resolve the named metric and its setup once. Keep the unfitted specification
# separately for the established explicit-deletion (loo) reference path.
.dnn_fitted_geometry <- function(distance, x_mat, ctx) {
  spec <- .resolve_distance(distance)
  fitted_ctx <- ctx
  if (is.null(fitted_ctx) && is.function(spec$setup)) fitted_ctx <- spec$setup(x_mat)
  frozen <- spec
  frozen$setup <- NULL
  list(distance = frozen, ctx = fitted_ctx, inference_distance = spec, inference_ctx = ctx,
       tie_policy = .dnn_tie_policy())
}

# Cached ranks/orders are authoritative. Reconstructed ranks use the current
# canonical covariate tie rule through the ordinary geometry resolver.
.dnn_retained_ranking <- function(object, allow_recompute = TRUE) {
  .dnn_check_tie_model(object$.model)
  model <- object$.model
  ranks <- if (is.list(model)) model$ranks else NULL
  if (is.null(ranks) && is.list(model)) ranks <- model$geometry$ranks
  if (is.null(ranks)) ranks <- object$ranks
  if (!is.null(ranks)) {
    ranks <- .rank_check_mat(ranks)
    ord <- matrix(apply(ranks, 2L, order), nrow(ranks), ncol(ranks))
    return(list(ranks = ranks, ord = ord, dx = object$dx))
  }
  if (!is.null(object$ord)) {
    ord <- .validate_permutation_matrix(object$ord, "object$ord")
    ranks <- matrix(apply(ord, 2L, ord_to_ranks), nrow(ord), ncol(ord))
    return(list(ranks = ranks, ord = ord, dx = object$dx))
  }
  if (!allow_recompute || !is.list(model) || is.null(model$x_mat) || is.null(object$x0)) {
    return(NULL)
  }
  ranking <- if (!is.null(object$dx)) list(dx = object$dx) else model$distance
  if (is.null(ranking)) return(NULL)
  .resolve_ranking(ranking, model$x_mat, object$x0, model$ctx)
}
