# Shared retained geometry for core DNN/TDNN estimators.

# Resolve the named metric and its setup once. Keep the unfitted specification
# separately for the established explicit-deletion (loo) reference path.
.dnn_fitted_geometry <- function(distance, x_mat, ctx) {
  spec <- .resolve_distance(distance)
  fitted_ctx <- ctx
  if (is.null(fitted_ctx) && is.function(spec$setup)) fitted_ctx <- spec$setup(x_mat)
  frozen <- spec
  frozen$setup <- NULL
  list(distance = frozen, ctx = fitted_ctx, inference_distance = spec, inference_ctx = ctx)
}
