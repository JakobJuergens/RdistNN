# Reporting-target resolution shared by NP, CATE and paired-ratio adapters.
# Adapters validate their own scientific state first; this boundary needs only
# canonical training covariates and the resolved tabular covariate identities.
.dnn_resolve_prediction_targets <- function(model, newdata = NULL, X0 = NULL, x0 = NULL) {
  if (!is.null(newdata) && (!is.null(X0) || !is.null(x0))) {
    .abort_bad_arg(
      "newdata/X0/x0",
      "must provide either `newdata` or `X0`/`x0`, not both",
      quoted = FALSE
    )
  }
  if (!is.null(newdata)) {
    x0_mat <- if (is.data.frame(newdata) && !is.null(model$x_cols)) {
      .dnn_col_mat(model$x_cols, data = newdata, arg = "x_cols", out_arg = "newdata")
    } else {
      newdata
    }
    x0_mat <- .validate_matrix(x0_mat, "newdata", min_rows = 1L,
      min_cols = 1L, numeric = TRUE, finite = TRUE)
    return(.validate_ncol_relation(model$x_mat, x0_mat,
      y_arg = "newdata", x_arg = "X", relation = "match"))
  }
  if (is.null(X0) && is.null(x0)) {
    .abort_bad_arg("newdata/X0/x0",
      "must provide `newdata`, `X0`, or a single `x0`", quoted = FALSE)
  }
  .dnn_resolve_targets(model$x_mat, X0 = X0, x0 = x0)
}
