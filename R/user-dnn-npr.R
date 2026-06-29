# DNN-based nonparametric regression wrapper (rank-weighted local mean)

#' Nonparametric regression via DNN rank weights
#'
#' Estimate a localized nonparametric mean by combining observed outcomes with
#' DNN-based local weights.
#'
#' @details
#' Two interfaces are supported.
#'
#' 1. Data-frame interface (recommended): pass `data`, `x_cols`, and either
#'    `y_col` or `y` (exactly one).
#' 2. Matrix interface: pass `X`, one of `X0` or `x0`, and `Y`.
#'
#' For `method = "exact"`, DNN weights are computed from the full rank
#' permutation at each target. In particular, for the paper-faithful `k = 1`
#' case the exact support extends through rank `n - s + 1`; this is not an
#' `s`-nearest-neighbor truncation. The `"fast"` jackknife option reuses the
#' full-sample rank matrix and updates ranks under deletion.
#'
#' @param data Optional data frame containing covariates and outcome columns.
#' @param x_cols Names or indices of covariate columns in `data`.
#' @param y_col Name or index of the outcome column in `data`.
#' @param y Optional numeric outcome vector of length `n`. In the data-frame
#'   interface, use this instead of `y_col`.
#' @param X Numeric matrix (`n x p`) for the matrix interface.
#' @param Y Numeric outcome vector of length `n` for the matrix interface.
#' @param X0 Numeric matrix (`m x p`) of target points; mutually exclusive with
#'   `x0`.
#' @param x0 Numeric vector of length `p` for a single target point.
#' @param distance Distance name or specification (for example `"l2"`).
#' @param method One of `"exact"` or `"asymp"` for DNN weights.
#' @param k,s Integer DNN tuning parameters.
#' @param tuning Optional kernel-order tuning specification. When supplied,
#'   `s` must be omitted and is resolved from `tuning`.
#' @param normalize Logical; renormalize columns of `W` to sum to 1.
#' @param se Logical; compute delete-1 jackknife standard errors.
#' @param jackknife One of `"fast"`, `"loo"`, or `"none"`. `"fast"` (default)
#'   uses cached full-sample ranks with delete-1 rank updates, `"loo"`
#'   recomputes leave-one-out weights directly, and `"none"` skips SEs.
#' @param ctx Optional context passed to distance and ranking components.
#' @param return Optional extras in addition to the always-returned `mu`, `se`,
#'   and `meta`. Allowed values:
#'   `c("W", "ranks", "ord", "dx", "y", "ess", "neighbors", "tuning")`.
#'
#' @return An object of class `np_result` with components:
#'   - `mu`: estimated local mean values at target points;
#'   - `se`: optional jackknife standard errors;
#'   - `meta`: metadata (`n`, `m`, `p`, `k`, `s`, `method`, `distance`);
#'   - optional extras requested via `return`.
#'
#' @seealso [weights_dnn()], [np_demo_data()].
#'
#' @examples
#' set.seed(1)
#' n <- 60
#' x_mat <- matrix(rnorm(n * 2), n, 2)
#' y_vec <- x_mat[, 1] - 0.5 * x_mat[, 2] + rnorm(n, sd = 0.3)
#' out <- np_dnn(
#'   X = x_mat, x0 = c(0, 0), Y = y_vec,
#'   distance = "l2", method = "exact", k = 1L, s = 3L,
#'   se = FALSE, return = c("W", "ess")
#' )
#' out$mu
#' @export
np_dnn <- function(
  data = NULL,
  x_cols = NULL,
  y_col = NULL,
  y = NULL,
  X = NULL,
  Y = NULL,
  X0 = NULL,
  x0 = NULL,
  distance = "l2",
  method = c("exact", "asymp"),
  k = 1L,
  s = 1L,
  tuning = NULL,
  normalize = TRUE,
  se = TRUE,
  jackknife = c("fast", "loo", "none"),
  ctx = NULL,
  return = c("W")
) {
  method <- match.arg(method)
  jackknife <- match.arg(jackknife)
  s_arg <- if (missing(s)) NULL else s
  k <- .validate_scalar_int(k, "k")
  .warn_unknown_return(return, valid = c("W", "ranks", "ord", "dx", "y", "ess", "neighbors", "tuning"))

  if (!is.null(data)) {
    if (!is.data.frame(data)) {
      .abort_bad_arg("data", "must be a data.frame")
    }
    if (is.null(x_cols)) {
      stop("With `data=`, `x_cols` is required.", call. = FALSE)
    }
    x_mat <- .dnn_col_mat(x_cols, data = data, arg = "x_cols", out_arg = "X")
    n <- nrow(x_mat)
    p <- ncol(x_mat)

    if (is.null(y_col) == is.null(y)) {
      stop("With `data=`, provide exactly one of `y_col` or `y`.", call. = FALSE)
    }
    if (!is.null(y_col)) {
      y_vec <- .dnn_col_vec(y_col, data, n, arg = "y_col", ref = "X")
    } else {
      y_vec <- .validate_vector_nrow(y, n, arg = "y", ref = "X", strict = TRUE)
    }
  } else {
    if (is.null(X)) {
      stop("Provide either `data` ... or matrix interface via `X`.", call. = FALSE)
    }
    x_mat <- .validate_matrix_nonempty(X, arg = "X")
    n <- nrow(x_mat)
    p <- ncol(x_mat)

    if (!is.null(y)) {
      stop("Matrix API uses `Y`; `y` is only supported with `data=`.", call. = FALSE)
    }
    if (is.null(Y)) {
      stop("Matrix API requires `Y`.", call. = FALSE)
    }
    y_vec <- .validate_vector_nrow(Y, n, arg = "Y", ref = "X", strict = TRUE)
  }

  s_resolved <- .dnn_resolve_s(
    s = s_arg,
    tuning = tuning,
    n = n,
    p = p,
    k = k,
    se = se,
    jackknife = jackknife,
    ctx = ctx
  )
  s <- s_resolved$s
  tuning_payload <- s_resolved$tuning
  tuning_used <- isTRUE(s_resolved$tuning_used)

  x0_mat <- .dnn_resolve_targets(x_mat, X0 = X0, x0 = x0)

  weight_return <- return
  if (isTRUE(se) && identical(jackknife, "fast")) {
    # Request ranks once so fast jackknife can reuse them.
    weight_return <- unique(c(return, "ranks"))
  }
  weight_fit <- .dnn_weights_with_extras(
    x_mat = x_mat,
    x0_mat = x0_mat,
    distance = distance,
    method = method,
    k = k,
    s = s,
    normalize = normalize,
    ctx = ctx,
    return = weight_return
  )
  wt <- weight_fit$wt
  weight_matrix <- weight_fit$W
  m_out <- weight_fit$m

  mu_hat <- as.numeric(.dnn_estimate_from_signal(y_vec, weight_matrix, n = n, arg = "Y")[1L, ])
  se_hat <- NULL
  if (isTRUE(se) && jackknife != "none") {
    se_hat <- .np_dnn_jackknife_se_loo(
      x_mat = x_mat,
      x0_mat = x0_mat,
      y = y_vec,
      distance = distance,
      method = method,
      k = k,
      s = s,
      normalize = normalize,
      ctx = ctx,
      jackknife = jackknife,
      precomputed_ranks = wt$ranks
    )
  }

  out <- list(
    mu = mu_hat,
    se = se_hat,
    meta = list(
      n = n, m = m_out, p = p,
      k = k, s = s,
      method = method, distance = .dnn_distance_name(distance)
    )
  )
  if (tuning_used) out$meta$tuning <- tuning_payload
  out <- .dnn_attach_common_extras(
    out,
    return = return,
    wt = wt,
    W = weight_matrix,
    k = k,
    tuning_payload = tuning_payload
  )
  if ("y" %in% return) out$y <- y_vec
  out$.model <- list(
    y_vec = y_vec
  )

  class(out) <- "np_result"
  out
}

.np_dnn_jackknife_se_loo <- function(
  x_mat, x0_mat, y, distance, method, k, s, normalize, ctx,
  jackknife = c("fast", "loo"),
  precomputed_ranks = NULL
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
      precomputed_ranks = precomputed_ranks
    )[1L, ]
  )
}
