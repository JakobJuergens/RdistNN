# DNN-based nonparametric regression wrapper (rank-weighted local mean)

#' Nonparametric regression via DNN rank weights
#'
#' Estimate a localized nonparametric mean by combining observed outcomes with
#' DNN-based local weights.
#'
#' @details
#' The documented sampling framework uses independent observations; cluster
#' and survey-design inference are not supplied. Observations with exactly
#' equal computed distance to a target share their block's average rank weight,
#' even when their covariates differ. Groups are target-specific and use exact
#' represented equality, without a tolerance. Ordering uses distance then input
#' row index; outcomes and nuisances do not define groups. Every jackknife
#' deletion averages the surviving blocks with the sample-size-n-1 rule.
#' A warning is emitted once per public call with distance ties: the accompanying
#' theory assumes a continuous distance distribution and does not establish
#' inference under persistent ties. Earlier retained fit schemas require refitting.
#'
#' Two interfaces are supported.
#'
#' 1. Data-frame interface (recommended): pass `data`, `x_cols`, and either
#'    `y_col` or `y` (exactly one).
#' 2. Matrix interface: pass `X`, one of `X0` or `x0`, and `Y`.
#'
#' With `data`, omit `X` and `Y`; use exactly one of `y_col` or `y`.
#' Column selectors require `data`. Conflicting routes are rejected.
#'
#' For `method = "exact"`, DNN weights are computed from the full rank
#' permutation at each target. The unaveraged rank-table support extends through rank
#' `n - s + k` (`n - s + 1` in the paper-faithful `k = 1` case); this is not
#' an `s`-nearest-neighbor truncation. The `"fast"` jackknife option reuses
#' the full-sample rank matrix and updates ranks under deletion.
#'
#' Averaging can extend nonzero observation weights through the end of a
#' distance block that straddles the nominal rank-table edge. Effective sample
#' sizes and weight diagnostics use the resulting observation weights.
#'
#' `tuning_mode = "default"` evaluates only the default kernel order and
#' preserves the ordinary scalar-fit return shape. `tuning_mode = "grid"`
#' evaluates the feasible candidates in `tuning`, or the paper default tuning
#' grid when `tuning` is omitted, and attaches a tidy `tuning_grid` payload.
#' The primary `mu`, `se`, `meta`, and optional extras such as `W`, `ranks`,
#' `ess`, and `neighbors` remain anchored to the default candidate. Full
#' per-candidate weight matrices are not stored. Candidate diagnostics include
#' `ess`, `max_weight`, `top_n`, `top5_weight_share`, and `degenerate`;
#' `top_n` is fixed at 5 and is independent of the estimator's `k`.
#'
#' @param data Optional data frame containing covariates and outcome columns.
#'   Column names must be non-missing, non-empty, and unique.
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
#' @param k Integer DNN neighbour count inside each size-`s` subsample. The
#'   paper theory uses `k = 1` (the default).
#' @param s Integer kernel order (undersmoothing parameter), or `NULL`. When
#'   both `s` and `tuning` are `NULL` (the default), `s` is resolved from the
#'   paper undersmoothing rule [tuning_dnn_paper()],
#'   \eqn{s = \lfloor 0.4 n^{(k_x+1)/(k_x+2)} \rfloor}. Note that small fixed
#'   values such as `s = 1` (or `k = s`) yield uniform weights and hence no
#'   localization; a warning is emitted in that case.
#' @param tuning Optional kernel-order tuning specification. When supplied,
#'   `s` must be omitted and is resolved from `tuning`.
#' @param tuning_mode One of `"default"` or `"grid"`. `"default"` evaluates
#'   the default `s` only. `"grid"` evaluates the feasible candidates from
#'   `tuning`, or from [tuning_dnn_paper()] when `tuning` is omitted, and
#'   stores candidate-level estimates, standard errors, diagnostics, and
#'   metadata in `tuning_grid`. Grid mode cannot be combined with explicit `s`.
#' @param normalize Logical; renormalize columns of `W` to sum to 1.
#' @param se Logical; compute delete-1 jackknife standard errors.
#' @param jackknife One of `"fast"`, `"loo"`, or `"none"`. `"fast"` (default)
#'   uses cached full-sample ranks with delete-1 rank updates, `"loo"`
#'   recomputes leave-one-out weights directly, and `"none"` skips SEs.
#' @param ctx Optional context passed to distance and ranking components.
#' @param return Optional extras in addition to the always-returned `mu`, `se`,
#'   and `meta`. Allowed values:
#'   `c("W", "ranks", "ord", "dx", "y", "ess", "neighbors", "tuning")`.
#'   `neighbors` is a `k x m` matrix of nearest training-observation indices,
#'   one column per target. It is a rank summary, not the full DNN weight support.
#'
#' @return An object of class `np_result` with components:
#'   - `mu`: estimated local mean values at target points;
#'   - `se`: optional jackknife standard errors;
#'   - `x0`: the `m x p` matrix of evaluation points, with covariate
#'     column names propagated from the training data when available;
#'   - `meta`: metadata (`n`, `m`, `p`, `k`, `s`, `method`, `distance`);
#'   - `tuning_grid`: in grid mode, a list with `candidates`, `estimates`,
#'     `diagnostics`, and `meta` tables. Candidate diagnostics report
#'     top-5 weight concentration separately from the estimator's `k`;
#'   - `.model`: retained canonical training state and fitted geometry used
#'     by prediction and covariance methods, independently of optional extras;
#'   - optional extras requested via `return`.
#'
#' @seealso [weights_dnn()], [np_demo_data()], [np_result] for summary,
#'   coefficient, joint covariance, confidence interval and prediction methods.
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
  s = NULL,
  tuning = NULL,
  tuning_mode = c("default", "grid"),
  normalize = TRUE,
  se = TRUE,
  jackknife = c("fast", "loo", "none"),
  ctx = NULL,
  return = c("W")
) {
  method <- match.arg(method)
  jackknife <- match.arg(jackknife)
  tuning_mode <- match.arg(tuning_mode)
  s_arg <- s
  k <- .validate_scalar_int(k, "k")
  normalize <- .validate_scalar_logical(normalize, "normalize")
  se <- .validate_scalar_logical(se, "se")
  .warn_unknown_return(return, valid = c("W", "ranks", "ord", "dx", "y", "ess", "neighbors", "tuning"))
  if (identical(tuning_mode, "grid") && !is.null(s_arg)) {
    .abort_bad_arg(
      "tuning_mode",
      paste0(
        "cannot be \"grid\" when explicit `s` is supplied; omit `s` ",
        "and use `tuning`, or omit both to use the paper default tuning grid"
      ),
      quoted = FALSE
    )
  }

  .dnn_validate_input_route(data, list(X = X, Y = Y),
    list(x_cols = x_cols, y_col = y_col))

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

  .dnn_distance_scope()
  geometry <- .dnn_fitted_geometry(distance, x_mat, ctx)
  weight_return <- return
  if (identical(tuning_mode, "grid") || (isTRUE(se) && identical(jackknife, "fast"))) {
    # Request ranks once so fast jackknife can reuse them.
    weight_return <- unique(c(return, "ranks"))
  }
  weight_fit <- .dnn_weights_with_extras(
    x_mat = x_mat,
    x0_mat = x0_mat,
    distance = geometry$distance,
    method = method,
    k = k,
    s = s,
    normalize = normalize,
    ctx = geometry$ctx,
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
      precomputed_ranks = wt$ranks, precomputed_groups = wt$groups
    )
  }

  out <- list(
    mu = mu_hat,
    se = se_hat,
    x0 = .dnn_x0_with_names(x0_mat, x_mat),
    meta = list(
      n = n, m = m_out, p = p,
      k = k, s = s,
      method = method, distance = .dnn_distance_name(distance),
      jackknife = if (isTRUE(se)) jackknife else "none",
      tie_rule = .dnn_tie_policy()
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
  if (identical(tuning_mode, "grid")) {
    out$tuning_grid <- .np_dnn_tuning_grid(
      candidates = tuning_payload$grid$candidates,
      ranks = wt$ranks, groups = wt$groups,
      mu_hat = mu_hat,
      se_hat = se_hat,
      x_mat = x_mat,
      x0_mat = x0_mat,
      y_vec = y_vec,
      distance = distance,
      method = method,
      k = k,
      normalize = normalize,
      ctx = ctx,
      se = se,
      jackknife = jackknife,
      tuning_payload = tuning_payload
    )
  }
  out$.model <- c(list(
    kind = "dnn", y_vec = y_vec, x_mat = x_mat,
    x_cols = if (!is.null(data)) .dnn_resolve_data_cols(x_cols, data, "x_cols") else NULL,
    k = k, s = s, method = method, normalize = normalize,
    ranks = wt$ranks, groups = wt$groups, dx = wt$dx
  ), geometry)

  class(out) <- "np_result"
  out
}
