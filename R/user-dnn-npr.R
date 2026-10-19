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
#' permutation at each target. The exact support extends through rank
#' `n - s + k` (`n - s + 1` in the paper-faithful `k = 1` case); this is not
#' an `s`-nearest-neighbor truncation. The `"fast"` jackknife option reuses
#' the full-sample rank matrix and updates ranks under deletion.
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
#'   The `neighbors` extra returns the first `k` nearest observations used
#'   inside each size-`s` subsample, not the full DNN weight support.
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

  geometry <- .np_fitted_geometry(distance, x_mat, ctx)
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
      precomputed_ranks = wt$ranks
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
      jackknife = if (isTRUE(se)) jackknife else "none"
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
      ranks = wt$ranks,
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
    x_cols = if (!is.null(data)) colnames(x_mat) else NULL,
    k = k, s = s, method = method, normalize = normalize,
    ranks = wt$ranks
  ), geometry)

  class(out) <- "np_result"
  out
}

.np_dnn_rank_weight_spec <- function(method, k, s, normalize) {
  .dnn_rank_weight_spec(method, k, s, normalize)
}

.np_dnn_grid_estimates <- function(mu, se) {
  out <- data.frame(
    target_id = seq_along(mu),
    component = NA_character_,
    estimate = as.numeric(mu),
    stringsAsFactors = FALSE
  )
  out$se <- if (is.null(se)) rep(NA_real_, length(mu)) else as.numeric(se)
  out
}

.np_dnn_grid_default_precomputed <- function(
  candidates, mu_hat, se_hat, n, k, method, normalize
) {
  default_row <- candidates[candidates$is_default, , drop = FALSE]
  spec <- .np_dnn_rank_weight_spec(
    method = method,
    k = k,
    s = default_row$s[[1L]],
    normalize = normalize
  )
  tab <- .rank_weight_table(spec, n)
  result <- list(
    estimates = .np_dnn_grid_estimates(mu = mu_hat, se = se_hat),
    diagnostics = .dnn_grid_diagnostics_from_table(tab, n = n, ctx = list(k = k))
  )
  stats::setNames(list(result), as.character(default_row$computed_key[[1L]]))
}

.np_dnn_grid_se_fn <- function(
  ranks, y_vec, x_mat, x0_mat, distance, method, k, normalize, ctx, jackknife
) {
  force(ranks)
  force(y_vec)
  force(x_mat)
  force(x0_mat)
  force(distance)
  force(method)
  force(k)
  force(normalize)
  force(ctx)
  force(jackknife)

  function(candidate_row, tab_jk) {
    s <- as.integer(candidate_row$s[[1L]])
    se <- if (identical(jackknife, "fast")) {
      spec <- .np_dnn_rank_weight_spec(method = method, k = k, s = s, normalize = normalize)
      loo <- if (is.null(tab_jk)) {
        .dnn_loo_loop(y_vec, ranks, spec, center = TRUE)
      } else {
        .dnn_loo_centered(
          y_vec,
          ranks,
          tab_jk,
          tab_full = attr(tab_jk, "tab_full", exact = TRUE)
        )$loo
      }
      as.numeric(.dnn_jackknife_se_from_loo(
        loo = loo,
        q = 1L,
        n = nrow(ranks),
        m = ncol(ranks)
      )[1L, ])
    } else {
      .np_dnn_jackknife_se_loo(
        x_mat = x_mat,
        x0_mat = x0_mat,
        y = y_vec,
        distance = distance,
        method = method,
        k = k,
        s = s,
        normalize = normalize,
        ctx = ctx,
        jackknife = "loo",
        precomputed_ranks = NULL
      )
    }

    data.frame(
      target_id = seq_along(se),
      component = NA_character_,
      se = as.numeric(se),
      stringsAsFactors = FALSE
    )
  }
}

.np_dnn_warn_infeasible_grid <- function(candidates) {
  skipped <- candidates[!candidates$feasible, c("candidate_id", "reason"), drop = FALSE]
  if (nrow(skipped) < 1L) {
    return(invisible(NULL))
  }
  details <- paste(
    sprintf("%s (%s)", skipped$candidate_id, skipped$reason),
    collapse = ", "
  )
  warning(
    sprintf("Skipped infeasible tuning candidate(s): %s.", details),
    call. = FALSE
  )
  invisible(NULL)
}

.np_dnn_grid_meta <- function(tuning_payload, engine_meta, n, p, k, method, distance, se, jackknife) {
  grid_meta <- tuning_payload$grid$meta
  name <- if (!is.null(grid_meta$name)) {
    as.character(grid_meta$name)
  } else if (!is.null(tuning_payload$default$meta$name)) {
    as.character(tuning_payload$default$meta$name)
  } else {
    as.character(tuning_payload$source)
  }
  list(
    name = name,
    n = as.integer(n),
    kx = as.integer(p),
    k = as.integer(k),
    method = method,
    distance = .dnn_distance_name(distance),
    jackknife = if (isTRUE(se)) jackknife else "none",
    n_candidates = engine_meta$n_candidates,
    n_feasible = engine_meta$n_feasible,
    n_skipped = engine_meta$n_skipped,
    n_failed = engine_meta$n_failed,
    n_unique = engine_meta$n_unique
  )
}

.np_dnn_tuning_grid <- function(
  candidates,
  ranks,
  mu_hat,
  se_hat,
  x_mat,
  x0_mat,
  y_vec,
  distance,
  method,
  k,
  normalize,
  ctx,
  se,
  jackknife,
  tuning_payload
) {
  .np_dnn_warn_infeasible_grid(candidates)
  se_fn <- if (isTRUE(se) && !identical(jackknife, "none")) {
    .np_dnn_grid_se_fn(
      ranks = ranks,
      y_vec = y_vec,
      x_mat = x_mat,
      x0_mat = x0_mat,
      distance = distance,
      method = method,
      k = k,
      normalize = normalize,
      ctx = ctx,
      jackknife = jackknife
    )
  } else {
    NULL
  }

  grid <- .dnn_grid_eval(
    candidates = candidates,
    ranks = ranks,
    weight_spec_fn = function(candidate_row) {
      .np_dnn_rank_weight_spec(
        method = method,
        k = k,
        s = candidate_row$s[[1L]],
        normalize = normalize
      )
    },
    estimate_fn = function(W, candidate_row) {
      mu <- as.numeric(.dnn_estimate_from_signal(y_vec, W, n = nrow(x_mat), arg = "Y")[1L, ])
      data.frame(
        target_id = seq_along(mu),
        component = NA_character_,
        estimate = mu,
        stringsAsFactors = FALSE
      )
    },
    se_fn = se_fn,
    ctx = list(k = k),
    precomputed = .np_dnn_grid_default_precomputed(
      candidates = candidates,
      mu_hat = mu_hat,
      se_hat = se_hat,
      n = nrow(x_mat),
      k = k,
      method = method,
      normalize = normalize
    )
  )
  grid$meta <- .np_dnn_grid_meta(
    tuning_payload = tuning_payload,
    engine_meta = grid$meta,
    n = nrow(x_mat),
    p = ncol(x_mat),
    k = k,
    method = method,
    distance = distance,
    se = se,
    jackknife = jackknife
  )
  grid
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
