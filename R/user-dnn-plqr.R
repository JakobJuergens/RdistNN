# DNN-based PLQR wrapper (rank-weighted local orthogonal moments)

.plqr_dnn_g_matrix <- function(d_mat, mD_mat, q_vec, y_vec, beta, tau) {
  thresh <- as.vector(d_mat %*% beta) + q_vec
  resid <- as.numeric(tau) - as.numeric(y_vec <= thresh)
  (d_mat - mD_mat) * resid
}

.plqr_dnn_jackknife_se_loo <- function(
  x_mat, x0_mat, g_mat, distance, method, k, s, normalize, ctx,
  jackknife = c("fast", "loo"),
  precomputed_ranks = NULL
) {
  jackknife <- match.arg(jackknife)
  .dnn_jackknife_se_loo_matrix(
    x_mat = x_mat,
    x0_mat = x0_mat,
    signal = g_mat,
    distance = distance,
    method = method,
    k = k,
    s = s,
    normalize = normalize,
    ctx = ctx,
    signal_arg = "g",
    jackknife = jackknife,
    precomputed_ranks = precomputed_ranks
  )
}

#' PLQR local orthogonal moments via DNN rank weights
#'
#' Compute localized orthogonal PLQR moments
#' \deqn{
#'   \mathbb{E}\!\left[
#'     \left(\tau - \mathbf{1}\{Y \le D^\top \beta + q_{\tau}(X)\}\right)
#'     \left(D - m_D(X)\right)
#'     \;\middle|\; X \approx x_0
#'   \right]
#' }
#' using DNN rank weights at one or multiple target points.
#'
#' @details
#' Two interfaces are supported.
#'
#' 1. Data-frame interface (recommended): pass `data` with `x_cols`, `D_cols`,
#'    `q_col`, `mD_cols`, and `y_col`.
#' 2. Matrix interface: pass `X`, `D`, `q`, `mD`, and `Y`.
#'
#' No nuisance models are fit internally; this wrapper assumes `q_tau(X)` and
#' `m_D(X)` are supplied.
#'
#' As in [weights_dnn()], exact DNN localization uses the full rank permutation
#' at each target rather than an `s`-nearest-neighbor truncation. The `"fast"`
#' jackknife option reuses the full-sample rank matrix and updates ranks under
#' deletion. For exact DNN fits, the finite-population support edge is rank
#' `n - s + k`, not rank `s`.
#'
#' `tuning_mode = "default"` evaluates only the default kernel order and
#' preserves the ordinary scalar-fit return shape. `tuning_mode = "grid"`
#' evaluates the feasible candidates in `tuning`, or the paper default tuning
#' grid when `tuning` is omitted, and attaches a tidy `tuning_grid` payload.
#' The payload stores one long-form row per candidate, target, and PLQR moment
#' component; the public `moment` and `se` matrix fields remain anchored to the
#' default candidate. Candidate diagnostics include `ess`, `max_weight`,
#' `top_n`, `top5_weight_share`, and `degenerate`; `top_n` is fixed at 5 and
#' is independent of the estimator's `k`.
#'
#' @param data Optional data frame containing covariates, regressors, outcome,
#'   and nuisance columns. Column names must be non-missing, non-empty, and
#'   unique.
#' @param x_cols Names or indices of covariate columns in `data`.
#' @param D_cols Names or indices of partially linear regressor columns in
#'   `data`.
#' @param q_col Name or index of the `q_tau(X)` nuisance column in `data`.
#' @param mD_cols Names or indices of nuisance columns for `m_D(X)`.
#' @param y_col Name or index of the outcome column in `data`.
#' @param beta Numeric coefficient vector with length `length(D_cols)` (or
#'   `ncol(D)` in matrix interface).
#' @param tau Quantile level in `(0,1)`.
#' @param X Numeric matrix (`n x p`) of covariates for the matrix interface.
#' @param D Numeric matrix (`n x p_d`) of partially linear regressors for the
#'   matrix interface.
#' @param q Numeric vector (`n`) with `q_tau(X)` values for the matrix
#'   interface.
#' @param mD Numeric matrix (`n x p_d`) with `m_D(X)` values for the matrix
#'   interface.
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
#'   \eqn{s = \lfloor n^{(k_x+1)/(k_x+2)} \rfloor}. Note that small fixed
#'   values such as `s = 1` (or `k = s`) yield uniform weights and hence no
#'   localization; a warning is emitted in that case.
#' @param tuning Optional kernel-order tuning specification. When supplied,
#'   `s` must be omitted and is resolved from `tuning`.
#' @param tuning_mode One of `"default"` or `"grid"`. `"default"` evaluates
#'   the default `s` only. `"grid"` evaluates the feasible candidates from
#'   `tuning` (or the paper default tuning grid when `tuning` is `NULL`) and
#'   stores candidate-level moments, standard errors, diagnostics, and metadata
#'   in `tuning_grid`. Grid mode cannot be combined with explicit `s`.
#' @param normalize Logical; renormalize columns of `W` to sum to 1.
#' @param se Logical; compute delete-1 jackknife standard errors.
#' @param jackknife One of `"fast"`, `"loo"`, or `"none"`. `"fast"` (default)
#'   uses cached full-sample ranks with delete-1 rank updates, `"loo"`
#'   recomputes leave-one-out weights directly, and `"none"` skips SEs.
#' @param ctx Optional context passed to distance and ranking components.
#' @param return Optional extras in addition to the always-returned `moment`,
#'   `se`, `beta`, `tau`, and `meta`. Allowed values:
#'   `c("W", "ranks", "ord", "dx", "g", "ess", "neighbors", "tuning")`.
#'   The `neighbors` extra returns the first `k` nearest observations used
#'   inside each size-`s` subsample, not the full DNN weight support.
#'
#' @return An object of class `plqr_result` with components:
#'   - `moment`: matrix (`p_d x m`) of localized PLQR moments;
#'   - `se`: optional jackknife standard errors (`p_d x m`);
#'   - `beta`: the coefficient vector used in the moment map;
#'   - `tau`: quantile level;
#'   - `x0`: the `m x p` matrix of evaluation points, with covariate
#'     column names propagated from the training data when available;
#'   - `meta`: metadata (`n`, `m`, `p`, `p_d`, `k`, `s`, `method`, `distance`);
#'   - `tuning_grid`: in grid mode, a list with `candidates`, `estimates`,
#'     `diagnostics`, and `meta` tables. Candidate diagnostics report
#'     top-5 weight concentration separately from the estimator's `k`;
#'   - optional extras requested via `return`.
#'
#' @seealso [weights_dnn()], [score_plqr_component()], [plqr_demo_data()].
#'
#' @examples
#' set.seed(1)
#' demo <- plqr_demo_data(n = 120, p = 2, p_d = 2, tau = 0.5, seed = 1)
#' out <- plqr_dnn(
#'   data = demo$data,
#'   x_cols = c("X1", "X2"),
#'   D_cols = c("D1", "D2"),
#'   q_col = "q_oracle",
#'   mD_cols = c("mD1", "mD2"),
#'   y_col = "Y",
#'   beta = demo$truth$beta_true,
#'   tau = 0.5,
#'   x0 = c(0, 0),
#'   method = "exact",
#'   k = 1L, s = 3L,
#'   se = FALSE
#' )
#' out$moment
#' @export
plqr_dnn <- function(
  data = NULL,
  x_cols = NULL,
  D_cols = NULL,
  q_col = NULL,
  mD_cols = NULL,
  y_col = "Y",
  beta,
  tau = 0.5,
  X = NULL,
  D = NULL,
  q = NULL,
  mD = NULL,
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
  .warn_unknown_return(return, valid = c("W", "ranks", "ord", "dx", "g", "ess", "neighbors", "tuning"))
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
  tau <- .validate_scalar_prob(tau, arg = "tau")

  if (!is.numeric(beta) || length(beta) < 1L || any(!is.finite(beta))) {
    .abort_bad_arg("beta", "must be a finite numeric vector")
  }
  beta <- as.numeric(beta)

  if (!is.null(data)) {
    if (!is.data.frame(data)) {
      .abort_bad_arg("data", "must be a data.frame")
    }
    if (is.null(x_cols) || is.null(D_cols) || is.null(q_col) || is.null(mD_cols) || is.null(y_col)) {
      missing_args <- c("x_cols", "D_cols", "q_col", "mD_cols", "y_col")[
        c(is.null(x_cols), is.null(D_cols), is.null(q_col), is.null(mD_cols), is.null(y_col))
      ]
      .abort_bad_arg(
        "x_cols/D_cols/q_col/mD_cols/y_col",
        "must be provided with `data=` interface",
        expected = "all of `x_cols`, `D_cols`, `q_col`, `mD_cols`, and `y_col` are non-NULL",
        actual = sprintf("missing: %s", paste(missing_args, collapse = ", ")),
        quoted = FALSE
      )
    }

    x_mat <- .dnn_col_mat(x_cols, data = data, arg = "x_cols", out_arg = "X")
    n <- nrow(x_mat)
    p <- ncol(x_mat)

    d_mat <- .dnn_col_mat(D_cols, data = data, arg = "D_cols", out_arg = "D")
    mD_mat <- .dnn_col_mat(mD_cols, data = data, arg = "mD_cols", out_arg = "mD")
    q_vec <- .dnn_col_vec(q_col, data, n, arg = "q_col", ref = "X")
    y_vec <- .dnn_col_vec(y_col, data, n, arg = "y_col", ref = "X")
  } else {
    if (is.null(X)) {
      .abort_bad_arg(
        "X",
        "must be provided for the matrix interface when `data` is NULL",
        expected = "non-NULL `X`, or use the `data=` interface"
      )
    }
    if (is.null(D) || is.null(q) || is.null(mD) || is.null(Y)) {
      missing_args <- c("D", "q", "mD", "Y")[c(is.null(D), is.null(q), is.null(mD), is.null(Y))]
      .abort_bad_arg(
        "D/q/mD/Y",
        "must all be provided for the matrix interface",
        expected = "non-NULL `D`, `q`, `mD`, and `Y`",
        actual = sprintf("missing: %s", paste(missing_args, collapse = ", ")),
        quoted = FALSE
      )
    }

    x_mat <- .validate_matrix_nonempty(X, arg = "X")
    n <- nrow(x_mat)
    p <- ncol(x_mat)

    d_mat <- .validate_matrix_nonempty(D, arg = "D")
    mD_mat <- .validate_matrix_nonempty(mD, arg = "mD")
    q_vec <- .validate_vector_nrow(q, n, arg = "q", ref = "X", strict = TRUE)
    y_vec <- .validate_vector_nrow(Y, n, arg = "Y", ref = "X", strict = TRUE)
  }

  if (nrow(d_mat) != n || nrow(mD_mat) != n) {
    .abort_bad_arg(
      "D/mD",
      "must have row counts equal to nrow(X)",
      expected = sprintf("nrow(X)=%d", n),
      actual = sprintf("nrow(D)=%d, nrow(mD)=%d", nrow(d_mat), nrow(mD_mat)),
      quoted = FALSE
    )
  }
  if (ncol(d_mat) != ncol(mD_mat)) {
    .abort_bad_arg(
      "ncol(mD)",
      "must match ncol(D)",
      expected = sprintf("ncol(D)=%d", ncol(d_mat)),
      actual = sprintf("ncol(mD)=%d", ncol(mD_mat)),
      quoted = FALSE
    )
  }
  p_d <- ncol(d_mat)
  if (length(beta) != p_d) {
    .abort_bad_arg(
      "beta",
      "length must equal ncol(D)",
      expected = sprintf("ncol(D)=%d", p_d),
      actual = sprintf("length(beta)=%d", length(beta))
    )
  }

  .abort_if_na_or_nonfinite(d_mat, arg = "D")
  .abort_if_na_or_nonfinite(mD_mat, arg = "m_D(X)")
  .abort_if_na_or_nonfinite(q_vec, arg = "q_tau(X)")
  .abort_if_na_or_nonfinite(y_vec, arg = "Y")

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
  if (identical(tuning_mode, "grid") || (isTRUE(se) && identical(jackknife, "fast"))) {
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

  d_resid_mat <- d_mat - mD_mat
  g_mat <- .plqr_dnn_g_matrix(d_mat, mD_mat, q_vec, y_vec, beta = beta, tau = tau)
  moment_hat <- .dnn_estimate_from_signal(g_mat, weight_matrix, n = n, arg = "g")
  rownames(moment_hat) <- paste0("comp_", seq_len(p_d))
  colnames(moment_hat) <- paste0("target_", seq_len(m_out))

  se_hat <- NULL
  if (isTRUE(se) && jackknife != "none") {
    se_hat <- .plqr_dnn_jackknife_se_loo(
      x_mat = x_mat,
      x0_mat = x0_mat,
      g_mat = g_mat,
      distance = distance,
      method = method,
      k = k,
      s = s,
      normalize = normalize,
      ctx = ctx,
      jackknife = jackknife,
      precomputed_ranks = wt$ranks
    )
    rownames(se_hat) <- rownames(moment_hat)
    colnames(se_hat) <- colnames(moment_hat)
  }

  out <- list(
    moment = moment_hat,
    se = se_hat,
    beta = beta,
    tau = tau,
    x0 = .dnn_x0_with_names(x0_mat, x_mat),
    meta = list(
      n = n, m = m_out, p = p, p_d = p_d,
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
  if ("g" %in% return) out$g <- g_mat
  if (identical(tuning_mode, "grid")) {
    out$tuning_grid <- .plqr_dnn_tuning_grid(
      candidates = tuning_payload$grid$candidates,
      ranks = wt$ranks,
      moment_hat = moment_hat,
      se_hat = se_hat,
      x_mat = x_mat,
      x0_mat = x0_mat,
      g_mat = g_mat,
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
  # Internal fit state used by diagnostics and calibration plotting.
  # below_threshold is the indicator 1{Y <= D'beta + q} entering the
  # quantile score; its weighted mean per target is the local quantile
  # calibration checked by plot_calibration().
  out$.model <- list(
    x_mat = x_mat,
    distance = distance,
    method = method,
    k = k,
    s = s,
    normalize = normalize,
    ctx = ctx,
    ranks = if (identical(tuning_mode, "grid")) wt$ranks else NULL,
    d_resid_mat = d_resid_mat,
    below_threshold = as.numeric(y_vec <= as.vector(d_mat %*% beta) + q_vec)
  )

  class(out) <- "plqr_result"
  out
}

.plqr_dnn_rank_weight_spec <- function(method, k, s, normalize) {
  if (identical(method, "exact")) {
    rank_weight_dnn_exact(k = as.integer(k), s = as.integer(s))
  } else {
    rank_weight_dnn_asymp(k = as.integer(k), s = as.integer(s), normalize = normalize)
  }
}

.plqr_dnn_grid_estimates <- function(moment, se) {
  component <- rownames(moment)
  if (is.null(component)) {
    component <- paste0("comp_", seq_len(nrow(moment)))
  }
  out <- data.frame(
    target_id = rep(seq_len(ncol(moment)), each = nrow(moment)),
    component = rep(component, times = ncol(moment)),
    estimate = as.numeric(moment),
    stringsAsFactors = FALSE
  )
  out$se <- if (is.null(se)) rep(NA_real_, nrow(out)) else as.numeric(se)
  out
}

.plqr_dnn_grid_default_precomputed <- function(
  candidates, moment_hat, se_hat, n, k, method, normalize
) {
  default_row <- candidates[candidates$is_default, , drop = FALSE]
  spec <- .plqr_dnn_rank_weight_spec(
    method = method,
    k = k,
    s = default_row$s[[1L]],
    normalize = normalize
  )
  tab <- .rank_weight_table(spec, n)
  result <- list(
    estimates = .plqr_dnn_grid_estimates(moment = moment_hat, se = se_hat),
    diagnostics = .dnn_grid_diagnostics_from_table(tab, n = n, ctx = list(k = k))
  )
  stats::setNames(list(result), as.character(default_row$computed_key[[1L]]))
}

.plqr_dnn_grid_se_fn <- function(
  ranks, g_mat, x_mat, x0_mat, distance, method, k, normalize, ctx, jackknife
) {
  force(ranks)
  force(g_mat)
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
      spec <- .plqr_dnn_rank_weight_spec(method = method, k = k, s = s, normalize = normalize)
      loo <- if (is.null(tab_jk)) {
        .dnn_loo_loop(g_mat, ranks, spec)
      } else {
        .dnn_loo_prefix(
          g_mat,
          ranks,
          tab_jk,
          tab_full = attr(tab_jk, "tab_full", exact = TRUE)
        )
      }
      .dnn_jackknife_se_from_loo(
        loo = loo,
        q = ncol(g_mat),
        n = nrow(ranks),
        m = ncol(ranks)
      )
    } else {
      .plqr_dnn_jackknife_se_loo(
        x_mat = x_mat,
        x0_mat = x0_mat,
        g_mat = g_mat,
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

    out <- .plqr_dnn_grid_estimates(moment = se, se = NULL)
    data.frame(
      target_id = out$target_id,
      component = out$component,
      se = out$estimate,
      stringsAsFactors = FALSE
    )
  }
}

.plqr_dnn_warn_infeasible_grid <- function(candidates) {
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

.plqr_dnn_grid_meta <- function(
  tuning_payload, engine_meta, n, p, p_d, k, method, distance, se, jackknife
) {
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
    p_d = as.integer(p_d),
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

.plqr_dnn_tuning_grid <- function(
  candidates,
  ranks,
  moment_hat,
  se_hat,
  x_mat,
  x0_mat,
  g_mat,
  distance,
  method,
  k,
  normalize,
  ctx,
  se,
  jackknife,
  tuning_payload
) {
  .plqr_dnn_warn_infeasible_grid(candidates)
  se_fn <- if (isTRUE(se) && !identical(jackknife, "none")) {
    .plqr_dnn_grid_se_fn(
      ranks = ranks,
      g_mat = g_mat,
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
      .plqr_dnn_rank_weight_spec(
        method = method,
        k = k,
        s = candidate_row$s[[1L]],
        normalize = normalize
      )
    },
    estimate_fn = function(W, candidate_row) {
      moment <- .dnn_estimate_from_signal(g_mat, W, n = nrow(x_mat), arg = "g")
      rownames(moment) <- paste0("comp_", seq_len(nrow(moment)))
      colnames(moment) <- paste0("target_", seq_len(ncol(moment)))
      .plqr_dnn_grid_estimates(moment = moment, se = NULL)
    },
    se_fn = se_fn,
    ctx = list(k = k),
    precomputed = .plqr_dnn_grid_default_precomputed(
      candidates = candidates,
      moment_hat = moment_hat,
      se_hat = se_hat,
      n = nrow(x_mat),
      k = k,
      method = method,
      normalize = normalize
    )
  )
  grid$meta <- .plqr_dnn_grid_meta(
    tuning_payload = tuning_payload,
    engine_meta = grid$meta,
    n = nrow(x_mat),
    p = ncol(x_mat),
    p_d = ncol(g_mat),
    k = k,
    method = method,
    distance = distance,
    se = se,
    jackknife = jackknife
  )
  grid
}
