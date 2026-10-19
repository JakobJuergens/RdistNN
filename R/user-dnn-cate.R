# DNN-based CATE wrapper (rank-weighted local DR learner)

#' CATE via DNN rank weights (orthogonalized local averaging)
#'
#' Estimate conditional average treatment effects by combining a pseudo-outcome
#' with DNN-based local weights.
#'
#' @details
#' Two interfaces are supported.
#'
#' 1. Data-frame interface (recommended): pass `data` and `x_cols`, then provide
#'    either `pseudo_outcome` directly or nuisance columns used to construct it.
#'    For nuisance columns, use either:
#'    - full AIPW form: `y_col`, `t_col`, `e_col`, `m0_col`, `m1_col`;
#'    - Robinson form: `y_col`, `t_col`, `e_col`, `m_col`.
#' 2. Matrix interface: pass `X`, one of `X0` or `x0`,
#'    and `pseudo_outcome`.
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
#' The pseudo-outcome and nuisance columns are resolved once before grid
#' evaluation. The primary `tau`, `se`, `meta`, and optional extras such as
#' `W`, `ranks`, `ess`, `psi`, and `nuisance` remain anchored to the default
#' candidate. Full per-candidate weight matrices are not stored. Candidate
#' diagnostics include `ess`, `max_weight`, `top_n`, `top5_weight_share`, and
#' `degenerate`; `top_n` is fixed at 5 and is independent of the estimator's
#' `k`.
#'
#' Neyman orthogonality controls first-order population nuisance drift, but it
#' does not by itself validate arbitrary nuisance predictions fitted on the
#' same observations used for localization and inference. Valid full-sample
#' reuse additionally requires the nuisance-rate and full-state stability
#' conditions in the accompanying paper. This function computes the estimator
#' and its requested jackknife procedure; it cannot diagnose whether an
#' arbitrary learner satisfies those conditions. Supply oracle quantities,
#' cross-fitted or out-of-fold predictions, or full-sample predictions from a
#' learner for which the paper's rate and stability conditions are credible.
#'
#' @param data Data frame with covariates and optional nuisance columns. Column
#'   names must be non-missing, non-empty, and unique.
#' @param x_cols Names or indices of covariate columns in `data`.
#' @param y_col Name or index of the outcome column in `data`.
#' @param t_col Name or index of the treatment column in `data` (must be 0/1).
#' @param m0_col Name or index for \eqn{E[Y(0)\mid X]} (full AIPW interface).
#' @param m1_col Name or index for \eqn{E[Y(1)\mid X]} (full AIPW interface).
#' @param m_col Name or index for \eqn{E[Y\mid X]} (Robinson interface).
#' @param e_col Name or index for propensity \eqn{P(T=1\mid X)}.
#' @param clip Finite propensity clipping constant in `(0, 0.5)`, default
#'   `1e-8`. Raw nuisance propensities must be finite and lie in `[0, 1]`;
#'   nuisance-based AIPW and Robinson scores use
#'   `pmin(pmax(e, clip), 1 - clip)`. This argument does not alter a supplied
#'   `pseudo_outcome`.
#' @param pseudo_outcome Numeric vector of length `n`. If supplied, nuisance
#'   columns are not required.
#' @param X Numeric matrix (`n x p`) for the matrix interface.
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
#' @param return Optional extras in addition to the always-returned `tau`, `se`,
#'   and `meta`. Allowed values:
#'   `c("W", "ranks", "ord", "dx", "psi", "nuisance", "ess", "neighbors", "tuning")`.
#'   The `neighbors` extra returns the first `k` nearest observations used
#'   inside each size-`s` subsample, not the full DNN weight support.
#'
#' @return An object of class `cate_result` with components:
#'   - `tau`: estimated CATE values at target points;
#'   - `se`: optional jackknife standard errors;
#'   - `x0`: the `m x p` matrix of evaluation points, with covariate
#'     column names propagated from the training data when available;
#'   - `meta`: metadata (`n`, `m`, `p`, `k`, `s`, `method`, `distance`);
#'   - `tuning_grid`: in grid mode, a list with `candidates`, `estimates`,
#'     `diagnostics`, and `meta` tables. Candidate diagnostics report
#'     top-5 weight concentration separately from the estimator's `k`;
#'   - optional extras requested via `return`.
#'
#' @seealso [weights_dnn()], [weights_matrix()], [score_cate_rho()],
#'   [cate_demo_data()].
#'
#' @examples
#' set.seed(1)
#' n <- 50
#' x_mat <- matrix(rnorm(n * 2), n, 2)
#' x0_mat <- matrix(rnorm(6), 3, 2)
#' psi <- rnorm(n)
#' out <- cate_dnn(
#'   X = x_mat, X0 = x0_mat, pseudo_outcome = psi,
#'   distance = "l2", method = "exact", k = 1L, s = 3L,
#'   se = FALSE, return = c("W", "ess")
#' )
#' out$tau
#' @export
cate_dnn <- function(
  data = NULL,
  x_cols = NULL, y_col = NULL, t_col = NULL,
  m0_col = NULL, m1_col = NULL,
  m_col = NULL, e_col = NULL,
  pseudo_outcome = NULL,
  X = NULL, X0 = NULL, x0 = NULL,
  clip = 1e-8,
  distance = "l2",
  method = c("exact", "asymp"),
  k = 1L, s = NULL,
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
  if (!is.numeric(clip) || length(clip) != 1L || is.na(clip) ||
      !is.finite(clip) || clip <= 0 || clip >= 0.5) {
    .abort_bad_arg("clip", "must be a finite numeric scalar in (0, 0.5)")
  }
  clip <- as.numeric(clip)

  .warn_unknown_return(return, valid = c("W", "ranks", "ord", "dx", "psi", "nuisance", "ess", "neighbors", "tuning"))
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
  t_vec <- e_used <- m0 <- m1 <- m <- NULL

  # Dispatch: data-frame API vs matrix API.
  if (!is.null(data)) {
    if (!is.data.frame(data)) {
      .abort_bad_arg("data", "must be a data.frame")
    }
    x_mat <- .dnn_col_mat(x_cols, data = data, arg = "x_cols", out_arg = "X")
    n <- nrow(x_mat)
    p <- ncol(x_mat)

    if (!is.null(pseudo_outcome)) {
      psi <- .validate_vector_nrow(pseudo_outcome, n, arg = "pseudo_outcome", ref = "X", strict = TRUE)
    } else {
      # Determine which outcome-model interface is used
      has_aipw <- !is.null(m0_col) && !is.null(m1_col)
      has_robinson <- !is.null(m_col)
      if (has_aipw && has_robinson) {
        .abort_bad_arg(
          "m_col/m0_col/m1_col",
          "must specify either {`m0_col`, `m1_col`} (full AIPW) or `m_col` (Robinson), but not both",
          quoted = FALSE
        )
      }
      if (!has_aipw && !has_robinson) {
        .abort_bad_arg(
          "pseudo_outcome",
          "must be provided, or supply nuisance columns for the `data=` interface",
          expected = paste0(
            "`pseudo_outcome`, or full AIPW columns {`y_col`,`t_col`,`m0_col`,`m1_col`,`e_col`}, ",
            "or Robinson columns {`y_col`,`t_col`,`m_col`,`e_col`}"
          )
        )
      }
      if (is.null(y_col) || is.null(t_col) || is.null(e_col)) {
        missing_args <- c("y_col", "t_col", "e_col")[c(is.null(y_col), is.null(t_col), is.null(e_col))]
        .abort_bad_arg(
          "y_col/t_col/e_col",
          "must all be provided with `data=` when constructing pseudo-outcomes from nuisance columns",
          expected = "`y_col`, `t_col`, and `e_col` are all non-NULL",
          actual = sprintf("missing: %s", paste(missing_args, collapse = ", ")),
          quoted = FALSE
        )
      }

      y_vec <- .dnn_col_vec(y_col,  data, n, arg = "y_col",  ref = "data")
      t_vec <- .dnn_col_vec(t_col,  data, n, arg = "t_col",  ref = "data")
      t_vec <- .validate_binary01(t_vec, arg = "t_col")
      e_raw <- .dnn_col_vec(e_col, data, n, arg = "e_col", ref = "data")
      e_used <- .cate_validate_and_clip_propensity(
        e_raw, n = n, clip = clip, arg = "e_col"
      )

      if (has_aipw) {
        # Full doubly-robust AIPW pseudo-outcome:
        # tau_hat_i = mu1_hat - mu0_hat
        #           + T_i (Y_i - mu1_hat) / pi_hat
        #           - (1-T_i) (Y_i - mu0_hat) / (1-pi_hat)
        m0  <- .dnn_col_vec(m0_col, data, n, arg = "m0_col", ref = "data")
        m1  <- .dnn_col_vec(m1_col, data, n, arg = "m1_col", ref = "data")
        m   <- NULL
        psi <- cate_rho_aipw(
          y = y_vec,
          d = t_vec,
          m1hat = m1,
          m0hat = m0,
          ehat = e_used,
          clip = clip
        )
      } else {
        # Robinson-form: (Y - m(X))(T - e) / (e(1-e))
        m0 <- m1 <- NULL
        m   <- .dnn_col_vec(m_col,  data, n, arg = "m_col",  ref = "data")
        psi <- cate_rho_rlearner(
          y = y_vec,
          d = t_vec,
          ehat = e_used,
          mhat = m,
          clip = clip
        )
      }
    }
  } else {
    # matrix API
    if (is.null(X)) {
      .abort_bad_arg(
        "X",
        "must be provided for the matrix interface when `data` is NULL",
        expected = "non-NULL `X`, or use the `data=` interface"
      )
    }
    x_mat <- .validate_matrix_nonempty(X, arg = "X")
    n <- nrow(x_mat)
    p <- ncol(x_mat)
    if (is.null(pseudo_outcome)) {
      .abort_bad_arg("pseudo_outcome", "is required for the matrix interface")
    }
    psi <- .validate_vector_nrow(pseudo_outcome, n, arg = "pseudo_outcome", ref = "X", strict = TRUE)
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

  # Targets
  x0_mat <- .dnn_resolve_targets(x_mat, X0 = X0, x0 = x0)

  # Weights (rank-only DNN).
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

  # Point estimates.
  tau <- as.numeric(.dnn_estimate_from_signal(psi, weight_matrix, n = n, arg = "pseudo_outcome")[1L, ])
  se_hat <- NULL
  if (isTRUE(se) && jackknife != "none") {
    se_hat <- .cate_dnn_jackknife_se_loo(
      x_mat = x_mat,
      x0_mat = x0_mat,
      psi = psi,
      distance = distance,
      method = method,
      k = as.integer(k),
      s = as.integer(s),
      normalize = normalize,
      ctx = ctx,
      jackknife = jackknife,
      precomputed_ranks = wt$ranks
    )
  }

  # Assemble output.
  out <- list(
    tau = tau,
    se = se_hat,
    x0 = .dnn_x0_with_names(x0_mat, x_mat),
    meta = list(
      n = n, m = m_out, p = p,
      k = as.integer(k), s = as.integer(s),
      method = method, distance = .dnn_distance_name(distance),
      jackknife = if (isTRUE(se)) jackknife else "none"
    )
  )
  if (!is.null(e_used)) {
    out$meta$propensity_clip <- clip
    out$meta$propensity_version <- "clipped"
  }
  if (tuning_used) out$meta$tuning <- tuning_payload
  out <- .dnn_attach_common_extras(
    out,
    return = return,
    wt = wt,
    W = weight_matrix,
    k = k,
    tuning_payload = tuning_payload
  )
  if ("psi" %in% return) out$psi <- psi
  if ("nuisance" %in% return) {
    if (!is.null(m0) && !is.null(m1) && !is.null(e_used)) {
      out$nuisance <- list(m0 = m0, m1 = m1, e = e_used)
    } else if (!is.null(m) && !is.null(e_used)) {
      out$nuisance <- list(m = m, e = e_used)
    }
  }
  if (identical(tuning_mode, "grid")) {
    out$tuning_grid <- .cate_dnn_tuning_grid(
      candidates = tuning_payload$grid$candidates,
      ranks = wt$ranks,
      tau_hat = tau,
      se_hat = se_hat,
      x_mat = x_mat,
      x0_mat = x0_mat,
      psi = psi,
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

  # Internal fit state used by predict.cate_result().
  out$.model <- list(
    x_mat = x_mat,
    psi = psi,
    t_vec = t_vec,
    e = e_used,
    propensity_clip = if (!is.null(e_used)) clip else NULL,
    propensity_version = if (!is.null(e_used)) "clipped" else NULL,
    x_cols = if (!is.null(data)) colnames(x_mat) else NULL,
    distance = distance,
    method = method,
    k = k,
    s = s,
    normalize = normalize,
    ctx = ctx,
    ranks = if (identical(tuning_mode, "grid")) wt$ranks else NULL
  )

  class(out) <- "cate_result"
  out
}

.cate_dnn_rank_weight_spec <- function(method, k, s, normalize) {
  .dnn_rank_weight_spec(method, k, s, normalize)
}

.cate_dnn_grid_estimates <- function(tau, se) {
  out <- data.frame(
    target_id = seq_along(tau),
    component = NA_character_,
    estimate = as.numeric(tau),
    stringsAsFactors = FALSE
  )
  out$se <- if (is.null(se)) rep(NA_real_, length(tau)) else as.numeric(se)
  out
}

.cate_dnn_grid_default_precomputed <- function(
  candidates, tau_hat, se_hat, n, k, method, normalize
) {
  default_row <- candidates[candidates$is_default, , drop = FALSE]
  spec <- .cate_dnn_rank_weight_spec(
    method = method,
    k = k,
    s = default_row$s[[1L]],
    normalize = normalize
  )
  tab <- .rank_weight_table(spec, n)
  result <- list(
    estimates = .cate_dnn_grid_estimates(tau = tau_hat, se = se_hat),
    diagnostics = .dnn_grid_diagnostics_from_table(tab, n = n, ctx = list(k = k))
  )
  stats::setNames(list(result), as.character(default_row$computed_key[[1L]]))
}

.cate_dnn_grid_se_fn <- function(
  ranks, psi, x_mat, x0_mat, distance, method, k, normalize, ctx, jackknife
) {
  force(ranks)
  force(psi)
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
      spec <- .cate_dnn_rank_weight_spec(method = method, k = k, s = s, normalize = normalize)
      loo <- if (is.null(tab_jk)) {
        .dnn_loo_loop(psi, ranks, spec, center = TRUE)
      } else {
        .dnn_loo_centered(
          psi,
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
      .cate_dnn_jackknife_se_loo(
        x_mat = x_mat,
        x0_mat = x0_mat,
        psi = psi,
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

.cate_dnn_warn_infeasible_grid <- function(candidates) {
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

.cate_dnn_grid_meta <- function(tuning_payload, engine_meta, n, p, k, method, distance, se, jackknife) {
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

.cate_dnn_tuning_grid <- function(
  candidates,
  ranks,
  tau_hat,
  se_hat,
  x_mat,
  x0_mat,
  psi,
  distance,
  method,
  k,
  normalize,
  ctx,
  se,
  jackknife,
  tuning_payload
) {
  .cate_dnn_warn_infeasible_grid(candidates)
  se_fn <- if (isTRUE(se) && !identical(jackknife, "none")) {
    .cate_dnn_grid_se_fn(
      ranks = ranks,
      psi = psi,
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
      .cate_dnn_rank_weight_spec(
        method = method,
        k = k,
        s = candidate_row$s[[1L]],
        normalize = normalize
      )
    },
    estimate_fn = function(W, candidate_row) {
      tau <- as.numeric(.dnn_estimate_from_signal(psi, W, n = nrow(x_mat), arg = "pseudo_outcome")[1L, ])
      data.frame(
        target_id = seq_along(tau),
        component = NA_character_,
        estimate = tau,
        stringsAsFactors = FALSE
      )
    },
    se_fn = se_fn,
    ctx = list(k = k),
    precomputed = .cate_dnn_grid_default_precomputed(
      candidates = candidates,
      tau_hat = tau_hat,
      se_hat = se_hat,
      n = nrow(x_mat),
      k = k,
      method = method,
      normalize = normalize
    )
  )
  grid$meta <- .cate_dnn_grid_meta(
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

# Internal: paper-consistent delete-1 jackknife for cate_dnn().
# Implements the delete-1 recomputation in Eq. (3.30)-(3.32) of
# Local_Z_Est_JJuergens.pdf for the scalar (p_alpha = 1) CATE target.
.cate_dnn_jackknife_se_loo <- function(
  x_mat, x0_mat, psi, distance, method, k, s, normalize, ctx,
  jackknife = c("fast", "loo"),
  precomputed_ranks = NULL
) {
  jackknife <- match.arg(jackknife)
  as.numeric(
    .dnn_jackknife_se_loo_matrix(
      x_mat = x_mat,
      x0_mat = x0_mat,
      signal = psi,
      distance = distance,
      method = method,
      k = k,
      s = s,
      normalize = normalize,
      ctx = ctx,
      signal_arg = "pseudo_outcome",
      jackknife = jackknife,
      precomputed_ranks = precomputed_ranks
    )[1L, ]
  )
}

#' Print a `cate_result` object
#'
#' @param x An object returned by [cate_dnn()].
#' @param ... Not used.
#' @return `x`, invisibly.
#' @examples
#' set.seed(1)
#' n <- 20
#' x_mat <- matrix(rnorm(n * 2), n, 2)
#' out <- cate_dnn(
#'   X = x_mat,
#'   x0 = c(0, 0),
#'   pseudo_outcome = rnorm(n),
#'   se = FALSE
#' )
#' print(out)
#' @export
print.cate_result <- function(x, ...) {
  cat("<cate_result>\n", sep = "")
  cat("  data  : n=", x$meta$n, ", p=", x$meta$p, ", m=", x$meta$m, "\n", sep = "")
  cat("  DNN   : method=", x$meta$method, ", k=", x$meta$k, ", s=", x$meta$s,
    ", distance=", x$meta$distance, "\n",
    sep = ""
  )
  if (!is.null(x$se)) {
    cat("  est   : tau[1:min(5,m)] = ",
      paste0(round(utils::head(x$tau, 5), 4), collapse = ", "),
      "  (SEs available)\n",
      sep = ""
    )
  } else {
    cat("  est   : tau[1:min(5,m)] = ",
      paste0(round(utils::head(x$tau, 5), 4), collapse = ", "),
      "  (SEs not requested)\n",
      sep = ""
    )
  }
  invisible(x)
}
