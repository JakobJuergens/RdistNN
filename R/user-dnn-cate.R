# DNN-based CATE wrapper (rank-weighted pseudo-outcome averaging)

#' CATE via DNN rank-weighted pseudo-outcome averaging
#'
#' Estimate conditional average treatment effects by combining a pseudo-outcome
#' with DNN-based local weights.
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
#' 1. Data-frame interface (recommended): pass `data` and `x_cols`, then provide
#'    either `pseudo_outcome` directly or nuisance columns used to construct it.
#'    For AIPW nuisance construction, use `y_col`, `t_col`, `e_col`,
#'    `m0_col`, and `m1_col`.
#' 2. Matrix interface: pass `X`, one of `X0` or `x0`,
#'    and `pseudo_outcome`.
#'
#' As in [weights_dnn()], exact DNN localization uses the full rank permutation
#' at each target rather than an `s`-nearest-neighbor truncation. The `"fast"`
#' jackknife option reuses the full-sample rank matrix and updates ranks under
#' deletion. For exact DNN fits, the unaveraged rank-table support edge is rank
#' `n - s + k`, not rank `s`.
#'
#' Averaging can extend nonzero observation weights through the end of a
#' distance block that straddles the nominal rank-table edge. Effective sample
#' sizes and weight diagnostics use the resulting observation weights.
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
#' For causal identification and pointwise inference, including conditions on
#' nuisance accuracy and deletion stability, see the accompanying
#' [paper](https://jakobjuergens.com/research/#orthogonal-inference-conditional-z).
#' This function accepts supplied oracle, out-of-fold or full-sample predictions;
#' it does not verify those conditions for a learner. The normalized Robinson
#' construction is not implemented; use AIPW inputs or a justified supplied signal.
#' A pooled outcome mean cannot replace the two arm-specific means.
#'
#' With `data`, omit `X`; a direct `pseudo_outcome` vector remains supported.
#' Column selectors require `data`. Prediction retains the fitted metric and
#' context even if a registered distance name is subsequently replaced.
#'
#' @section Signal provenance:
#' `meta$signal_source` records the actual input branch: `"aipw"` for built-in
#' nuisance construction, `"supplied"` for `pseudo_outcome` (which takes
#' precedence over nuisance selectors). Neither signal values nor column names
#' establish provenance.
#' `meta$created_with_version` records the package version loaded when the
#' original fit was created. Summaries and candidate views preserve these fields.
#' Provenance does not certify causal identification, nuisance accuracy/stability,
#' localization bias, or pointwise inference.
#'
#' @param data Data frame with covariates and optional nuisance columns. Column
#'   names must be non-missing, non-empty, and unique.
#' @param x_cols Names or indices of covariate columns in `data`.
#' @param y_col Name or index of the outcome column in `data`.
#' @param t_col Name or index of the treatment column in `data` (must be 0/1).
#' @param m0_col Name or index of the supplied control-arm conditional outcome mean.
#' @param m1_col Name or index of the supplied treated-arm conditional outcome mean.
#' @param e_col Name or index of the supplied treatment propensity predictions.
#' @param clip Finite propensity clipping constant in `(0, 0.5)` with
#'   `1 - clip < 1` in double precision, default
#'   `1e-8`. Raw nuisance propensities must be finite and lie in `[0, 1]`;
#'   nuisance-based AIPW scores use
#'   `pmin(pmax(e, clip), 1 - clip)`. This argument does not alter a supplied
#'   `pseudo_outcome`.
#' @param pseudo_outcome Numeric vector of length `n`. If supplied, nuisance
#'   columns are not required and are ignored. The signal is localized unchanged;
#'   its construction and inferential justification are the caller's responsibility.
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
#'   `neighbors` is a `k x m` matrix of nearest training-observation indices,
#'   one column per target. It is a rank summary, not the full DNN weight support.
#'
#' @return An object of class `cate_result` with components:
#'   - `tau`: estimated CATE values at target points;
#'   - `se`: optional jackknife standard errors;
#'   - `x0`: the `m x p` matrix of evaluation points, with covariate
#'     column names propagated from the training data when available;
#'   - `meta`: metadata (`n`, `m`, `p`, `k`, `s`, `method`, `distance`,
#'     `signal_source`, `created_with_version`); see Signal provenance;
#'   - `tuning_grid`: in grid mode, a list with `candidates`, `estimates`,
#'     `diagnostics`, and `meta` tables. Candidate diagnostics report
#'     top-5 weight concentration separately from the estimator's `k`;
#'   - `.model`: retained observations, supplied signals and fitted geometry
#'     used by prediction and joint covariance, independently of optional extras;
#'   - optional extras requested via `return`.
#'
#' @seealso [weights_dnn()], [weights_matrix()], [score_cate_rho()],
#'   [cate_demo_data()], [vcov.cate_result()], [predict.cate_result()].
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
  e_col = NULL,
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
  clip <- .validate_propensity_clip(clip)

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
  t_vec <- e_used <- m0 <- m1 <- NULL
  signal_source <- "supplied"

  .dnn_validate_input_route(data, list(X = X),
    list(x_cols = x_cols, y_col = y_col, t_col = t_col,
      m0_col = m0_col, m1_col = m1_col, e_col = e_col))

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
      if (xor(is.null(m0_col), is.null(m1_col))) {
        .abort_bad_arg("m0_col/m1_col", "must be supplied together for full AIPW",
          quoted = FALSE)
      }
      has_aipw <- !is.null(m0_col) && !is.null(m1_col)
      if (!has_aipw) {
        .abort_bad_arg(
          "pseudo_outcome",
          "must be provided, or supply nuisance columns for the `data=` interface",
          expected = paste0(
            "`pseudo_outcome`, or full AIPW columns {`y_col`,`t_col`,`m0_col`,`m1_col`,`e_col`}"
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

      # Full doubly-robust AIPW pseudo-outcome:
      # tau_hat_i = mu1_hat - mu0_hat
      #           + T_i (Y_i - mu1_hat) / pi_hat
      #           - (1-T_i) (Y_i - mu0_hat) / (1-pi_hat)
      m0  <- .dnn_col_vec(m0_col, data, n, arg = "m0_col", ref = "data")
      m1  <- .dnn_col_vec(m1_col, data, n, arg = "m1_col", ref = "data")
      psi <- cate_rho_aipw(
        y = y_vec,
        d = t_vec,
        m1hat = m1,
        m0hat = m0,
        ehat = e_used,
        clip = clip
      )
      signal_source <- "aipw"
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

  .dnn_distance_scope()
  geometry <- .dnn_fitted_geometry(distance, x_mat, ctx)
  # Retain the resolved unfitted specification for geometry-refitting LOO.
  distance <- geometry$inference_distance

  # Weights (rank-only DNN).
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
      precomputed_ranks = wt$ranks, precomputed_groups = wt$groups
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
      jackknife = if (isTRUE(se)) jackknife else "none",
      tie_rule = .dnn_tie_policy(),
      signal_source = signal_source,
      created_with_version = as.character(getNamespaceVersion(environment(cate_dnn)))
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
    }
  }
  if (identical(tuning_mode, "grid")) {
    out$tuning_grid <- .cate_dnn_tuning_grid(
      candidates = tuning_payload$grid$candidates,
      ranks = wt$ranks, groups = wt$groups,
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
  out$.model <- c(list(
    x_mat = x_mat,
    psi = psi,
    t_vec = t_vec,
    e = e_used,
    propensity_clip = if (!is.null(e_used)) clip else NULL,
    propensity_version = if (!is.null(e_used)) "clipped" else NULL,
    x_cols = if (!is.null(data)) .dnn_resolve_data_cols(x_cols, data, "x_cols") else NULL,
    method = method,
    k = k,
    s = s,
    normalize = normalize,
    ranks = wt$ranks, groups = wt$groups, dx = wt$dx
  ), geometry)

  class(out) <- "cate_result"
  out
}
