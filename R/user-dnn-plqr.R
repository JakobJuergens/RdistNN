# DNN-based PLQR wrapper (rank-weighted local orthogonal moments)

#' Observation-level orthogonal PLQR scores at a supplied coefficient (experimental)
#'
#' Computes, for each observation \eqn{i}, the orthogonal PLQR score
#' \deqn{
#'   \psi_{i}
#'   \;=\;
#'   (\tau - \mathbf{1}\{\,Y_{i} \le D_{i}^{\top}\theta + g_{0}(U_{i})\,\})\,
#'   (D_{i} - m_{D}(U_{i})) \, ,
#' }
#' from supplied nuisance columns. No weights are applied and no localization
#' is performed; [plqr_dnn()] aggregates these scores with DNN rank weights.
#'
#' @section Lifecycle:
#' **Experimental.** The PLQR API and inference interfaces may change.
#' @param data A \code{data.frame} containing the outcome, regressor, and
#'   nuisance columns.
#' @param y_col Name or index of the outcome column in `data`.
#' @param D_cols Names or indices of partially linear regressor columns in
#'   `data`.
#' @param g_col Name or index of the shape-nuisance column in `data` holding
#'   the values of \eqn{g_0(U)}.
#' @param mD_cols Names or indices of nuisance columns for `m_D(U)`; must have
#'   the same length and order as `D_cols`.
#' @param theta Numeric coefficient vector with length `length(D_cols)`.
#' @param tau Quantile level in `(0,1)`.
#'
#' @return A numeric matrix of size \eqn{n \times p_d} of observation-level
#'   scores.
#' @seealso [plqr_dnn()], [plqr_demo_data()].
#' @examples
#' demo <- plqr_demo_data(n = 100, p = 2, p_d = 2, tau = 0.5, seed = 1)
#' scores <- plqr_score_matrix(
#'   data = demo$data,
#'   y_col = "Y",
#'   D_cols = c("D1", "D2"),
#'   g_col = "g_oracle",
#'   mD_cols = c("mD1", "mD2"),
#'   theta = demo$truth$theta_true,
#'   tau = 0.5
#' )
#' head(scores)
#' @export
plqr_score_matrix <- function(data, y_col, D_cols, g_col, mD_cols, theta, tau) {
  if (!is.data.frame(data)) {
    .abort_bad_arg("data", "must be a data.frame")
  }
  if (!is.numeric(theta) || length(theta) < 1L) {
    .abort_bad_arg("theta", "must be a non-empty numeric vector")
  }
  .abort_if_na_or_nonfinite(theta, arg = "theta")
  theta <- as.numeric(theta)
  tau <- .validate_scalar_prob(tau, arg = "tau")

  d_mat <- .dnn_col_mat(D_cols, data = data, arg = "D_cols", out_arg = "D")
  n <- nrow(d_mat)
  mD_mat <- .dnn_col_mat(mD_cols, data = data, arg = "mD_cols", out_arg = "mD")
  g_vec <- .dnn_col_vec(g_col, data, n, arg = "g_col", ref = "D_cols")
  y_vec <- .dnn_col_vec(y_col, data, n, arg = "y_col", ref = "D_cols")

  if (ncol(mD_mat) != ncol(d_mat)) {
    .abort_bad_arg(
      "mD_cols",
      "must have the same length and order as `D_cols`",
      expected = sprintf("%d columns", ncol(d_mat)),
      actual = sprintf("%d columns", ncol(mD_mat)),
      quoted = FALSE
    )
  }
  if (length(theta) != ncol(d_mat)) {
    .abort_bad_arg(
      "theta",
      "length must equal the number of `D_cols` columns",
      expected = sprintf("length %d", ncol(d_mat)),
      actual = sprintf("length(theta)=%d", length(theta))
    )
  }

  .abort_if_na_or_nonfinite(d_mat, arg = "D_cols")
  .abort_if_na_or_nonfinite(mD_mat, arg = "mD_cols")
  .abort_if_na_or_nonfinite(g_vec, arg = "g_col")
  .abort_if_na_or_nonfinite(y_vec, arg = "y_col")

  .plqr_dnn_score_matrix(d_mat, mD_mat, g_vec, y_vec, theta_eval = theta, tau = tau)
}

#' PLQR local moment evaluation and direct coefficient inference via DNN rank weights (experimental)
#'
#' Compute localized orthogonal PLQR moments
#' \deqn{
#'   \mathbb{E}\!\left[
#'     \left(\tau - \mathbf{1}\{Y \le D^\top \theta + g_0(U)\}\right)
#'     \left(D - m_D(U)\right)
#'     \;\middle|\; X \approx x_0
#'   \right]
#' }
#' at a supplied coefficient value `theta_eval`, using DNN rank weights at one
#' or multiple target points. When positive threshold density values are
#' supplied through `density` or `density_col`, the same fit reports the
#' direct PLQR coefficient sandwich used for standard errors and confidence
#' intervals. No local root solving is performed: `theta` in the result is the
#' supplied evaluation value (`meta$theta_source = "supplied"`).
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
#' inference under persistent ties. Refit saved models after package updates during experimental development.
#'
#' Two interfaces are supported.
#'
#' 1. Data-frame interface (recommended): pass `data` with `x_cols`, `D_cols`,
#'    `g_col`, `mD_cols`, and `y_col`.
#' 2. Matrix interface: pass `X`, `D`, `g`, `mD`, and `Y`.
#'
#' With `data`, omit `X`, `D`, `g`, `mD`, and `Y`; these direct inputs are
#' reserved for the matrix interface. Column selectors require `data`.
#' A direct `density` vector is supported with either interface.
#'
#' Supply the shape component `g`, residualized-regressor means `mD` and,
#' for coefficient uncertainty, positive threshold densities. These can use
#' adjustment variables beyond the reporting covariates. The threshold is
#' `D' theta_eval + g`; `g` alone is not a conditional quantile.
#'
#' Coefficient uncertainty requires `normalize = TRUE` and a density-backed
#' Jacobian. `theta_meat` chooses analytic score covariance or the full
#' jackknife covariance of the localized moment. `moment_se` alone is not a
#' coefficient standard error. For statistical background and conditions,
#' see the [paper](https://jakobjuergens.com/research/#orthogonal-inference-conditional-z).
#'
#' A Jacobian, its inverse, or the assembled sandwich covariance outside the
#' numeric range leaves coefficient uncertainty unavailable with rescaling
#' guidance. Supplied coefficients and evaluated moments remain available.
#' Finite singular Jacobians also leave coefficient uncertainty unavailable.
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
#' The payload stores one long-form row per candidate, target, and PLQR moment
#' component; the public `moment` and `moment_se` matrix fields remain anchored
#' to the default candidate. Candidate diagnostics include `ess`, `max_weight`,
#' `top_n`, `top5_weight_share`, and `degenerate`; `top_n` is fixed at 5 and
#' is independent of the estimator's `k`. Direct PLQR sandwich inference from
#' `density` or `density_col` currently uses `tuning_mode = "default"` only.
#'
#' @section Lifecycle:
#' **Experimental.** The PLQR API and inference interfaces may change.
#' @param data Optional data frame containing covariates, regressors, outcome,
#'   and nuisance columns. Column names must be non-missing, non-empty, and
#'   unique.
#' @param x_cols Names or indices of covariate columns in `data`.
#' @param D_cols Names or indices of partially linear regressor columns in
#'   `data`.
#' @param g_col Name or index of the shape-nuisance column in `data` holding
#'   the values of \eqn{g_0(U)}.
#' @param mD_cols Names or indices of nuisance columns for `m_D(U)`.
#' @param y_col Name or index of the outcome column in `data`.
#' @param theta_eval Numeric coefficient vector at which the PLQR moment is
#'   evaluated, with length `length(D_cols)` (or `ncol(D)` in the matrix
#'   interface). No root solving is performed; the result's `theta` repeats
#'   this supplied value across targets.
#' @param tau Quantile level in `(0,1)`.
#' @param X Numeric matrix (`n x p`) of covariates for the matrix interface.
#' @param D Numeric matrix (`n x p_d`) of partially linear regressors for the
#'   matrix interface.
#' @param g Numeric vector (`n`) with \eqn{g_0(U)} values for the matrix
#'   interface.
#' @param mD Numeric matrix (`n x p_d`) with `m_D(U)` values for the matrix
#'   interface.
#' @param Y Numeric outcome vector of length `n` for the matrix interface.
#' @param density Optional positive numeric vector (`n`) estimating the
#'   conditional density of `Y` at the threshold `D' theta_eval + g`. When
#'   supplied with `se = TRUE`, PLQR uses the direct coefficient sandwich
#'   route and returns coefficient standard errors in `theta_se`. A direct
#'   vector is allowed with either input interface; do not also supply
#'   `density_col`.
#' @param density_col Optional column name or index in `data` containing the
#'   positive density values used for the direct coefficient sandwich.
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
#'   `tuning` (or the paper default tuning grid when `tuning` is `NULL`) and
#'   stores candidate-level moments, standard errors, diagnostics, and metadata
#'   in `tuning_grid`. Grid mode cannot be combined with explicit `s`.
#' @param normalize Logical; renormalize columns of `W` to sum to 1.
#' @param se Logical; compute moment or density-backed coefficient uncertainty.
#' @param jackknife One of `"fast"`, `"loo"`, or `"none"`. `"fast"` (default)
#'   uses cached full-sample ranks with delete-1 rank updates, `"loo"`
#'   recomputes leave-one-out weights directly, and `"none"` skips moment
#'   jackknifing. Analytic coefficient uncertainty is controlled by `se` and
#'   the density input, independently of this selector.
#' @param theta_meat Meat estimator for the coefficient sandwich when
#'   `density`/`density_col` is supplied. `"analytic"` (default) uses the
#'   analytic score covariance
#'   \eqn{S(x) = \tau(1-\tau)\sum_i w_i \tilde D_i \tilde D_i'};
#'   `"jackknife"` uses the full delete-1 jackknife covariance of the
#'   localized vector moment (`moment_vcov`) as the meat, so
#'   \eqn{\widehat V_\theta = J^{-1} \widehat S_{JK} J^{-T}}. The bread is
#'   density-backed in both routes. `"jackknife"` requires
#'   `density`/`density_col`, `se = TRUE`, and `jackknife != "none"`.
#' @param ctx Optional context passed to distance and ranking components.
#' @param return Optional extras in addition to the always-returned `moment`,
#'   `theta`, `theta_eval`, `tau`, and `meta`. Allowed values:
#'   `c("W", "ranks", "ord", "dx", "score", "ess", "neighbors", "tuning")`.
#'   `neighbors` is a `k x m` matrix of nearest training-observation indices,
#'   one column per target. It is a rank summary, not the full DNN weight support.
#'
#' @return An object of class `rdistnn_plqr_result` with components:
#'   - `moment`: matrix (`p_d x m`) of localized PLQR moments;
#'   - `theta`: matrix (`p_d x m`) of the supplied coefficient values, repeated
#'     across targets; no root solving is performed;
#'   - `theta_eval`: the supplied coefficient vector used in the moment map;
#'   - `theta_se`: coefficient standard errors (`p_d x m`), present only when
#'     `density` or `density_col` supplies the direct PLQR sandwich route;
#'   - `moment_se`: componentwise jackknife standard errors for `moment`.
#'     Present when no density is supplied and `se = TRUE` with
#'     `jackknife != "none"` (moment-level diagnostic only), and when
#'     `theta_meat = "jackknife"` (as the diagonal of `moment_vcov`). A
#'     componentwise `moment_se` is never by itself coefficient inference;
#'   - `moment_vcov`: with `theta_meat = "jackknife"`, the full delete-1
#'     jackknife covariance of the localized vector moment, one `p_d x p_d`
#'     matrix per target — the sandwich meat;
#'   - `tau`: quantile level;
#'   - `x0`: the `m x p` matrix of evaluation points, with covariate
#'     column names propagated from the training data when available;
#'   - `meta`: metadata (`n`, `m`, `p`, `p_d`, `k`, `s`, `method`, `distance`,
#'     `theta_source = "supplied"`, `theta_se_method`, `theta_meat_method`,
#'     `moment_se_method`);
#'   - when direct coefficient inference is available, `theta_vcov`,
#'     `theta_omega`, `jacobian`, and `neff` (the paper's
#'     effective sample size \eqn{n_{\mathrm{eff}}(x) = 1/\sum_i w_i(x)^2},
#'     identical to the `ess` return extra), plus `score_covariance` (the
#'     analytic meat) when `theta_meat = "analytic"`;
#'   - `tuning_grid`: in grid mode, a list with `candidates`, `estimates`,
#'     `diagnostics`, and `meta` tables. Candidate diagnostics report
#'     top-5 weight concentration separately from the estimator's `k`;
#'   - optional extras requested via `return`.
#'
#' @seealso [weights_dnn()], [plqr_score_matrix()], [plqr_demo_data()].
#'
#' @examples
#' set.seed(1)
#' demo <- plqr_demo_data(n = 120, p = 2, p_d = 2, tau = 0.5, seed = 1)
#' out <- plqr_dnn(
#'   data = demo$data,
#'   x_cols = c("X1", "X2"),
#'   D_cols = c("D1", "D2"),
#'   g_col = "g_oracle",
#'   mD_cols = c("mD1", "mD2"),
#'   y_col = "Y",
#'   density_col = "density_oracle",
#'   theta_eval = demo$truth$theta_true,
#'   tau = 0.5,
#'   x0 = c(0, 0),
#'   method = "exact",
#'   k = 1L, s = 3L,
#'   se = TRUE
#' )
#' summary(out)
#' @export
plqr_dnn <- function(
  data = NULL,
  x_cols = NULL,
  D_cols = NULL,
  g_col = NULL,
  mD_cols = NULL,
  y_col = NULL,
  theta_eval,
  tau = 0.5,
  X = NULL,
  D = NULL,
  g = NULL,
  mD = NULL,
  Y = NULL,
  density = NULL,
  density_col = NULL,
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
  theta_meat = c("analytic", "jackknife"),
  ctx = NULL,
  return = c("W")
) {
  method <- .dnn_match_arg(method)
  jackknife <- .dnn_match_arg(jackknife)
  theta_meat <- .dnn_match_arg(theta_meat)
  tuning_mode <- .dnn_match_arg(tuning_mode)
  s_arg <- s
  k <- .validate_scalar_int(k, "k")
  normalize <- .validate_scalar_logical(normalize, "normalize")
  se <- .validate_scalar_logical(se, "se")
  .warn_unknown_return(return, valid = c("W", "ranks", "ord", "dx", "score", "ess", "neighbors", "tuning"))
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

  if (!is.numeric(theta_eval) || length(theta_eval) < 1L || any(!is.finite(theta_eval))) {
    .abort_bad_arg("theta_eval", "must be a finite numeric vector")
  }
  theta_eval <- as.numeric(theta_eval)

  .dnn_validate_input_route(data, list(X = X, D = D, g = g, mD = mD, Y = Y),
    list(x_cols = x_cols, D_cols = D_cols, g_col = g_col,
      mD_cols = mD_cols, y_col = y_col))

  if (!is.null(data)) {
    if (!is.data.frame(data)) {
      .abort_bad_arg("data", "must be a data.frame")
    }
    if (is.null(x_cols) || is.null(D_cols) || is.null(g_col) || is.null(mD_cols) || is.null(y_col)) {
      missing_args <- c("x_cols", "D_cols", "g_col", "mD_cols", "y_col")[
        c(is.null(x_cols), is.null(D_cols), is.null(g_col), is.null(mD_cols), is.null(y_col))
      ]
      .abort_bad_arg(
        "x_cols/D_cols/g_col/mD_cols/y_col",
        "must be provided with `data=` interface",
        expected = "all of `x_cols`, `D_cols`, `g_col`, `mD_cols`, and `y_col` are non-NULL",
        actual = sprintf("missing: %s", paste(missing_args, collapse = ", ")),
        quoted = FALSE
      )
    }

    x_mat <- .dnn_col_mat(x_cols, data = data, arg = "x_cols", out_arg = "X")
    n <- nrow(x_mat)
    p <- ncol(x_mat)

    d_mat <- .dnn_col_mat(D_cols, data = data, arg = "D_cols", out_arg = "D")
    mD_mat <- .dnn_col_mat(mD_cols, data = data, arg = "mD_cols", out_arg = "mD")
    g_vec <- .dnn_col_vec(g_col, data, n, arg = "g_col", ref = "X")
    y_vec <- .dnn_col_vec(y_col, data, n, arg = "y_col", ref = "X")
    if (!is.null(density_col) && !is.null(density)) {
      .abort_bad_arg(
        "density/density_col",
        "must provide at most one direct-sandwich density input",
        quoted = FALSE
      )
    }
    density_vec <- if (!is.null(density_col)) {
      .dnn_col_vec(density_col, data, n, arg = "density_col", ref = "X")
    } else if (!is.null(density)) {
      .validate_vector_nrow(density, n, arg = "density", ref = "X", strict = TRUE)
    } else {
      NULL
    }
  } else {
    if (is.null(X)) {
      .abort_bad_arg(
        "X",
        "must be provided for the matrix interface when `data` is NULL",
        expected = "non-NULL `X`, or use the `data=` interface"
      )
    }
    if (is.null(D) || is.null(g) || is.null(mD) || is.null(Y)) {
      missing_args <- c("D", "g", "mD", "Y")[c(is.null(D), is.null(g), is.null(mD), is.null(Y))]
      .abort_bad_arg(
        "D/g/mD/Y",
        "must all be provided for the matrix interface",
        expected = "non-NULL `D`, `g`, `mD`, and `Y`",
        actual = sprintf("missing: %s", paste(missing_args, collapse = ", ")),
        quoted = FALSE
      )
    }

    x_mat <- .validate_matrix_nonempty(X, arg = "X")
    n <- nrow(x_mat)
    p <- ncol(x_mat)

    d_mat <- .validate_matrix_nonempty(D, arg = "D")
    mD_mat <- .validate_matrix_nonempty(mD, arg = "mD")
    g_vec <- .validate_vector_nrow(g, n, arg = "g", ref = "X", strict = TRUE)
    y_vec <- .validate_vector_nrow(Y, n, arg = "Y", ref = "X", strict = TRUE)
    if (!is.null(density_col)) {
      .abort_bad_arg(
        "density_col",
        "can only be used with the `data=` interface",
        quoted = FALSE
      )
    }
    density_vec <- if (!is.null(density)) {
      .validate_vector_nrow(density, n, arg = "density", ref = "X", strict = TRUE)
    } else {
      NULL
    }
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
  if (length(theta_eval) != p_d) {
    .abort_bad_arg(
      "theta_eval",
      "length must equal ncol(D)",
      expected = sprintf("ncol(D)=%d", p_d),
      actual = sprintf("length(theta_eval)=%d", length(theta_eval))
    )
  }

  .abort_if_na_or_nonfinite(d_mat, arg = "D")
  .abort_if_na_or_nonfinite(mD_mat, arg = "mD")
  .abort_if_na_or_nonfinite(g_vec, arg = "g")
  .abort_if_na_or_nonfinite(y_vec, arg = "Y")
  if (!is.null(density_vec)) {
    density_vec <- .plqr_validate_density(density_vec, n, arg = "density")
  }
  if (!is.null(density_vec) && identical(tuning_mode, "grid")) {
    .abort_bad_arg(
      "tuning_mode",
      paste0(
        "cannot be \"grid\" when direct PLQR sandwich density values are ",
        "supplied; use `tuning_mode = \"default\"` or omit ",
        "`density`/`density_col`"
      ),
      quoted = FALSE
    )
  }
  if (!is.null(density_vec) && isTRUE(se) && !isTRUE(normalize)) {
    .abort_bad_arg(
      "normalize",
      paste0(
        "must be TRUE for direct PLQR coefficient inference; the sandwich ",
        "formulas assume weights that sum to one per target. Use ",
        "`normalize = TRUE`, or omit `density`/`density_col`"
      ),
      quoted = FALSE
    )
  }
  if (identical(theta_meat, "jackknife")) {
    if (is.null(density_vec)) {
      .abort_bad_arg(
        "theta_meat",
        paste0(
          "= \"jackknife\" requires threshold-density values via `density` ",
          "or `density_col`; the jackknife supplies the sandwich meat, but ",
          "the bread is always density-backed"
        ),
        quoted = FALSE
      )
    }
    if (!isTRUE(se)) {
      .abort_bad_arg(
        "theta_meat",
        "= \"jackknife\" requires `se = TRUE`",
        quoted = FALSE
      )
    }
    if (identical(jackknife, "none")) {
      .abort_bad_arg(
        "theta_meat",
        paste0(
          "= \"jackknife\" requires `jackknife = \"fast\"` or ",
          "`jackknife = \"loo\"` to compute the moment covariance"
        ),
        quoted = FALSE
      )
    }
  }

  moment_jackknife <- isTRUE(se) &&
    (is.null(density_vec) || identical(theta_meat, "jackknife"))
  s_resolved <- .dnn_resolve_s(
    s = s_arg,
    tuning = tuning,
    n = n,
    p = p,
    k = k,
    se = moment_jackknife,
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

  d_resid_mat <- d_mat - mD_mat
  below_threshold <- .plqr_below_threshold(d_mat, theta_eval, g_vec, y_vec)
  score_mat <- .plqr_dnn_score_matrix(d_mat, mD_mat, g_vec, y_vec,
    theta_eval = theta_eval, tau = tau, below_threshold = below_threshold)
  moment_hat <- .dnn_estimate_from_signal(score_mat, weight_matrix, n = n, arg = "score")
  rownames(moment_hat) <- paste0("comp_", seq_len(p_d))
  colnames(moment_hat) <- paste0("target_", seq_len(m_out))
  target_names <- colnames(moment_hat)
  coef_names <- colnames(d_mat)
  if (is.null(coef_names)) {
    coef_names <- paste0("theta_", seq_len(p_d))
  }
  theta_hat <- matrix(rep(theta_eval, times = m_out), nrow = p_d, ncol = m_out)
  theta_hat <- .plqr_dimnames_matrix(theta_hat, coef_names, target_names)

  theta_se <- NULL
  moment_se <- NULL
  direct_fit <- NULL
  direct_available <- !is.null(density_vec)
  if (isTRUE(se) && isTRUE(direct_available)) {
    if (identical(theta_meat, "jackknife")) {
      direct_fit <- .plqr_jackknife_meat_sandwich(
        x_mat = x_mat,
        x0_mat = x0_mat,
        score_mat = score_mat,
        d_mat = d_mat,
        mD_mat = mD_mat,
        density_vec = density_vec,
        W = weight_matrix,
        distance = distance,
        method = method,
        k = k,
        s = s,
        normalize = normalize,
        ctx = ctx,
        jackknife = jackknife,
        precomputed_ranks = wt$ranks, precomputed_groups = wt$groups,
        coef_names = coef_names,
        target_names = target_names
      )
      moment_se <- direct_fit$moment_se
    } else {
      direct_fit <- .plqr_direct_sandwich(
        d_mat = d_mat,
        mD_mat = mD_mat,
        density_vec = density_vec,
        W = weight_matrix,
        tau = tau,
        coef_names = coef_names,
        target_names = target_names
      )
    }
    theta_se <- direct_fit$theta_se
  } else if (isTRUE(se) && jackknife != "none") {
    moment_se <- .plqr_dnn_jackknife_se_loo(
      x_mat = x_mat,
      x0_mat = x0_mat,
      score_mat = score_mat,
      distance = distance,
      method = method,
      k = k,
      s = s,
      normalize = normalize,
      ctx = ctx,
      jackknife = jackknife,
      precomputed_ranks = wt$ranks, precomputed_groups = wt$groups
    )
    rownames(moment_se) <- rownames(moment_hat)
    colnames(moment_se) <- colnames(moment_hat)
  }

  theta_meat_method <- if (isTRUE(se) && isTRUE(direct_available)) theta_meat else "none"
  moment_se_available <- !is.null(moment_se)
  out <- list(
    moment = moment_hat,
    theta = theta_hat,
    theta_se = theta_se,
    moment_se = moment_se,
    theta_eval = theta_eval,
    tau = tau,
    x0 = .dnn_x0_with_names(x0_mat, x_mat),
    meta = list(
      n = n, m = m_out, p = p, p_d = p_d,
      k = k, s = s,
      method = method, distance = .dnn_distance_name(distance),
      jackknife = if (isTRUE(se) &&
        (!isTRUE(direct_available) || identical(theta_meat, "jackknife"))) {
        jackknife
      } else {
        "none"
      },
      tie_rule = .dnn_tie_policy(),
      theta_source = "supplied",
      theta_se_method = if (isTRUE(se) && isTRUE(direct_available)) "direct_sandwich" else "none",
      theta_meat_method = theta_meat_method,
      moment_se_method = if (moment_se_available) "jackknife" else "none",
      coef_names = coef_names
    )
  )
  if (!is.null(direct_fit)) {
    out$theta_vcov <- direct_fit$theta_vcov
    out$theta_omega <- direct_fit$theta_omega
    out$jacobian <- direct_fit$jacobian
    out$score_covariance <- direct_fit$score_covariance
    out$moment_vcov <- direct_fit$moment_vcov
    out$neff <- direct_fit$neff
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
  if ("score" %in% return) out$score <- score_mat
  if (identical(tuning_mode, "grid")) {
    out$tuning_grid <- .plqr_dnn_tuning_grid(
      candidates = tuning_payload$grid$candidates,
      ranks = wt$ranks, groups = wt$groups,
      moment_hat = moment_hat,
      se_hat = moment_se,
      x_mat = x_mat,
      x0_mat = x0_mat,
      score_mat = score_mat,
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
  # below_threshold is the indicator 1{Y <= D'theta_eval + g} entering the
  # quantile score; its weighted mean per target is the local quantile
  # calibration checked by plot_calibration().
  out$.model <- c(list(
    x_mat = x_mat,
    d_mat = d_mat,
    mD_mat = mD_mat,
    method = method,
    k = k,
    s = s,
    normalize = normalize,
    ranks = wt$ranks, groups = wt$groups, dx = wt$dx,
    d_resid_mat = d_resid_mat,
    below_threshold = as.numeric(below_threshold),
    density = density_vec,
    theta_eval = theta_eval,
    tau = tau,
    x_cols = if (!is.null(data)) .dnn_resolve_data_cols(x_cols, data, "x_cols") else NULL,
    coef_names = coef_names
  ), geometry)

  class(out) <- "rdistnn_plqr_result"
  out
}
