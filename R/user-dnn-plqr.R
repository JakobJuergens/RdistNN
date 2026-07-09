# DNN-based PLQR wrapper (rank-weighted local orthogonal moments)

.plqr_dnn_score_matrix <- function(d_mat, mD_mat, g_vec, y_vec, theta_eval, tau) {
  thresh <- as.vector(d_mat %*% theta_eval) + g_vec
  resid <- as.numeric(tau) - as.numeric(y_vec <= thresh)
  (d_mat - mD_mat) * resid
}

.plqr_dnn_jackknife_se_loo <- function(
  x_mat, x0_mat, score_mat, distance, method, k, s, normalize, ctx,
  jackknife = c("fast", "loo"),
  precomputed_ranks = NULL
) {
  jackknife <- match.arg(jackknife)
  .dnn_jackknife_se_loo_matrix(
    x_mat = x_mat,
    x0_mat = x0_mat,
    signal = score_mat,
    distance = distance,
    method = method,
    k = k,
    s = s,
    normalize = normalize,
    ctx = ctx,
    signal_arg = "score",
    jackknife = jackknife,
    precomputed_ranks = precomputed_ranks
  )
}

#' Observation-level orthogonal PLQR scores at a supplied coefficient
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

.plqr_validate_density <- function(density, n, arg = "density") {
  density_vec <- .validate_vector_nrow(
    density,
    n,
    arg = arg,
    ref = "X",
    strict = TRUE,
    finite = TRUE
  )
  if (any(density_vec <= 0)) {
    .abort_bad_arg(arg, "must contain positive density values", quoted = FALSE)
  }
  density_vec
}

.plqr_safe_solve <- function(A) {
  tryCatch(solve(A), error = function(e) NULL)
}

.plqr_dimnames_matrix <- function(x, coef_names, target_names) {
  rownames(x) <- coef_names
  colnames(x) <- target_names
  x
}

# The direct-sandwich theorem convention assumes weights that sum to one per
# target; unnormalized weights would make coefficient SEs scale-dependent.
.plqr_validate_normalized_weights <- function(W, tol = 1e-8) {
  col_sums <- colSums(as.matrix(W))
  if (any(!is.finite(col_sums)) || any(abs(col_sums - 1) > tol)) {
    .abort_bad_state(
      "Direct PLQR sandwich requires normalized weight columns (each summing to 1)."
    )
  }
  invisible(NULL)
}

# Density-backed bread of the PLQR coefficient sandwich: per-target Jacobian
# J(x) = sum_i w_i f_i tilde_D_i D_i' plus effective sample sizes.
.plqr_direct_bread <- function(d_mat, mD_mat, density_vec, W,
                               coef_names, target_names) {
  d_resid_mat <- d_mat - mD_mat
  m <- ncol(W)
  jacobian <- vector("list", m)
  for (j in seq_len(m)) {
    J <- crossprod(d_resid_mat, d_mat * as.numeric(density_vec * W[, j]))
    dimnames(J) <- list(coef_names, coef_names)
    jacobian[[j]] <- J
  }
  names(jacobian) <- target_names
  list(
    jacobian = jacobian,
    neff = stats::setNames(as.numeric(1 / colSums(W^2)), target_names)
  )
}

# Sandwich J^-1 meat J^-T per target. `meat_is_vcov = FALSE` treats the meat
# as an asymptotic score covariance (analytic route): omega = J^-1 S J^-T and
# vcov = omega / neff. `meat_is_vcov = TRUE` treats the meat as a direct
# estimate of Var(moment) (jackknife route): vcov = J^-1 meat J^-T and
# omega = vcov * neff.
.plqr_sandwich_from_meat <- function(jacobian, meat, neff,
                                     coef_names, target_names,
                                     meat_is_vcov) {
  m <- length(jacobian)
  p_d <- length(coef_names)
  theta_omega <- vector("list", m)
  theta_vcov <- vector("list", m)
  theta_se <- matrix(NA_real_, nrow = p_d, ncol = m)

  for (j in seq_len(m)) {
    J_inv <- .plqr_safe_solve(jacobian[[j]])
    if (is.null(J_inv) || any(!is.finite(J_inv))) {
      warning(sprintf(
        "plqr_dnn: direct PLQR Jacobian is singular at target %d; coefficient SEs set to NA.",
        j
      ), call. = FALSE)
      omega_j <- matrix(NA_real_, nrow = p_d, ncol = p_d)
      vcov_j <- matrix(NA_real_, nrow = p_d, ncol = p_d)
    } else {
      sandwiched <- J_inv %*% meat[[j]] %*% t(J_inv)
      sandwiched <- (sandwiched + t(sandwiched)) / 2
      if (isTRUE(meat_is_vcov)) {
        vcov_j <- sandwiched
        omega_j <- sandwiched * neff[[j]]
      } else {
        omega_j <- sandwiched
        vcov_j <- sandwiched / neff[[j]]
      }
      theta_se[, j] <- sqrt(pmax(diag(vcov_j), 0))
    }
    dimnames(omega_j) <- list(coef_names, coef_names)
    dimnames(vcov_j) <- list(coef_names, coef_names)
    theta_omega[[j]] <- omega_j
    theta_vcov[[j]] <- vcov_j
  }

  names(theta_omega) <- names(theta_vcov) <- target_names
  list(
    theta_se = .plqr_dimnames_matrix(theta_se, coef_names, target_names),
    theta_vcov = theta_vcov,
    theta_omega = theta_omega
  )
}

.plqr_direct_sandwich <- function(d_mat, mD_mat, density_vec, W, tau,
                                  coef_names, target_names) {
  .plqr_validate_normalized_weights(W)
  d_resid_mat <- d_mat - mD_mat
  m <- ncol(W)

  bread <- .plqr_direct_bread(
    d_mat = d_mat,
    mD_mat = mD_mat,
    density_vec = density_vec,
    W = W,
    coef_names = coef_names,
    target_names = target_names
  )

  score_covariance <- vector("list", m)
  for (j in seq_len(m)) {
    S <- as.numeric(tau) * (1 - as.numeric(tau)) *
      crossprod(d_resid_mat, d_resid_mat * as.numeric(W[, j]))
    S <- (S + t(S)) / 2
    dimnames(S) <- list(coef_names, coef_names)
    score_covariance[[j]] <- S
  }
  names(score_covariance) <- target_names

  sandwich <- .plqr_sandwich_from_meat(
    jacobian = bread$jacobian,
    meat = score_covariance,
    neff = bread$neff,
    coef_names = coef_names,
    target_names = target_names,
    meat_is_vcov = FALSE
  )

  c(
    sandwich,
    list(
      jacobian = bread$jacobian,
      score_covariance = score_covariance,
      neff = bread$neff
    )
  )
}

# Jackknife-meat coefficient sandwich: full delete-1 covariance of the
# localized vector moment as the meat, density-backed bread unchanged.
.plqr_jackknife_meat_sandwich <- function(
  x_mat, x0_mat, score_mat, d_mat, mD_mat, density_vec, W,
  distance, method, k, s, normalize, ctx, jackknife,
  precomputed_ranks, coef_names, target_names
) {
  .plqr_validate_normalized_weights(W)

  loo <- .dnn_jackknife_loo_array(
    x_mat = x_mat,
    x0_mat = x0_mat,
    signal = score_mat,
    distance = distance,
    method = method,
    k = k,
    s = s,
    normalize = normalize,
    ctx = ctx,
    signal_arg = "score",
    jackknife = jackknife,
    precomputed_ranks = precomputed_ranks
  )
  jk <- .dnn_jackknife_vcov_from_loo(
    loo = loo,
    q = ncol(score_mat),
    n = nrow(x_mat),
    m = ncol(W)
  )
  moment_vcov <- lapply(jk$vcov, function(v) {
    dimnames(v) <- list(coef_names, coef_names)
    v
  })
  names(moment_vcov) <- target_names
  moment_se <- .plqr_dimnames_matrix(jk$se, coef_names, target_names)

  bread <- .plqr_direct_bread(
    d_mat = d_mat,
    mD_mat = mD_mat,
    density_vec = density_vec,
    W = W,
    coef_names = coef_names,
    target_names = target_names
  )
  sandwich <- .plqr_sandwich_from_meat(
    jacobian = bread$jacobian,
    meat = moment_vcov,
    neff = bread$neff,
    coef_names = coef_names,
    target_names = target_names,
    meat_is_vcov = TRUE
  )

  c(
    sandwich,
    list(
      jacobian = bread$jacobian,
      neff = bread$neff,
      moment_vcov = moment_vcov,
      moment_se = moment_se
    )
  )
}

#' PLQR local moment evaluation and direct coefficient inference via DNN rank weights
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
#' Two interfaces are supported.
#'
#' 1. Data-frame interface (recommended): pass `data` with `x_cols`, `D_cols`,
#'    `g_col`, `mD_cols`, and `y_col`.
#' 2. Matrix interface: pass `X`, `D`, `g`, `mD`, and `Y`.
#'
#' No nuisance models are fit internally; this wrapper assumes the
#' partially-linear shape nuisance \eqn{g_0(U)} and the residualizing nuisance
#' \eqn{m_D(U) = E[D \mid U]} are supplied. The threshold is
#' \eqn{D'\theta + g_0(U)}; \eqn{g_0} alone is not itself a quantile function.
#' Both `g` and `mD` may depend on covariates beyond `x_cols` (the paper's
#' \eqn{U = (X, Z)}); the package does not restrict their domain, since they
#' are supplied pre-computed. Density values for direct coefficient inference
#' are also supplied by the user rather than estimated internally; they play
#' the role of the paper's local conditional density at the threshold,
#' entering only the coefficient Jacobian. See
#' `markdown_notes/design/15_plqr_variance_estimation_rationale.md` for the
#' full nuisance-notation mapping and for why coefficient-level inference
#' requires this direct-sandwich route rather than the jackknife SE used by
#' [cate_dnn()]/[np_dnn()]. Direct coefficient inference requires
#' `normalize = TRUE`, matching the theorem convention of weights summing to
#' one per target.
#'
#' With a density input, the coefficient sandwich supports two meat
#' estimators (`theta_meat`): the analytic score covariance (default) and
#' the full delete-1 jackknife covariance of the localized vector moment.
#' The bread is density-backed in both routes; only the meat changes. A
#' componentwise `moment_se` alone is never coefficient inference — see
#' `markdown_notes/design/15_plqr_variance_estimation_rationale.md` for the
#' route table.
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
#' component; the public `moment` and `moment_se` matrix fields remain anchored
#' to the default candidate. Candidate diagnostics include `ess`, `max_weight`,
#' `top_n`, `top5_weight_share`, and `degenerate`; `top_n` is fixed at 5 and
#' is independent of the estimator's `k`. Direct PLQR sandwich inference from
#' `density` or `density_col` currently uses `tuning_mode = "default"` only.
#'
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
#'   route and returns coefficient standard errors in `theta_se`.
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
#'   The `neighbors` extra returns the first `k` nearest observations used
#'   inside each size-`s` subsample, not the full DNN weight support.
#'
#' @return An object of class `plqr_result` with components:
#'   - `moment`: matrix (`p_d x m`) of localized PLQR moments;
#'   - `theta`: matrix (`p_d x m`) of the supplied coefficient values, repeated
#'     across targets until a root-solving path is requested in future versions;
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
#' demo$data$density_oracle <- 2
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
  method <- match.arg(method)
  jackknife <- match.arg(jackknife)
  theta_meat <- match.arg(theta_meat)
  tuning_mode <- match.arg(tuning_mode)
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
  score_mat <- .plqr_dnn_score_matrix(d_mat, mD_mat, g_vec, y_vec, theta_eval = theta_eval, tau = tau)
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
        precomputed_ranks = wt$ranks,
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
      precomputed_ranks = wt$ranks
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
      ranks = wt$ranks,
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
  out$.model <- list(
    x_mat = x_mat,
    d_mat = d_mat,
    mD_mat = mD_mat,
    distance = distance,
    method = method,
    k = k,
    s = s,
    normalize = normalize,
    ctx = ctx,
    ranks = if (identical(tuning_mode, "grid")) wt$ranks else NULL,
    d_resid_mat = d_resid_mat,
    below_threshold = as.numeric(y_vec <= as.vector(d_mat %*% theta_eval) + g_vec),
    density = density_vec,
    theta_eval = theta_eval,
    tau = tau,
    x_cols = if (!is.null(data)) x_cols else NULL,
    coef_names = coef_names
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
  ranks, score_mat, x_mat, x0_mat, distance, method, k, normalize, ctx, jackknife
) {
  force(ranks)
  force(score_mat)
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
        .dnn_loo_loop(score_mat, ranks, spec)
      } else {
        .dnn_loo_prefix(
          score_mat,
          ranks,
          tab_jk,
          tab_full = attr(tab_jk, "tab_full", exact = TRUE)
        )
      }
      .dnn_jackknife_se_from_loo(
        loo = loo,
        q = ncol(score_mat),
        n = nrow(ranks),
        m = ncol(ranks)
      )
    } else {
      .plqr_dnn_jackknife_se_loo(
        x_mat = x_mat,
        x0_mat = x0_mat,
        score_mat = score_mat,
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
  score_mat,
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
      score_mat = score_mat,
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
      moment <- .dnn_estimate_from_signal(score_mat, W, n = nrow(x_mat), arg = "score")
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
    p_d = ncol(score_mat),
    k = k,
    method = method,
    distance = distance,
    se = se,
    jackknife = jackknife
  )
  grid
}
