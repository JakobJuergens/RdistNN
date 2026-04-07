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
#' deletion.
#'
#' @param data Optional data frame containing covariates, regressors, outcome,
#'   and nuisance columns.
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
#' @param k,s Integer DNN tuning parameters.
#' @param normalize Logical; renormalize columns of `W` to sum to 1.
#' @param se Logical; compute delete-1 jackknife standard errors.
#' @param jackknife One of `"fast"`, `"loo"`, or `"none"`. `"fast"` (default)
#'   uses cached full-sample ranks with delete-1 rank updates, `"loo"`
#'   recomputes leave-one-out weights directly, and `"none"` skips SEs.
#' @param ctx Optional context passed to distance and ranking components.
#' @param return Optional extras in addition to the always-returned `moment`,
#'   `se`, `beta`, `tau`, and `meta`. Allowed values:
#'   `c("W", "ranks", "ord", "dx", "g", "ess", "neighbors")`.
#'
#' @return An object of class `plqr_result` with components:
#'   - `moment`: matrix (`p_d x m`) of localized PLQR moments;
#'   - `se`: optional jackknife standard errors (`p_d x m`);
#'   - `beta`: the coefficient vector used in the moment map;
#'   - `tau`: quantile level;
#'   - `meta`: metadata (`n`, `m`, `p`, `p_d`, `k`, `s`, `method`, `distance`);
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
#'   k = 2L, s = 3L,
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
  s = 1L,
  normalize = TRUE,
  se = TRUE,
  jackknife = c("fast", "loo", "none"),
  ctx = NULL,
  return = c("W")
) {
  method <- match.arg(method)
  jackknife <- match.arg(jackknife)
  k <- .validate_scalar_int(k, "k")
  s <- .validate_scalar_int(s, "s")
  if (k > s) .abort_bad_arg("k", sprintf("must satisfy k <= s (got k=%d, s=%d)", k, s))
  .warn_unknown_return(return, valid = c("W", "ranks", "ord", "dx", "g", "ess", "neighbors"))
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
      stop(
        "With `data=`, provide `x_cols`, `D_cols`, `q_col`, `mD_cols`, and `y_col`.",
        call. = FALSE
      )
    }

    x_mat <- .validate_matrix_nonempty(data[, x_cols, drop = FALSE], arg = "X")
    n <- nrow(x_mat)
    p <- ncol(x_mat)

    d_mat <- .validate_matrix_nonempty(data[, D_cols, drop = FALSE], arg = "D")
    mD_mat <- .validate_matrix_nonempty(data[, mD_cols, drop = FALSE], arg = "mD")
    q_vec <- .validate_vector_nrow(data[[q_col]], n, arg = "q_col", ref = "X")
    y_vec <- .validate_vector_nrow(data[[y_col]], n, arg = "y_col", ref = "X")
  } else {
    if (is.null(X)) {
      stop("Provide either `data` ... or matrix interface via `X`.", call. = FALSE)
    }
    if (is.null(D) || is.null(q) || is.null(mD) || is.null(Y)) {
      stop("Matrix API requires `D`, `q`, `mD`, and `Y`.", call. = FALSE)
    }

    x_mat <- .validate_matrix_nonempty(X, arg = "X")
    n <- nrow(x_mat)
    p <- ncol(x_mat)

    d_mat <- .validate_matrix_nonempty(D, arg = "D")
    mD_mat <- .validate_matrix_nonempty(mD, arg = "mD")
    q_vec <- .validate_vector_nrow(q, n, arg = "q", ref = "X")
    y_vec <- .validate_vector_nrow(Y, n, arg = "Y", ref = "X")
  }

  if (nrow(d_mat) != n || nrow(mD_mat) != n) {
    stop("`D` and `mD` must have nrow equal to nrow(X).", call. = FALSE)
  }
  if (ncol(d_mat) != ncol(mD_mat)) {
    stop("`D` and `mD` must have the same number of columns.", call. = FALSE)
  }
  p_d <- ncol(d_mat)
  if (length(beta) != p_d) {
    stop(sprintf("`beta` length (%d) must equal ncol(D) (%d).", length(beta), p_d), call. = FALSE)
  }

  if (!all(is.finite(d_mat))) stop("D must be finite.", call. = FALSE)
  if (!all(is.finite(mD_mat))) stop("m_D(X) must be finite.", call. = FALSE)
  if (!all(is.finite(q_vec))) stop("q_tau(X) must be finite.", call. = FALSE)
  if (!all(is.finite(y_vec))) stop("Y must be finite.", call. = FALSE)

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
    meta = list(
      n = n, m = m_out, p = p, p_d = p_d,
      k = k, s = s,
      method = method, distance = .dnn_distance_name(distance)
    )
  )
  out <- .dnn_attach_common_extras(out, return = return, wt = wt, W = weight_matrix, k = k)
  if ("g" %in% return) out$g <- g_mat

  class(out) <- "plqr_result"
  out
}
