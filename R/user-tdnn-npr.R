# TDNN-based nonparametric regression wrapper (rank-weighted local mean)
#
# TDNN = Demirkaya two-scale DNN (k = 1 fixed).  The TDNN weight matrix is a
# signed linear combination of two single-scale DNN tables at orders s1 < s2.
# The combination coefficients cancel the leading bias term in the DNN MSE
# expansion, yielding a faster rate.  Negative weights are by design.
#
# This wrapper mirrors np_dnn() with (k, s) replaced by (s1, s2, dim).

#' Nonparametric regression via TDNN rank weights
#'
#' Estimate a localized nonparametric mean by combining observed outcomes with
#' Demirkaya two-scale DNN (TDNN) local weights.  TDNN is a signed linear
#' combination of two DNN tables at orders \code{s1 < s2}; the combination
#' cancels the leading bias term and achieves a faster MSE rate than
#' single-scale DNN under sufficient smoothness.
#'
#' @details
#' The documented sampling framework uses independent observations; cluster
#' and survey-design inference are not supplied. Exact distance ties use
#' lexicographic order of the original covariates, in their supplied column
#' order. Identical covariate rows retain input-order ties; permuting those
#' duplicates can change estimates. Ordering uses no outcomes or nuisances.
#'
#' Two interfaces are supported.
#'
#' 1. Data-frame interface (recommended): pass `data`, `x_cols`, and exactly
#'    one of `y_col` or `y`.
#' 2. Matrix interface: pass `X`, one of `X0` or `x0`, and `Y`.
#'
#' With `data`, omit `X` and `Y`; use exactly one of `y_col` or `y`.
#' Column selectors require `data`. Conflicting routes are rejected.
#'
#' **Signed weights.** The columns of the TDNN weight matrix sum to 1, but
#' individual entries can be negative.  This is expected and reflects the
#' bias-correction structure; it does not indicate numerical problems.
#'
#' **Localization dimension.** `dim` governs the coefficients combining the
#' two DNN scales. It defaults to `ncol(X)` (the covariate
#' dimension) and should be overridden only when a custom distance has a
#' different effective dimension.
#'
#' **Jackknife.** `"fast"` (default) uses the prefix/suffix rank-order engine
#' ([jackknife_rank_full()]), which precomputes the TDNN rank tables at `n`
#' and `n-1` once and runs an \eqn{O(n \cdot m)} loop without re-evaluating
#' the distance.  `"loo"` recomputes weights from scratch for each deletion
#' and is useful for validation.
#'
#' @param data Optional data frame (or tibble) containing covariates and the
#'   outcome column. Column names must be non-missing, non-empty, and unique.
#' @param x_cols Names or integer indices of covariate columns in `data`.
#' @param y_col Name or integer index of the outcome column in `data`.
#' @param y Optional numeric outcome vector of length `n`. In the data-frame
#'   interface, use this instead of `y_col`.
#' @param X Numeric matrix (`n x p`) for the matrix interface.
#' @param Y Numeric outcome vector of length `n` for the matrix interface.
#' @param X0 Numeric matrix (`m x p`) of target points; mutually exclusive
#'   with `x0`.
#' @param x0 Numeric vector of length `p` for a single target point.
#' @param distance Distance name or specification (for example `"l2"`).
#' @param method One of `"exact"` or `"asymp"` for TDNN rank weights.
#' @param s1 Positive integer; smaller U-statistic order.  Must satisfy
#'   \code{s1 < s2}.
#' @param s2 Positive integer; larger U-statistic order.  Must satisfy
#'   \code{s1 < s2 <= n}.
#' @param dim Positive integer; localization dimension.  Defaults to
#'   \code{ncol(X)} when \code{NULL}.
#' @param normalize Logical; renormalize the asymptotic rank table to sum 1
#'   (only used when `method = "asymp"`; default \code{TRUE}).
#' @param se Logical; compute delete-1 jackknife standard errors.
#' @param jackknife One of `"fast"`, `"loo"`, or `"none"`.  `"fast"` (default)
#'   uses the prefix/suffix rank-order engine; `"loo"` recomputes leave-one-out
#'   weights from scratch; `"none"` skips SEs entirely.
#' @param ctx Optional context passed through to distance and rank-weight
#'   evaluators.
#' @param return Optional extras in addition to the always-returned `mu`, `se`,
#'   and `meta`.  Allowed values:
#'   `c("W", "ranks", "ord", "dx", "y", "ess", "neighbors")`.
#'   The `neighbors` extra returns the nearest observation for each target; it
#'   is not the TDNN signed weight support.
#'
#' @return An object of class `tdnn_np_result` with components:
#'   - `mu`: estimated local mean values at target points;
#'   - `se`: optional jackknife standard errors;
#'   - `x0`: the `m x p` matrix of evaluation points, with covariate
#'     column names propagated from the training data when available;
#'   - `meta`: metadata (`n`, `m`, `p`, `s1`, `s2`, `dim`, `method`,
#'     `distance`);
#'   - `.model`: retained canonical training state and fitted geometry used
#'     by the inherited prediction and covariance methods;
#'   - optional extras requested via `return`.
#'
#' @seealso [weights_tdnn()], [np_dnn()], [np_demo_data()], [np_result] for
#'   inherited summary, covariance, confidence interval and prediction methods.
#'
#' @examples
#' set.seed(1)
#' n <- 80
#' x_mat <- matrix(rnorm(n * 2), n, 2)
#' y_vec <- x_mat[, 1] - 0.5 * x_mat[, 2] + rnorm(n, sd = 0.3)
#' out <- np_tdnn(
#'   X = x_mat, x0 = c(0, 0), Y = y_vec,
#'   distance = "l2", method = "exact", s1 = 2L, s2 = 4L,
#'   se = FALSE, return = c("W", "ess")
#' )
#' out$mu
#' @export
np_tdnn <- function(
  data     = NULL,
  x_cols   = NULL,
  y_col    = NULL,
  y        = NULL,
  X        = NULL,
  Y        = NULL,
  X0       = NULL,
  x0       = NULL,
  distance = "l2",
  method   = c("exact", "asymp"),
  s1       = 1L,
  s2       = 2L,
  dim      = NULL,
  normalize = TRUE,
  se       = TRUE,
  jackknife = c("fast", "loo", "none"),
  ctx      = NULL,
  return   = c("W")
) {
  # ------------------------------------------------------------------
  # Block 1: option validation (before any data is touched)
  # ------------------------------------------------------------------
  method    <- match.arg(method)
  jackknife <- match.arg(jackknife)
  s1 <- .validate_scalar_int(s1, "s1")
  s2 <- .validate_scalar_int(s2, "s2")
  if (s1 >= s2) {
    .abort_bad_arg("s1", sprintf("must satisfy s1 < s2 (got s1=%d, s2=%d)", s1, s2))
  }
  normalize <- .validate_scalar_logical(normalize, "normalize")
  se <- .validate_scalar_logical(se, "se")
  .warn_unknown_return(return, valid = c("W", "ranks", "ord", "dx", "y", "ess", "neighbors"))

  # ------------------------------------------------------------------
  # Block 2: dispatch and input resolution
  # ------------------------------------------------------------------
  .dnn_validate_input_route(data, list(X = X, Y = Y),
    list(x_cols = x_cols, y_col = y_col))

  if (!is.null(data)) {
    if (!is.data.frame(data)) {
      .abort_bad_arg("data", "must be a data.frame")
    }
    if (is.null(x_cols)) {
      .abort_bad_arg("x_cols", "is required with `data=` interface")
    }
    x_mat <- .dnn_col_mat(x_cols, data = data, arg = "x_cols", out_arg = "X")
    n     <- nrow(x_mat)
    p     <- ncol(x_mat)

    if (is.null(y_col) == is.null(y)) {
      .abort_bad_arg(
        "y_col/y",
        "must provide exactly one of `y_col` or `y` with `data=` interface",
        quoted = FALSE
      )
    }
    if (!is.null(y_col)) {
      y_vec <- .dnn_col_vec(y_col, data, n, arg = "y_col", ref = "X")
    } else {
      y_vec <- .validate_vector_nrow(y, n, arg = "y", ref = "X", strict = TRUE)
    }
  } else {
    if (is.null(X)) {
      .abort_bad_arg(
        "X",
        "must be provided for the matrix interface when `data` is NULL",
        expected = "non-NULL `X`, or use `data=` with `x_cols`"
      )
    }
    x_mat <- .validate_matrix_nonempty(X, arg = "X")
    n     <- nrow(x_mat)
    p     <- ncol(x_mat)

    if (!is.null(y)) {
      .abort_bad_arg("y", "is only supported with `data=`; matrix interface uses `Y`")
    }
    if (is.null(Y)) {
      .abort_bad_arg("Y", "is required for the matrix interface")
    }
    y_vec <- .validate_vector_nrow(Y, n, arg = "Y", ref = "X", strict = TRUE)
  }

  # ------------------------------------------------------------------
  # Block 3: canonicalize — resolve dim, jackknife feasibility, targets
  # ------------------------------------------------------------------
  if (is.null(dim)) {
    dim <- p
  } else {
    dim <- .validate_scalar_int(dim, "dim")
  }

  # Jackknife feasibility: delete-1 requires s2 < n so the same pair
  # (s1, s2) is valid at both sizes n and n-1.
  if (isTRUE(se) && jackknife != "none" && s2 >= n) {
    .abort_bad_arg(
      "s2",
      "must satisfy s2 < n when delete-1 jackknife is requested",
      expected = sprintf("s2 < n=%d", n),
      actual = sprintf("s2=%d (with `se`=%s, `jackknife`='%s')", s2, as.character(se), jackknife)
    )
  }

  x0_mat <- .dnn_resolve_targets(x_mat, X0 = X0, x0 = x0)

  # ------------------------------------------------------------------
  # Block 4: computation
  # ------------------------------------------------------------------

  # Request ranks on the first weights call so the fast jackknife can
  # reuse them without recomputing distances.
  geometry <- .dnn_fitted_geometry(distance, x_mat, ctx)
  weight_return <- return
  if (isTRUE(se) && identical(jackknife, "fast")) {
    weight_return <- unique(c(return, "ranks"))
  }

  weight_fit    <- .tdnn_weights_with_extras(
    x_mat     = x_mat,
    x0_mat    = x0_mat,
    distance  = geometry$distance,
    method    = method,
    s1        = s1,
    s2        = s2,
    dim       = dim,
    normalize = normalize,
    ctx       = geometry$ctx,
    return    = weight_return
  )
  wt            <- weight_fit$wt
  weight_matrix <- weight_fit$W
  m_out         <- weight_fit$m

  mu_hat <- as.numeric(
    .dnn_estimate_from_signal(y_vec, weight_matrix, n = n, arg = "Y")[1L, ]
  )

  se_hat <- NULL
  if (isTRUE(se) && jackknife != "none") {
    se_hat <- .np_tdnn_jackknife_se(
      x_mat             = x_mat,
      x0_mat            = x0_mat,
      y                 = y_vec,
      distance          = distance,
      method            = method,
      s1                = s1,
      s2                = s2,
      dim               = dim,
      normalize         = normalize,
      ctx               = ctx,
      jackknife         = jackknife,
      precomputed_ranks = wt$ranks
    )
  }

  # ------------------------------------------------------------------
  # Block 5: object assembly
  # ------------------------------------------------------------------
  out <- list(
    mu  = mu_hat,
    se  = se_hat,
    x0  = .dnn_x0_with_names(x0_mat, x_mat),
    meta = list(
      n        = n,
      m        = m_out,
      p        = p,
      s1       = s1,
      s2       = s2,
      dim      = dim,
      method   = method,
      distance = .dnn_distance_name(distance),
      jackknife = if (isTRUE(se)) jackknife else "none",
      tie_rule = "distance_then_covariates_then_input_position"
    )
  )
  # k = 1L is TDNN's fixed single-neighbor order; used only by the
  # "neighbors" extra in .dnn_attach_common_extras to cap the return size.
  out <- .dnn_attach_common_extras(out, return = return, wt = wt, W = weight_matrix, k = 1L)
  if ("y" %in% return) out$y <- y_vec

  out$.model <- c(list(
    kind = "tdnn", y_vec = y_vec, x_mat = x_mat,
    x_cols = if (!is.null(data)) .dnn_resolve_data_cols(x_cols, data, "x_cols") else NULL,
    k = 1L, s1 = s1, s2 = s2, dim = dim, method = method,
    normalize = normalize, ranks = wt$ranks
  ), geometry)

  class(out) <- c("tdnn_np_result", "np_result")
  out
}

# ---- internal jackknife dispatcher -------------------------------------------

# Thin dispatcher that routes to the fast prefix/suffix engine or the
# plain LOO recomputation, matching the signature expected by np_tdnn().
.np_tdnn_jackknife_se <- function(
  x_mat, x0_mat, y, distance, method, s1, s2, dim, normalize, ctx,
  jackknife         = c("fast", "loo"),
  precomputed_ranks = NULL
) {
  jackknife <- match.arg(jackknife)

  if (identical(jackknife, "fast")) {
    # Prefix/suffix engine: precomputes rank tables at n and n-1 once,
    # then runs an O(n * m) sweep without re-evaluating distances.
    as.numeric(
      .tdnn_jackknife_se_loo_full(
        x_mat             = x_mat,
        x0_mat            = x0_mat,
        signal            = y,
        distance          = distance,
        method            = method,
        s1                = s1,
        s2                = s2,
        dim               = dim,
        normalize         = normalize,
        ctx               = ctx,
        signal_arg        = "Y",
        precomputed_ranks = precomputed_ranks
      )[1L, ]
    )
  } else {
    # LOO fallback: recomputes weights from scratch for every deletion.
    as.numeric(
      .tdnn_jackknife_se_loo_matrix(
        x_mat             = x_mat,
        x0_mat            = x0_mat,
        signal            = y,
        distance          = distance,
        method            = method,
        s1                = s1,
        s2                = s2,
        dim               = dim,
        normalize         = normalize,
        ctx               = ctx,
        signal_arg        = "Y",
        jackknife         = "loo",
        precomputed_ranks = NULL
      )[1L, ]
    )
  }
}
