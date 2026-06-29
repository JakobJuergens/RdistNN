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
#' deletion.
#'
#' @param data Data frame with covariates and optional nuisance columns.
#' @param x_cols Names or indices of covariate columns in `data`.
#' @param y_col Name or index of the outcome column in `data`.
#' @param t_col Name or index of the treatment column in `data` (must be 0/1).
#' @param m0_col Name or index for \eqn{E[Y(0)\mid X]} (full AIPW interface).
#' @param m1_col Name or index for \eqn{E[Y(1)\mid X]} (full AIPW interface).
#' @param m_col Name or index for \eqn{E[Y\mid X]} (Robinson interface).
#' @param e_col Name or index for propensity \eqn{P(T=1\mid X)}.
#' @param pseudo_outcome Numeric vector of length `n`. If supplied, nuisance
#'   columns are not required.
#' @param X Numeric matrix (`n x p`) for the matrix interface.
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
#' @param return Optional extras in addition to the always-returned `tau`, `se`,
#'   and `meta`. Allowed values:
#'   `c("W", "ranks", "ord", "dx", "psi", "nuisance", "ess", "neighbors", "tuning")`.
#'
#' @return An object of class `cate_result` with components:
#'   - `tau`: estimated CATE values at target points;
#'   - `se`: optional jackknife standard errors;
#'   - `meta`: metadata (`n`, `m`, `p`, `k`, `s`, `method`, `distance`);
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
  distance = "l2",
  method = c("exact", "asymp"),
  k = 1L, s = 1L,
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

  .warn_unknown_return(return, valid = c("W", "ranks", "ord", "dx", "psi", "nuisance", "ess", "neighbors", "tuning"))
  t_vec <- e <- m0 <- m1 <- m <- NULL

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
      e     <- .dnn_col_vec(e_col,  data, n, arg = "e_col",  ref = "data")
      den   <- pmax(e * (1 - e), 1e-8)

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
          ehat = e
        )
      } else {
        # Robinson-form: (Y - m(X))(T - e) / (e(1-e))
        m0 <- m1 <- NULL
        m   <- .dnn_col_vec(m_col,  data, n, arg = "m_col",  ref = "data")
        psi <- cate_rho_rlearner(
          y = y_vec,
          d = t_vec,
          ehat = e,
          mhat = m
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
    meta = list(
      n = n, m = m_out, p = p,
      k = as.integer(k), s = as.integer(s),
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
  if ("psi" %in% return) out$psi <- psi
  if ("nuisance" %in% return) {
    if (!is.null(m0) && !is.null(m1) && !is.null(e)) {
      out$nuisance <- list(m0 = m0, m1 = m1, e = e)
    } else if (!is.null(m) && !is.null(e)) {
      out$nuisance <- list(m = m, e = e)
    }
  }

  # Internal fit state used by predict.cate_result().
  out$.model <- list(
    x_mat = x_mat,
    psi = psi,
    t_vec = t_vec,
    e = e,
    x_cols = if (!is.null(data)) x_cols else NULL,
    distance = distance,
    method = method,
    k = k,
    s = s,
    normalize = normalize,
    ctx = ctx
  )

  class(out) <- "cate_result"
  out
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
