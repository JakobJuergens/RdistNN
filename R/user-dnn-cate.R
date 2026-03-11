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
#' 2. Matrix interface (backward compatible): pass `X`, one of `X0` or `x0`,
#'    and `pseudo_outcome`.
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
#' @param normalize Logical; renormalize columns of `W` to sum to 1.
#' @param se Logical; compute delete-1 jackknife standard errors.
#' @param jackknife One of `"fast"`, `"loo"`, or `"none"`. Currently `"fast"`
#'   and `"loo"` call the same delete-1 recomputation engine.
#' @param ctx Optional context passed to distance and ranking components.
#' @param return Optional extras in addition to the always-returned `tau`, `se`,
#'   and `meta`. Allowed values:
#'   `c("W", "ranks", "ord", "dx", "psi", "nuisance", "ess", "neighbors")`.
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
#'   distance = "l2", method = "exact", k = 2L, s = 3L,
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
  X = NULL, X0 = NULL, x0 = NULL, # matrix API (back-compat)
  distance = "l2",
  method = c("exact", "asymp"),
  k = 1L, s = 1L,
  normalize = TRUE,
  se = TRUE,
  jackknife = c("fast", "loo", "none"),
  ctx = NULL,
  return = c("W")
) {
  method <- match.arg(method)
  jackknife <- match.arg(jackknife)

  # Dispatch: data-frame API vs matrix API.
  if (!is.null(data)) {
    if (!is.data.frame(data)) {
      .abort_bad_arg("data", "must be a data.frame")
    }
    x_mat <- .validate_matrix_nonempty(data[, x_cols, drop = FALSE], arg = "X")
    n <- nrow(x_mat)
    p <- ncol(x_mat)

    if (!is.null(pseudo_outcome)) {
      psi <- .validate_vector_nrow(pseudo_outcome, n, arg = "`pseudo_outcome`", ref = "X")
      m0 <- m1 <- m <- e <- NULL
    } else {
      # Determine which outcome-model interface is used
      has_aipw <- !is.null(m0_col) && !is.null(m1_col)
      has_robinson <- !is.null(m_col)
      if (has_aipw && has_robinson) {
        stop("Provide either {m0_col, m1_col} (full AIPW) OR m_col (Robinson form), not both.",
          call. = FALSE)
      }
      if (!has_aipw && !has_robinson) {
        stop(
          "With `data=`, provide `pseudo_outcome` OR outcome model columns:\n",
          "  - preferred: {y_col, t_col, m0_col, m1_col, e_col} for full AIPW pseudo-outcome\n",
          "  - alternative: {y_col, t_col, m_col, e_col} for Robinson-form pseudo-outcome",
          call. = FALSE
        )
      }
      if (is.null(y_col) || is.null(t_col) || is.null(e_col)) {
        stop("With `data=`, `y_col`, `t_col`, and `e_col` are required.", call. = FALSE)
      }

      y_vec <- .validate_vector_nrow(data[[y_col]], n, arg = "y_col", ref = "X")
      t_vec <- .validate_vector_nrow(data[[t_col]], n, arg = "t_col", ref = "X")
      t_vec <- .validate_binary01(t_vec, arg = "t_col")
      e     <- .validate_vector_nrow(data[[e_col]], n, arg = "e_col", ref = "X")
      den   <- pmax(e * (1 - e), 1e-8)

      if (has_aipw) {
        # Full doubly-robust AIPW pseudo-outcome (paper formula, Algorithm 1 / Eq. 5.1):
        # φ̂_i = T_i(Y_i - μ̂₁) / (π̂_i(1-π̂_i)) - (1-T_i)(Y_i - μ̂₀) / (π̂_i(1-π̂_i)) + μ̂₁ - μ̂₀
        m0  <- .validate_vector_nrow(data[[m0_col]], n, arg = "m0_col", ref = "X")
        m1  <- .validate_vector_nrow(data[[m1_col]], n, arg = "m1_col", ref = "X")
        m   <- NULL
        psi <- (t_vec * (y_vec - m1) - (1 - t_vec) * (y_vec - m0)) / den + (m1 - m0)
      } else {
        # Robinson-form: (Y - m(X))(T - e) / (e(1-e))
        m0 <- m1 <- NULL
        m   <- .validate_vector_nrow(data[[m_col]], n, arg = "m_col", ref = "X")
        psi <- (y_vec - m) * (t_vec - e) / den
      }
    }
  } else {
    # matrix API
    if (is.null(X)) stop("Provide either `data` ... or matrix interface via `X`.", call. = FALSE)
    x_mat <- .validate_matrix_nonempty(X, arg = "X")
    n <- nrow(x_mat)
    p <- ncol(x_mat)
    if (is.null(pseudo_outcome)) stop("Matrix API requires `pseudo_outcome`.", call. = FALSE)
    psi <- .validate_vector_nrow(pseudo_outcome, n, arg = "`pseudo_outcome`", ref = "X")
    m0 <- m1 <- m <- e <- NULL
  }

  # Targets
  if (!is.null(X0) && !is.null(x0)) stop("Provide either `X0` or `x0`, not both.", call. = FALSE)
  if (is.null(X0)) {
    if (is.null(x0)) stop("Provide `X0` (m x p) or a single `x0` (length p).", call. = FALSE)
    x0_vec <- .validate_target_vector(x0, p, x0_arg = "x0", x_arg = "X")
    x0_mat <- matrix(x0_vec, nrow = 1L)
  } else {
    x0_mat <- .validate_ncol_relation(x_mat, X0, y_arg = "X0", x_arg = "X", relation = "equal")
  }

  # Weights (rank-only DNN).
  keep_extras <- intersect(return, c("ranks", "ord", "dx"))
  wt <- weights_dnn(
    X = x_mat, X0 = x0_mat,
    distance = distance,
    k = as.integer(k), s = as.integer(s),
    method = method, normalize = normalize, ctx = ctx,
    return = unique(c("W", keep_extras))
  )
  weight_matrix <- wt$W
  m_out <- ncol(weight_matrix)

  # Point estimates.
  tau <- as.numeric(crossprod(psi, weight_matrix)) # length m_out
  se_hat <- NULL
  if (isTRUE(se) && jackknife != "none") {
    if (as.integer(s) >= n) {
      stop("Delete-1 jackknife requires `s < n` so the same order can be used after deletion.", call. = FALSE)
    }
    se_hat <- .cate_dnn_jackknife_se_loo(
      x_mat = x_mat,
      x0_mat = x0_mat,
      psi = psi,
      distance = distance,
      method = method,
      k = as.integer(k),
      s = as.integer(s),
      normalize = normalize,
      ctx = ctx
    )
  }

  # Assemble output.
  distance_name <- if (is.character(distance)) {
    distance
  } else if (is.list(distance) && !is.null(distance$name)) {
    distance$name
  } else {
    "user_distance"
  }

  out <- list(
    tau = tau,
    se = se_hat,
    meta = list(
      n = n, m = m_out, p = p,
      k = as.integer(k), s = as.integer(s),
      method = method, distance = distance_name
    )
  )
  if ("W" %in% return) out$W <- weight_matrix
  if ("ranks" %in% return && !is.null(wt$ranks)) out$ranks <- wt$ranks
  if ("ord" %in% return && !is.null(wt$ord)) out$ord <- wt$ord
  if ("dx" %in% return && !is.null(wt$dx)) out$dx <- wt$dx
  if ("psi" %in% return) out$psi <- psi
  if ("nuisance" %in% return) {
    if (!is.null(m0) && !is.null(m1) && !is.null(e)) {
      out$nuisance <- list(m0 = m0, m1 = m1, e = e)
    } else if (!is.null(m) && !is.null(e)) {
      out$nuisance <- list(m = m, e = e)
    }
  }
  if ("ess" %in% return) {
    out$ess <- 1 / colSums(weight_matrix^2)
  }
  if ("neighbors" %in% return) {
    k_eff <- max(1L, min(nrow(weight_matrix), as.integer(k)))
    if (!is.null(wt$ord)) {
      out$neighbors <- wt$ord[seq_len(k_eff), , drop = FALSE]
    } else if (!is.null(wt$ranks)) {
      out$neighbors <- apply(wt$ranks, 2L, function(r) order(r)[seq_len(k_eff)])
    }
  }

  class(out) <- "cate_result"
  out
}

# Internal: paper-consistent delete-1 jackknife for cate_dnn().
# Implements the delete-1 recomputation in Eq. (3.30)-(3.32) of
# Local_Z_Est_JJuergens.pdf for the scalar (p_alpha = 1) CATE target.
.cate_dnn_jackknife_se_loo <- function(
  x_mat, x0_mat, psi, distance, method, k, s, normalize, ctx
) {
  n <- nrow(x_mat)
  m <- nrow(x0_mat)
  tau_loo <- matrix(NA_real_, nrow = n, ncol = m)

  for (i in seq_len(n)) {
    keep <- (seq_len(n) != i)
    wt_i <- weights_dnn(
      X = x_mat[keep, , drop = FALSE],
      X0 = x0_mat,
      distance = distance,
      k = k,
      s = s,
      method = method,
      normalize = normalize,
      ctx = ctx,
      return = "W"
    )
    tau_loo[i, ] <- as.numeric(crossprod(psi[keep], wt_i$W))
  }

  tau_loo_bar <- colMeans(tau_loo)
  centered <- sweep(tau_loo, 2L, tau_loo_bar, FUN = "-")
  var_jk <- ((n - 1) / n) * colSums(centered^2)
  sqrt(pmax(var_jk, 0))
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

#' Print a `cate_dnn_result` object
#'
#' @inheritParams print.cate_result
#' @return `x`, invisibly.
#' @export
print.cate_dnn_result <- function(x, ...) {
  cat("<cate_dnn_result>\n", sep = "")
  cat("  data  : n=", x$meta$n, ", p=", x$meta$p, ", m=", x$meta$m, "\n", sep = "")
  cat("  dnn   : k=", x$meta$k, ", s=", x$meta$s, ", method=", x$meta$method, "\n", sep = "")
  cat("  rank  : ", x$meta$weight_name, "; distance: ", x$meta$distance_name, "\n", sep = "")
  cat("  components: ", paste(names(x), collapse = ", "), "\n", sep = "")
  invisible(x)
}
