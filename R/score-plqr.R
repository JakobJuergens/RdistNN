# R/score-plqr-component.R

.validate_plqr_factory <- function(D_cols, mD_cols, beta, j, tau, s, k) {
  errs <- character()

  # scalar/basic
  if (!is.numeric(tau) || length(tau) != 1L || !(tau > 0 && tau < 1)) {
    errs <- c(errs, "`tau` must be a numeric scalar in (0,1).")
  }
  if (!is.numeric(s) || length(s) != 1L || s < 1L) {
    errs <- c(errs, "`s` must be a positive integer.")
  }
  if (!is.numeric(k) || length(k) != 1L || k < 1L) {
    errs <- c(errs, "`k` must be a positive integer.")
  }

  # shapes by interface (no data needed)
  if (length(D_cols) < 1L) {
    errs <- c(errs, "`D_cols` must be a non-empty vector of columns.")
  }
  if (length(mD_cols) != length(D_cols)) {
    errs <- c(errs, "`mD_cols` must have the same length and order as `D_cols`.")
  }
  if (!is.numeric(beta) || length(beta) != length(D_cols)) {
    errs <- c(errs, sprintf(
      "`beta` length (%d) must equal length(D_cols) (%d).",
      length(beta), length(D_cols)
    ))
  }
  if (!is.numeric(j) || length(j) != 1L || j != as.integer(j) ||
    j < 1L || j > max(1L, length(D_cols))) {
    errs <- c(errs, sprintf(
      "`j` must be an integer in 1..length(D_cols) (=%d).",
      length(D_cols)
    ))
  }

  if (length(errs)) stop(paste(errs, collapse = "  "), call. = FALSE)
  invisible(TRUE)
}

#' PLQR diagnostic score (component \eqn{j}); pure score (no weights, no fitting)
#'
#' Constructs a **pure-score** kernel for PLQR diagnostics that stores the
#' per-observation contributions
#' \deqn{
#'   g_{i}
#'   \;=\;
#'   \Bigl(\tau - \mathbf{1}\{\,Y_{i} \le D_{i}^{\top}\beta + q_{\tau}(X_{i})\,\}\Bigr)\,
#'   \bigl(D_{i,j} - m_{D,ij}(X_{i})\bigr)\,,
#' }
#' and returns subset sums \eqn{\sum_{i \in S} g_i} in \code{eval}/\code{eval_batch}.
#' This score performs **no** nearest-neighbor ordering or L-stat weighting;
#' any aggregation/weights are handled outside the score layer.
#'
#' @param x0 Numeric length-\eqn{p} target point. Used only for basic dimension
#'   checks against \code{x_cols}; not used in the kernel.
#' @param x_cols Integer or character indices of the \eqn{X} columns in
#'   \code{data}; used for dimension checks with \code{x0}.
#' @param D_cols Integer or character indices of the \eqn{D} columns in
#'   \code{data} (length \eqn{p_d}).
#' @param q_col Character; name of the \eqn{q_{\tau}(X)} column in \code{data}.
#' @param mD_cols Integer or character indices of the \eqn{m_{D}(X)} columns in
#'   \code{data}; must have the same length and order as \code{D_cols}. The
#'   \eqn{j}-th column represents \eqn{m_{D,ij}(X)}.
#' @param beta Numeric vector of length \code{length(D_cols)} giving the fixed
#'   coefficient vector \eqn{\beta}.
#' @param j Integer component index in \eqn{1,\ldots,\mathrm{length}(D\_cols)}.
#' @param tau Quantile level in \eqn{(0,1)}.
#' @param s Integer kernel order (subset size). Recorded in the score metadata;
#'   not used internally by this pure score.
#' @param k Integer kept for API compatibility with other scores; ignored here.
#'
#' @return An \code{ustat_score} object whose \code{$setup} stores the vector
#'   \code{g} of contributions in the context, and whose \code{$eval}/\code{$eval_batch}
#'   return subset sums of \code{g} (no weights applied).
#'
#' @examples
#' demo <- plqr_demo_data(n = 200, p = 3, p_d = 2, tau = 0.5, seed = 1)
#' sc <- score_plqr_component(
#'   x0      = rep(0, 3),
#'   x_cols  = paste0("X", 1:3),
#'   D_cols  = paste0("D", 1:2),
#'   q_col   = "q_oracle",
#'   mD_cols = paste0("mD", 1:2),
#'   beta    = demo$truth$beta_true,
#'   j       = 1L,
#'   tau     = 0.5,
#'   s       = 5,   # recorded only
#'   k       = 1    # ignored here
#' )
#' ctx  <- score_setup(sc, demo$data)
#' subs <- cate_demo_subsets(n = nrow(demo$data), K = 3, s = 5, seed = 1)
#' score_eval_batch(sc, subs, demo$data, ctx)
#'
#' @export
score_plqr_component <- function(
  x0, x_cols, D_cols, q_col, mD_cols,
  beta, j = 1L, tau = 0.5,
  s, k = 1L
) {
  # ---- factory-time checks (order matters for tests) ----
  .validate_plqr_factory(D_cols, mD_cols, beta, j, tau, s, k)

  new_score(
    name = paste0("plqr_comp_j", as.integer(j), "_purescore"),
    order = function(n, data, ctx) as.integer(s),
    setup = function(data, ...) {
      X <- as.matrix(data[, x_cols, drop = FALSE])
      D <- as.matrix(data[, D_cols, drop = FALSE])
      mD <- as.matrix(data[, mD_cols, drop = FALSE])
      qx <- as.numeric(data[[q_col]])
      Y <- as.numeric(data$Y)

      # data-time checks
      if (length(x0) != ncol(X)) stop("x0 length must equal ncol(X).", call. = FALSE)
      if (!all(is.finite(X))) stop("X must be finite.", call. = FALSE)
      if (!all(is.finite(D))) stop("D must be finite.", call. = FALSE)
      if (!all(is.finite(mD))) stop("m_D(X) must be finite.", call. = FALSE)
      if (!all(is.finite(qx))) stop("q_tau(X) must be finite.", call. = FALSE)
      if (!all(is.finite(Y))) stop("Y must be finite.", call. = FALSE)

      n <- nrow(X)
      if (any(c(nrow(D), nrow(mD)) != n) || length(qx) != n || length(Y) != n) {
        stop("X, D, mD, q_tau(X), and Y must share the same n.", call. = FALSE)
      }
      if (ncol(D) != ncol(mD)) {
        stop("D and m_D(X) must have the same number of columns.", call. = FALSE)
      }
      if (length(beta) != ncol(D)) {
        stop(sprintf(
          "`beta` length (%d) must equal ncol(D) (%d).",
          length(beta), ncol(D)
        ), call. = FALSE)
      }
      if (!(j >= 1L && j <= ncol(D))) {
        stop(sprintf("`j` must be in 1..ncol(D) (=%d).", ncol(D)), call. = FALSE)
      }

      # ---- pure-score payload: compute g and stash it ----
      thr <- as.vector(D %*% as.numeric(beta)) + qx
      resid <- (as.numeric(tau) - as.numeric(Y <= thr))
      g <- resid * (D[, j] - mD[, j])

      if (!all(is.finite(g))) stop("Computed g contains non-finite values.", call. = FALSE)

      list(
        g = as.numeric(g),
        n = n, p = ncol(X),
        s = as.integer(s), j = as.integer(j),
        # compiled path intentionally omitted in pure layer
        ptr_ctx = NULL
      )
    },

    # R fallback: one subset → sum of g over that subset
    eval = function(idxs, data, ctx, rng) {
      gi <- ctx$g[idxs]
      if (any(!is.finite(gi))) stop("Non-finite g encountered.", call. = FALSE)
      as.numeric(sum(gi))
    },

    # R fast path: many subsets → rowwise sums of g[subs]
    eval_batch = function(subs, data, ctx, rng) {
      gi <- ctx$g
      if (any(!is.finite(gi))) stop("Non-finite g encountered.", call. = FALSE)
      K <- nrow(subs)
      if (K == 0L) {
        return(numeric(0L))
      }
      mat <- matrix(gi[subs], nrow = K, ncol = ncol(subs))
      rowSums(mat)
    },
    randomized = FALSE,
    symmetric = TRUE,
    ptr = NULL
  )
}
