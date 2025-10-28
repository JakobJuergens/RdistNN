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
    errs <- c(errs, sprintf("`beta` length (%d) must equal length(D_cols) (%d).",
                            length(beta), length(D_cols)))
  }
  if (!is.numeric(j) || length(j) != 1L || j != as.integer(j) ||
      j < 1L || j > max(1L, length(D_cols))) {
    errs <- c(errs, sprintf("`j` must be an integer in 1..length(D_cols) (=%d).",
                            length(D_cols)))
  }

  if (length(errs)) stop(paste(errs, collapse = "  "), call. = FALSE)
  invisible(TRUE)
}

# R/score-plqr-component.R
#' PLQR diagnostic score (component j), pure score (no weights, no fitting)
#'
#' Computes per-observation contributions
#'   g_i = (tau - 1{ Y_i <= D_i^T beta + q_tau(X_i) }) * (D_{i,j} - m_{D,ij}),
#' stores them in `ctx$g`, and returns subset sums of `g` in `eval`/`eval_batch`.
#' No ordering or L-stat weights here; any aggregation/weights are external.
#'
#' @param x0 numeric length-p target point (checked for dimension only; not used).
#' @param x_cols names/indices of X columns in `data`.
#' @param D_cols names/indices of D columns in `data` (length p_d).
#' @param q_col name of `q_tau(X)` column in `data`.
#' @param mD_cols names/indices of `m_D(X)` columns (same length/order as `D_cols`).
#' @param beta numeric vector of length `length(D_cols)` (fixed coefficients).
#' @param j integer component index in `1:length(D_cols)`.
#' @param tau quantile level in (0,1).
#' @param s integer kernel order (reported via `score_order`; not used internally).
#' @param k integer kept for API compatibility (ignored by the pure score).
#' @return `ustat_score` with `ctx$g` and simple index/sum evals.
#' @export
score_plqr_component <- function(
    x0, x_cols, D_cols, q_col, mD_cols,
    beta, j = 1L, tau = 0.5,
    s, k = 1L
) {
  # ---- factory-time checks (order matters for tests) ----
  .validate_plqr_factory(D_cols, mD_cols, beta, j, tau, s, k)

  new_score(
    name  = paste0("plqr_comp_j", as.integer(j), "_purescore"),
    order = function(n, data, ctx) as.integer(s),

    setup = function(data, ...) {
      X   <- as.matrix(data[, x_cols,  drop = FALSE])
      D   <- as.matrix(data[, D_cols,  drop = FALSE])
      mD  <- as.matrix(data[, mD_cols, drop = FALSE])
      qx  <- as.numeric(data[[q_col]])
      Y   <- as.numeric(data$Y)

      # data-time checks
      if (length(x0) != ncol(X)) stop("x0 length must equal ncol(X).", call. = FALSE)
      if (!all(is.finite(X)))    stop("X must be finite.",               call. = FALSE)
      if (!all(is.finite(D)))    stop("D must be finite.",               call. = FALSE)
      if (!all(is.finite(mD)))   stop("m_D(X) must be finite.",          call. = FALSE)
      if (!all(is.finite(qx)))   stop("q_tau(X) must be finite.",        call. = FALSE)
      if (!all(is.finite(Y)))    stop("Y must be finite.",               call. = FALSE)

      n <- nrow(X)
      if (any(c(nrow(D), nrow(mD)) != n) || length(qx) != n || length(Y) != n) {
        stop("X, D, mD, q_tau(X), and Y must share the same n.", call. = FALSE)
      }
      if (ncol(D) != ncol(mD)) {
        stop("D and m_D(X) must have the same number of columns.", call. = FALSE)
      }
      if (length(beta) != ncol(D)) {
        stop(sprintf("`beta` length (%d) must equal ncol(D) (%d).",
                     length(beta), ncol(D)), call. = FALSE)
      }
      if (!(j >= 1L && j <= ncol(D))) {
        stop(sprintf("`j` must be in 1..ncol(D) (=%d).", ncol(D)), call. = FALSE)
      }

      # ---- pure-score payload: compute g and stash it ----
      thr   <- as.vector(D %*% as.numeric(beta)) + qx
      resid <- (as.numeric(tau) - as.numeric(Y <= thr))
      g     <- resid * (D[, j] - mD[, j])

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
      if (K == 0L) return(numeric(0L))
      mat <- matrix(gi[subs], nrow = K, ncol = ncol(subs))
      rowSums(mat)
    },

    randomized = FALSE,
    symmetric  = TRUE,
    ptr        = NULL
  )
}
