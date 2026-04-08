.validate_plqr_factory <- function(D_cols, mD_cols, beta, j, tau, s, k) {
  errs <- character()

  # scalar/basic
  tryCatch(.validate_scalar_prob(tau, arg = "tau"), error = function(e) errs <<- c(errs, conditionMessage(e)))
  tryCatch(.validate_scalar_int(s, arg = "s", min = 1L), error = function(e) errs <<- c(errs, conditionMessage(e)))
  tryCatch(.validate_scalar_int(k, arg = "k", min = 1L), error = function(e) errs <<- c(errs, conditionMessage(e)))

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
#' @examplesIf TRUE
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
#'   s       = 5,
#'   k       = 1
#' )
#' ctx <- score_setup(sc, demo$data)
#' subs <- cate_demo_subsets(n = nrow(demo$data), K = 3, s = 5, seed = 1)
#' score_eval_batch(sc, subs, demo$data, ctx)
#'
#' @export
score_plqr_component <- function(
  x0, x_cols, D_cols, q_col, mD_cols,
  beta, j = 1L, tau = 0.5,
  s, k = 1L
) {
  # ---- factory-time checks ----
  .validate_plqr_factory(D_cols, mD_cols, beta, j, tau, s, k)

  new_score(
    name = paste0("plqr_comp_j", as.integer(j), "_purescore"),
    order = function(n, data, ctx) as.integer(s),
    setup = function(data, ...) {
      X <- .dnn_col_mat(x_cols, data = data, arg = "x_cols", out_arg = "X")
      n <- nrow(X)
      D <- .dnn_col_mat(D_cols, data = data, arg = "D_cols", out_arg = "D")
      mD <- .dnn_col_mat(mD_cols, data = data, arg = "mD_cols", out_arg = "m_D(X)")
      qx <- .dnn_col_vec(q_col, data = data, n = n, arg = "q_col", ref = "data")
      Y <- .dnn_col_vec("Y", data = data, n = n, arg = "Y", ref = "data")

      # data-time checks
      .validate_target_vector(x0, ncol(X), x0_arg = "x0", x_arg = "X")
      if (.has_na_only(X)) stop("X cannot contain `NA` values.", call. = FALSE)
      if (.has_nonfinite_non_na(X)) {
        stop("X cannot contain non-finite values (`Inf`, `-Inf`, or `NaN`).", call. = FALSE)
      }
      if (.has_na_only(D)) stop("D cannot contain `NA` values.", call. = FALSE)
      if (.has_nonfinite_non_na(D)) {
        stop("D cannot contain non-finite values (`Inf`, `-Inf`, or `NaN`).", call. = FALSE)
      }
      if (.has_na_only(mD)) stop("m_D(X) cannot contain `NA` values.", call. = FALSE)
      if (.has_nonfinite_non_na(mD)) {
        stop("m_D(X) cannot contain non-finite values (`Inf`, `-Inf`, or `NaN`).", call. = FALSE)
      }
      if (.has_na_only(qx)) stop("q_tau(X) cannot contain `NA` values.", call. = FALSE)
      if (.has_nonfinite_non_na(qx)) {
        stop("q_tau(X) cannot contain non-finite values (`Inf`, `-Inf`, or `NaN`).", call. = FALSE)
      }
      if (.has_na_only(Y)) stop("Y cannot contain `NA` values.", call. = FALSE)
      if (.has_nonfinite_non_na(Y)) stop("Y cannot contain non-finite values (`Inf`, `-Inf`, or `NaN`).", call. = FALSE)

      if (any(c(nrow(D), nrow(mD)) != n) || length(qx) != n || length(Y) != n) {
        stop("X, D, mD, q_tau(X), and Y must share the same n.", call. = FALSE)
      }
      if (ncol(D) != ncol(mD)) {
        stop("D and m_D(X) must have the same number of columns.", call. = FALSE)
      }
      if (length(beta) != ncol(D)) {
        stop(sprintf("`beta` length (%d) must equal ncol(D) (%d).", length(beta), ncol(D)), call. = FALSE)
      }
      if (!(j >= 1L && j <= ncol(D))) {
        stop(sprintf("`j` must be in 1..ncol(D) (=%d).", ncol(D)), call. = FALSE)
      }

      # ---- compute and stash g ----
      g <- NULL
      ptr_ctx <- NULL
      if (exists("rustat_plqr_ctx_create", mode = "function") &&
        exists("rustat_plqr_g_all", mode = "function")) {
        # CORRECT argument order + include j
        ptr_ctx <- rustat_plqr_ctx_create(
          Y, # NumericVector
          D, # NumericMatrix
          qx, # NumericVector
          mD, # NumericMatrix
          as.numeric(beta),
          as.numeric(tau),
          as.integer(j)
        )
        g <- as.numeric(rustat_plqr_g_all(ptr_ctx))
      } else {
        # R fallback
        thr <- as.vector(D %*% as.numeric(beta)) + qx
        resid <- (as.numeric(tau) - as.numeric(Y <= thr))
        g <- resid * (D[, j] - mD[, j])
      }

      if (.has_na_only(g)) stop("Computed g contains `NA` values.", call. = FALSE)
      if (.has_nonfinite_non_na(g)) {
        stop("Computed g contains non-finite values (`Inf`, `-Inf`, or `NaN`).", call. = FALSE)
      }

      list(
        g = as.numeric(g),
        n = nrow(X), p = ncol(X),
        s = as.integer(s), j = as.integer(j),
        ptr_ctx = ptr_ctx
      )
    },
    eval = function(idxs, data, ctx, rng) {
      gi <- ctx$g[idxs]
      if (.has_na_only(gi)) stop("g contains `NA` values.", call. = FALSE)
      if (.has_nonfinite_non_na(gi)) stop("g contains non-finite values (`Inf`, `-Inf`, or `NaN`).", call. = FALSE)
      as.numeric(sum(gi))
    },
    eval_batch = function(subs, data, ctx, rng) {
      gi <- ctx$g
      if (.has_na_only(gi)) stop("g contains `NA` values.", call. = FALSE)
      if (.has_nonfinite_non_na(gi)) stop("g contains non-finite values (`Inf`, `-Inf`, or `NaN`).", call. = FALSE)
      if (!is.matrix(subs)) subs <- as.matrix(subs)
      if (nrow(subs) == 0L) {
        return(numeric(0L))
      }

      if (exists("rustat_subset_row_sums_g", mode = "function")) {
        return(rustat_subset_row_sums_g(subs, gi))
      }
      rowSums(matrix(gi[subs], nrow = nrow(subs), ncol = ncol(subs)))
    },
    randomized = FALSE,
    symmetric = TRUE,
    ptr = NULL
  )
}
