.validate_plqr_factory <- function(D_cols, mD_cols, theta, j, tau, s) {
  errs <- character()

  # scalar/basic
  tryCatch(.validate_scalar_prob(tau, arg = "tau"),
    rdistnn_input_error = function(e) errs <<- c(errs, conditionMessage(e)))
  tryCatch(.validate_scalar_int(s, arg = "s", min = 1L),
    rdistnn_input_error = function(e) errs <<- c(errs, conditionMessage(e)))

  # shapes by interface (no data needed)
  if (length(D_cols) < 1L) {
    errs <- c(errs, "`D_cols` must be a non-empty vector of columns.")
  }
  if (length(mD_cols) != length(D_cols)) {
    errs <- c(errs, "`mD_cols` must have the same length and order as `D_cols`.")
  }
  if (!is.numeric(theta) || length(theta) != length(D_cols)) {
    errs <- c(errs, sprintf(
      "`theta` length (%d) must equal length(D_cols) (%d).",
      length(theta), length(D_cols)
    ))
  } else {
    tryCatch(.abort_if_na_or_nonfinite(theta, arg = "theta", quoted = TRUE),
      rdistnn_input_error = function(e) errs <<- c(errs, conditionMessage(e)))
  }
  if (!is.numeric(j) || length(j) != 1L || !is.finite(j) || j != floor(j) ||
    j < 1L || j > max(1L, length(D_cols))) {
    errs <- c(errs, sprintf(
      "`j` must be an integer in 1..length(D_cols) (=%d).",
      length(D_cols)
    ))
  }

  if (length(errs)) .dnn_abort(
    "argument_invalid",
    paste(errs, collapse = "  "),
    argument = "D_cols/mD_cols/theta/j/tau/s"
  )
  invisible(TRUE)
}

#' PLQR diagnostic score (component \eqn{j}); pure score (no weights, no fitting) (experimental)
#'
#' Constructs a **pure-score** kernel for PLQR diagnostics that stores the
#' per-observation score contributions
#' \deqn{
#'   \psi_{i}
#'   \;=\;
#'   \Bigl(\tau - \mathbf{1}\{\,Y_{i} \le D_{i}^{\top}\theta + g_{0}(U_{i})\,\}\Bigr)\,
#'   \bigl(D_{i,j} - m_{D,ij}(U_{i})\bigr)\,,
#' }
#' and returns subset sums \eqn{\sum_{i \in S} \psi_i} in \code{eval}/\code{eval_batch}.
#' This score performs **no** nearest-neighbor ordering or L-stat weighting;
#' any aggregation/weights are handled outside the score layer.
#'
#' @section Lifecycle:
#' **Experimental.** The PLQR API and inference interfaces may change.
#' @param x0 Numeric length-\eqn{p} target point. Used only for basic dimension
#'   checks against \code{x_cols}; not used in the kernel.
#' @param x_cols Integer or character indices of the \eqn{X} columns in
#'   \code{data}; used for dimension checks with \code{x0}.
#' @param D_cols Integer or character indices of the \eqn{D} columns in
#'   \code{data} (length \eqn{p_d}).
#' @param g_col Character; name of the column in \code{data} holding the
#'   shape-nuisance values \eqn{g_0(U)}.
#' @param mD_cols Integer or character indices of the \eqn{m_{D}(U)} columns in
#'   \code{data}; must have the same length and order as \code{D_cols}. The
#'   \eqn{j}-th column represents \eqn{m_{D,ij}(U)}.
#' @param y_col Name or index of the outcome column in \code{data}.
#' @param theta Numeric vector of length \code{length(D_cols)} giving the fixed
#'   coefficient vector \eqn{\theta}.
#' @param j Integer component index in \eqn{1,\ldots,\mathrm{length}(D\_cols)}.
#' @param tau Quantile level in \eqn{(0,1)}.
#' @param s Integer kernel order (subset size). Recorded in the score metadata;
#'   not used internally by this pure score.
#'
#' @return An \code{rdistnn_ustat_score} object whose \code{$setup} stores the vector
#'   \code{score} of contributions in the context, and whose
#'   \code{$eval}/\code{$eval_batch} return subset sums of \code{score}
#'   (no weights applied).
#'
#' @examples
#' demo <- plqr_demo_data(n = 200, p = 3, p_d = 2, tau = 0.5, seed = 1)
#' sc <- score_plqr_component(
#'   x0      = rep(0, 3),
#'   x_cols  = paste0("X", 1:3),
#'   D_cols  = paste0("D", 1:2),
#'   g_col   = "g_oracle",
#'   mD_cols = paste0("mD", 1:2),
#'   y_col   = "Y",
#'   theta   = demo$truth$theta_true,
#'   j       = 1L,
#'   tau     = 0.5,
#'   s       = 5
#' )
#' ctx <- score_setup(sc, demo$data)
#' subs <- cate_demo_subsets(n = nrow(demo$data), K = 3, s = 5, seed = 1)
#' score_eval_batch(sc, subs, demo$data, ctx)
#'
#' @export
score_plqr_component <- function(
  x0, x_cols, D_cols, g_col, mD_cols, y_col,
  theta, j = 1L, tau = 0.5,
  s
) {
  # ---- factory-time checks ----
  .validate_plqr_factory(D_cols, mD_cols, theta, j, tau, s)
  force(x0)
  force(x_cols)
  force(g_col)
  force(y_col)

  new_score(
    name = paste0("plqr_comp_j", as.integer(j), "_purescore"),
    order = function(n, data, ctx) as.integer(s),
    setup = function(data, ...) {
      X <- .dnn_col_mat(x_cols, data = data, arg = "x_cols", out_arg = "X")
      n <- nrow(X)
      D <- .dnn_col_mat(D_cols, data = data, arg = "D_cols", out_arg = "D")
      mD <- .dnn_col_mat(mD_cols, data = data, arg = "mD_cols", out_arg = "mD")
      gx <- .dnn_col_vec(g_col, data = data, n = n, arg = "g_col", ref = "data")
      Y <- .dnn_col_vec(y_col, data = data, n = n, arg = "y_col", ref = "data")

      # data-time checks
      .validate_target_vector(x0, ncol(X), x0_arg = "x0", x_arg = "X")
      if (.has_na_only(X)) .dnn_abort("argument_invalid", "X cannot contain `NA` values.", argument = "X")
      if (.has_nonfinite_non_na(X)) {
        .dnn_abort("argument_invalid", "X cannot contain non-finite values (`Inf`, `-Inf`, or `NaN`).", argument = "X")
      }
      if (.has_na_only(D)) .dnn_abort("argument_invalid", "D cannot contain `NA` values.", argument = "D")
      if (.has_nonfinite_non_na(D)) {
        .dnn_abort("argument_invalid", "D cannot contain non-finite values (`Inf`, `-Inf`, or `NaN`).", argument = "D")
      }
      if (.has_na_only(mD)) .dnn_abort("argument_invalid", "mD cannot contain `NA` values.", argument = "mD")
      if (.has_nonfinite_non_na(mD)) {
        .dnn_abort(
          "argument_invalid",
          "mD cannot contain non-finite values (`Inf`, `-Inf`, or `NaN`).",
          argument = "mD"
        )
      }
      if (.has_na_only(gx)) .dnn_abort("argument_invalid", "g cannot contain `NA` values.", argument = "g")
      if (.has_nonfinite_non_na(gx)) {
        .dnn_abort("argument_invalid", "g cannot contain non-finite values (`Inf`, `-Inf`, or `NaN`).", argument = "g")
      }
      if (.has_na_only(Y)) .dnn_abort("argument_invalid", "Y cannot contain `NA` values.", argument = "Y")
      if (.has_nonfinite_non_na(Y)) .dnn_abort(
        "argument_invalid",
        "Y cannot contain non-finite values (`Inf`, `-Inf`, or `NaN`).",
        argument = "Y"
      )

      if (any(c(nrow(D), nrow(mD)) != n) || length(gx) != n || length(Y) != n) {
        .dnn_abort("argument_invalid", "X, D, mD, g, and Y must share the same n.", argument = "data")
      }
      if (ncol(D) != ncol(mD)) {
        .dnn_abort("argument_invalid", "D and mD must have the same number of columns.", argument = "mD")
      }
      if (length(theta) != ncol(D)) {
        .dnn_abort("argument_invalid",
          sprintf("`theta` length (%d) must equal ncol(D) (%d).", length(theta), ncol(D)),
          argument = "theta"
        )
      }
      if (!(j >= 1L && j <= ncol(D))) {
        .dnn_abort("argument_invalid", sprintf("`j` must be in 1..ncol(D) (=%d).", ncol(D)), argument = "j")
      }

      # ---- compute and stash the score contributions ----
      score <- NULL
      ctx_create <- .rdistnn_compiled_fun("rustat_plqr_ctx_create")
      g_all <- .rdistnn_compiled_fun("rustat_plqr_g_all")
      if (!is.null(ctx_create) && !is.null(g_all)) {
        # CORRECT argument order + include j
        ptr_ctx <- ctx_create(
          Y, # NumericVector
          D, # NumericMatrix
          gx, # NumericVector
          mD, # NumericMatrix
          as.numeric(theta),
          as.numeric(tau),
          as.integer(j)
        )
        score <- as.numeric(g_all(ptr_ctx))
      } else {
        # R fallback
        below <- .plqr_below_threshold(D, as.numeric(theta), gx, Y)
        resid <- as.numeric(tau) - as.numeric(below)
        score <- resid * (D[, j] - mD[, j])
      }

      if (.has_na_only(score)) {
        .dnn_abort("numerical_range",
          "Computed score contains `NA` values; rescale regressors, coefficients and outcome inputs."
        )
      }
      if (.has_nonfinite_non_na(score)) {
        .dnn_abort("numerical_range",
          "Computed score contains non-finite values (`Inf`, `-Inf`, or `NaN`); rescale regressor and nuisance inputs.")
      }

      list(
        score = as.numeric(score),
        n = nrow(X), p = ncol(X),
        s = as.integer(s), j = as.integer(j)
      )
    },
    eval = function(idxs, data, ctx, rng) {
      psi <- ctx$score[idxs]
      if (.has_na_only(psi)) .dnn_abort("numerical_range", "score contains `NA` values.")
      if (.has_nonfinite_non_na(psi)) .dnn_abort(
        "numerical_range",
        "score contains non-finite values (`Inf`, `-Inf`, or `NaN`)."
      )
      .compensated_sum(psi)
    },
    eval_batch = function(subs, data, ctx, rng) {
      psi <- ctx$score
      if (.has_na_only(psi)) .dnn_abort("numerical_range", "score contains `NA` values.")
      if (.has_nonfinite_non_na(psi)) .dnn_abort(
        "numerical_range",
        "score contains non-finite values (`Inf`, `-Inf`, or `NaN`)."
      )
      if (!is.matrix(subs)) subs <- as.matrix(subs)
      if (nrow(subs) == 0L) {
        return(numeric(0L))
      }

      cpp <- .rdistnn_compiled_fun("rustat_subset_row_sums_g")
      if (!is.null(cpp)) {
        return(cpp(subs, psi))
      }
      vapply(seq_len(nrow(subs)), function(i) .compensated_sum_R(psi[subs[i, ]]), numeric(1L))
    },
    randomized = FALSE,
    symmetric = TRUE,
    ptr = NULL
  )
}
