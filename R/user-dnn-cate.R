# DNN-based CATE wrapper (rank-weighted local DR learner)

#' CATE via DNN rank weights (orthogonalized local averaging)
#'
#' Two equivalent interfaces are supported:
#' **(A) data-frame API** (recommended):
#'   - `data`, `x_cols`, `y_col`, `t_col`, `e_col`, and either `pseudo_outcome`
#'     **or** outcome model columns to build the pseudo-outcome internally:
#'     - `m0_col` + `m1_col` for the full AIPW pseudo-outcome (preferred), or
#'     - `m_col` alone for the Robinson-form pseudo-outcome.
#' **(B) matrix API** (back-compat):
#'   - `X`, `X0`/`x0`, and `pseudo_outcome`.
#'
#' @param data data.frame with at least X columns (and optionally Y,T,m,e).
#' @param x_cols names/indices of covariate columns in `data` (data-frame API).
#' @param y_col,t_col,e_col column names/indices (data-frame API).
#' @param m0_col,m1_col column names for E[Y(0)|X] and E[Y(1)|X] (preferred;
#'   enables the full doubly-robust AIPW pseudo-outcome from the paper).
#' @param m_col column name for a single outcome model (Robinson form; alternative
#'   to `m0_col`/`m1_col` when only E[Y|X] is available).
#' @param pseudo_outcome numeric length n; if provided, Y/T/m/e are not required.
#' @param X numeric n x p (matrix API).
#' @param X0 numeric m x p matrix of targets; mutually exclusive with `x0`.
#' @param x0 numeric length p single target (m = 1).
#' @param distance distance name or spec (e.g. `"l2"`).
#' @param method `"exact"` or `"asymp"`.
#' @param k,s integers for DNN weighting.
#' @param normalize logical; renormalize columns of W to sum 1.
#' @param se logical; compute jackknife SEs (fast path if available).
#' @param jackknife one of `"fast"`, `"loo"`, `"none"`.
#' @param ctx optional context forwarded to distance/ranking.
#' @param return extras from `weights_dnn()`: any of
#'   `c("W","ranks","ord","dx","psi","nuisance","ess","neighbors")`.
#'
#' @return A `cate_result` with fields `tau`, optional `se`, `meta`, and requested extras.
#' @examples
#' set.seed(1)
#' n <- 50
#' x_mat <- matrix(rnorm(n * 2), n, 2)
#' x0_mat <- matrix(rnorm(6), 3, 2)
#' psi <- rnorm(n)
#' out <- cate_dnn(
#'   X = x_mat, X0 = x0_mat, pseudo_outcome = psi,
#'   distance = "l2", method = "exact", k = 2L, s = 3L,
#'   se = FALSE, return = c("W", "tau")
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
    # Placeholder for a future jackknife implementation.
    se_hat <- NULL
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


#' @export
print.cate_dnn_result <- function(x, ...) {
  cat("<cate_dnn_result>\n", sep = "")
  cat("  data  : n=", x$meta$n, ", p=", x$meta$p, ", m=", x$meta$m, "\n", sep = "")
  cat("  dnn   : k=", x$meta$k, ", s=", x$meta$s, ", method=", x$meta$method, "\n", sep = "")
  cat("  rank  : ", x$meta$weight_name, "; distance: ", x$meta$distance_name, "\n", sep = "")
  cat("  components: ", paste(names(x), collapse = ", "), "\n", sep = "")
  invisible(x)
}
