# DNN-based CATE wrapper (rank-weighted local DR learner)

#' CATE via DNN rank weights (orthogonalized local averaging)
#'
#' Two equivalent interfaces are supported:
#' **(A) data-frame API** (recommended):
#'   - `data`, `x_cols`, `y_col`, `t_col`, and either `pseudo_outcome` **or**
#'     both `m_col` and `e_col` to build the pseudo-outcome internally.
#' **(B) matrix API** (back-compat):
#'   - `X`, `X0`/`x0`, and `pseudo_outcome`.
#'
#' @param data data.frame with at least X columns (and optionally Y,T,m,e).
#' @param x_cols names/indices of covariate columns in `data` (data-frame API).
#' @param y_col,t_col,m_col,e_col column names/indices (data-frame API).
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
    stopifnot(is.data.frame(data))
    x_mat <- as.matrix(data[, x_cols, drop = FALSE])
    n <- nrow(x_mat)
    p <- ncol(x_mat)

    if (!is.null(pseudo_outcome)) {
      psi <- as.numeric(pseudo_outcome)
      if (length(psi) != n) stop("`pseudo_outcome` length must equal nrow(X).", call. = FALSE)
      m <- e <- NULL
    } else {
      # Need y,t and nuisance m,e to build ψ
      if (is.null(y_col) || is.null(t_col) || is.null(m_col) || is.null(e_col)) {
        stop("With `data=`, provide `pseudo_outcome` OR {y_col, t_col, m_col, e_col}.", call. = FALSE)
      }
      y_vec <- as.numeric(data[[y_col]])
      t_vec <- as.numeric(data[[t_col]])
      if (any(!t_vec %in% c(0, 1))) stop("`t_col` must be binary in {0,1}.", call. = FALSE)
      m <- as.numeric(data[[m_col]])
      e <- as.numeric(data[[e_col]])
      if (length(y_vec) != n || length(t_vec) != n || length(m) != n || length(e) != n) {
        stop("y_col/t_col/m_col/e_col must all have length n.", call. = FALSE)
      }
      den <- pmax(e * (1 - e), 1e-8)
      psi <- (y_vec - m) * (t_vec - e) / den
    }
  } else {
    # matrix API
    if (is.null(X)) stop("Provide either `data` ... or matrix interface via `X`.", call. = FALSE)
    x_mat <- as.matrix(X)
    n <- nrow(x_mat)
    p <- ncol(x_mat)
    if (is.null(pseudo_outcome)) stop("Matrix API requires `pseudo_outcome`.", call. = FALSE)
    psi <- as.numeric(pseudo_outcome)
    if (length(psi) != n) stop("`pseudo_outcome` length must equal nrow(X).", call. = FALSE)
    m <- e <- NULL
  }

  # Targets
  if (!is.null(X0) && !is.null(x0)) stop("Provide either `X0` or `x0`, not both.", call. = FALSE)
  if (is.null(X0)) {
    if (is.null(x0)) stop("Provide `X0` (m x p) or a single `x0` (length p).", call. = FALSE)
    x0_mat <- matrix(as.numeric(x0), nrow = 1L)
  } else {
    x0_mat <- as.matrix(X0)
  }
  if (ncol(x0_mat) != p) stop("ncol(X0) must equal ncol(X).", call. = FALSE)

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
  if ("nuisance" %in% return && !is.null(m) && !is.null(e)) {
    out$nuisance <- list(m = m, e = e)
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
      paste0(round(head(x$tau, 5), 4), collapse = ", "),
      "  (SEs available)\n",
      sep = ""
    )
  } else {
    cat("  est   : tau[1:min(5,m)] = ",
      paste0(round(head(x$tau, 5), 4), collapse = ", "),
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
