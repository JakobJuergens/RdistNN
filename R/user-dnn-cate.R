# R/cate-dnn.R
# -----------------------------------------------------------------------------
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
#' @param return extras from `weights_dnn()`: any of `c("W","ranks","ord","dx")`,
#'   plus `"psi"`, `"nuisance"`, `"ess"`, `"neighbors"`, `"loo"`.
#'
#' @return A `cate_result` with fields `tau`, optional `se`, `meta`, and requested extras.
#' @examples
#' set.seed(1)
#' n <- 50
#' X <- matrix(rnorm(n * 2), n, 2)
#' X0 <- matrix(rnorm(6), 3, 2)
#' psi <- rnorm(n)
#' out <- cate_dnn(
#'   X = X, X0 = X0, pseudo_outcome = psi,
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

  # ---- Dispatch: data-frame API vs matrix API --------------------------------
  if (!is.null(data)) {
    stopifnot(is.data.frame(data))
    X <- as.matrix(data[, x_cols, drop = FALSE])
    n <- nrow(X)
    p <- ncol(X)

    if (!is.null(pseudo_outcome)) {
      psi <- as.numeric(pseudo_outcome)
      if (length(psi) != n) stop("`pseudo_outcome` length must equal nrow(X).", call. = FALSE)
      m <- e <- NULL
    } else {
      # Need y,t and nuisance m,e to build ψ
      if (is.null(y_col) || is.null(t_col) || is.null(m_col) || is.null(e_col)) {
        stop("With `data=`, provide `pseudo_outcome` OR {y_col, t_col, m_col, e_col}.", call. = FALSE)
      }
      Y <- as.numeric(data[[y_col]])
      Tt <- as.numeric(data[[t_col]])
      if (any(!Tt %in% c(0, 1))) stop("`t_col` must be binary in {0,1}.", call. = FALSE)
      m <- as.numeric(data[[m_col]])
      e <- as.numeric(data[[e_col]])
      if (length(Y) != n || length(Tt) != n || length(m) != n || length(e) != n) {
        stop("y_col/t_col/m_col/e_col must all have length n.", call. = FALSE)
      }
      den <- pmax(e * (1 - e), 1e-8)
      psi <- (Y - m) * (Tt - e) / den
    }
  } else {
    # matrix API
    if (is.null(X)) stop("Provide either `data` ... or matrix interface via `X`.", call. = FALSE)
    X <- as.matrix(X)
    n <- nrow(X)
    p <- ncol(X)
    if (is.null(pseudo_outcome)) stop("Matrix API requires `pseudo_outcome`.", call. = FALSE)
    psi <- as.numeric(pseudo_outcome)
    if (length(psi) != n) stop("`pseudo_outcome` length must equal nrow(X).", call. = FALSE)
    m <- e <- NULL
  }

  # Targets
  if (!is.null(X0) && !is.null(x0)) stop("Provide either `X0` or `x0`, not both.", call. = FALSE)
  if (is.null(X0)) {
    if (is.null(x0)) stop("Provide `X0` (m x p) or a single `x0` (length p).", call. = FALSE)
    X0 <- matrix(as.numeric(x0), nrow = 1L)
  } else {
    X0 <- as.matrix(X0)
  }
  if (ncol(X0) != p) stop("ncol(X0) must equal ncol(X).", call. = FALSE)

  # ---- Weights (rank-only DNN) -----------------------------------------------
  keep_extras <- intersect(return, c("ranks", "ord", "dx"))
  wt <- weights_dnn(
    X = X, X0 = X0,
    distance = distance,
    k = as.integer(k), s = as.integer(s),
    method = method, normalize = normalize, ctx = ctx,
    return = unique(c("W", keep_extras))
  )
  W <- wt$W
  m_out <- ncol(W)

  # ---- Point estimates --------------------------------------------------------
  tau <- as.numeric(crossprod(psi, W)) # length m_out

  # ---- Jackknife SEs ----------------------------------------------------------
  se_hat <- NULL
  # if (isTRUE(se) && m_out > 0L) {
  #   se_try_fast <- function() {
  #     if (!exists("jackknife_fast_rank", mode = "function")) return(NULL)
  #     R <- wt$ranks
  #     if (is.null(R)) return(NULL)
  #     out <- try(jackknife_fast_rank(R, matrix(psi, ncol = 1L)), silent = TRUE)
  #     if (inherits(out, "try-error")) return(NULL)
  #     if (!is.null(out$theta0) && !is.null(out$theta_loo)) {
  #       theta0 <- as.numeric(out$theta0)
  #       loo    <- out$theta_loo
  #       if (is.matrix(loo) && length(theta0) == m_out) {
  #         v <- (n - 1) / n * colSums((loo - matrix(theta0, n, m_out, byrow = TRUE))^2)
  #         return(sqrt(pmax(v, 0)))
  #       }
  #     }
  #     NULL
  #   }
  #
  #   se_hat <- se_try_fast()
  #   if (is.null(se_hat) && jackknife != "none") {
  #     # LOO fallback
  #     loo <- matrix(NA_real_, n, m_out)
  #     for (i in seq_len(n)) {
  #       Xi   <- X[-i, , drop = FALSE]
  #       psii <- psi[-i]
  #       Wi <- weights_dnn(
  #         X = Xi, X0 = X0, distance = distance,
  #         k = as.integer(k), s = as.integer(s),
  #         method = method, normalize = normalize, ctx = ctx,
  #         return = "W"
  #       )
  #       loo[i, ] <- as.numeric(crossprod(psii, Wi$W))
  #     }
  #     v <- (n - 1) / n * colSums((loo - matrix(tau, n, m_out, byrow = TRUE))^2)
  #     se_hat <- sqrt(pmax(v, 0))
  #   }
  # }

  # ---- Assemble ---------------------------------------------------------------
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
  if ("W" %in% return) out$W <- W
  if ("ranks" %in% return && !is.null(wt$ranks)) out$ranks <- wt$ranks
  if ("ord" %in% return && !is.null(wt$ord)) out$ord <- wt$ord
  if ("dx" %in% return && !is.null(wt$dx)) out$dx <- wt$dx
  if ("psi" %in% return) out$psi <- psi
  if ("nuisance" %in% return && !is.null(m) && !is.null(e)) {
    out$nuisance <- list(m = m, e = e)
  }
  if ("ess" %in% return) {
    out$ess <- 1 / colSums(W^2)
  }
  if ("neighbors" %in% return) {
    k_eff <- max(1L, min(nrow(W), as.integer(k)))
    if (!is.null(wt$ord)) {
      out$neighbors <- wt$ord[seq_len(k_eff), , drop = FALSE]
    } else if (!is.null(wt$ranks)) {
      out$neighbors <- apply(wt$ranks, 2L, function(r) order(r)[seq_len(k_eff)])
    }
  }
  if (exists("loo", inherits = FALSE) && "loo" %in% return) out$loo <- loo

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
