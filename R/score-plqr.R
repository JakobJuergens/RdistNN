#' PLQR diagnostic score: one component of the orthogonal moment (no fitting)
#'
#' For subset S, returns sum_{i in S} w_i(x0) * (tau - 1(Y_i <= D_i^T beta + q(X_i)) * (D_ij - mD_ij)
#' using a localizer around x0. No nuisances are fit inside; all columns must be supplied.
#'
#' @param x0 numeric length-p target point.
#' @param x_cols names/indices of X columns in `data`.
#' @param D_cols names/indices of D columns in `data` (length p_d).
#' @param q_col name of oracle q_tau(X) column in `data`.
#' @param mD_cols names of oracle m_D(X) columns, same length/order as D_cols.
#' @param tau quantile level in (0,1).
#' @param j integer component index in 1..length(D_cols) to extract a scalar score.
#' @param beta numeric length = length(D_cols), fixed coefficients.
#' @param localizer "dnn1","dnnk","kernel".
#' @param k integer for "dnnk".
#' @param h bandwidth for "kernel" (NULL => auto median of positive distances).
#' @export
score_plqr_component <- function(
    x0, x_cols, D_cols, q_col, mD_cols,
    tau = 0.5, j = 1L, beta,
    localizer = c("dnn1","dnnk","kernel"),
    k = 5L, h = NULL
) {
  localizer <- match.arg(localizer)

  # ---- EARLY, EXPLICIT ARG VALIDATION (factory-time) ----------------
  if (!is.numeric(tau) || length(tau) != 1L || !(tau > 0 && tau < 1)) {
    stop("`tau` must be a numeric scalar in (0,1).", call. = FALSE)
  }
  if (!is.numeric(beta)) {
    stop("`beta` must be numeric (length must match length(D_cols)).", call. = FALSE)
  }
  if (length(beta) != length(D_cols)) {
    stop(sprintf("`beta` length (%d) must equal length(D_cols) (%d).",
                 length(beta), length(D_cols)), call. = FALSE)
  }
  if (!is.numeric(j) || length(j) != 1L || j != as.integer(j) || j < 1L || j > length(D_cols)) {
    stop(sprintf("`j` must be an integer in 1..length(D_cols) (%d).",
                 length(D_cols)), call. = FALSE)
  }
  if (length(mD_cols) != length(D_cols)) {
    stop("`mD_cols` must have the same length and order as `D_cols`.", call. = FALSE)
  }

  .plqr_weights <- function(dx, localizer, k, h) {
    s <- length(dx)
    if (localizer == "dnn1") {
      w <- numeric(s); w[which.min(dx)] <- 1; return(w)
    }
    if (localizer == "dnnk") {
      ord <- order(dx)                 # stable tie-break
      kk  <- max(1L, min(as.integer(k), s))
      w <- numeric(s); w[ord[seq_len(kk)]] <- 1/kk; return(w)
    }
    # kernel: Gaussian with exact-median auto bandwidth
    hh <- h
    if (is.null(hh) || !(is.finite(hh) && hh > 0)) {
      pos <- dx[is.finite(dx) & dx > 0]
      if (length(pos) == 0L) hh <- 1.0 else {
        pos <- sort(pos); m <- length(pos)
        hh  <- if (m %% 2L) pos[(m+1L)%/%2L] else 0.5*(pos[m/2L] + pos[m/2L+1L])
      }
    }
    w <- exp(-0.5 * (dx/hh)^2); sw <- sum(w)
    if (!is.finite(sw) || sw <= 0) w[] <- 1/s else w <- w/sw
    w
  }

  new_score(
    name  = paste0("plqr_comp_j", as.integer(j), "_", localizer),
    order = function(n, data, ctx) ctx$s,

    setup = function(data, ...) {
      X   <- as.matrix(data[, x_cols, drop = FALSE])
      D   <- as.matrix(data[, D_cols, drop = FALSE])
      qx  <- as.numeric(data[[q_col]])
      mD  <- as.matrix(data[, mD_cols, drop = FALSE])
      Y   <- as.numeric(data$Y)

      if (length(x0) != ncol(X)) stop("x0 length must equal ncol(X).", call. = FALSE)
      if (!all(is.finite(X)))    stop("X must be finite.",               call. = FALSE)
      if (!all(is.finite(D)))    stop("D must be finite.",               call. = FALSE)
      if (!all(is.finite(qx)))   stop("q_tau(X) must be finite.",        call. = FALSE)
      if (!all(is.finite(mD)))   stop("m_D(X) must be finite.",          call. = FALSE)
      if (!all(is.finite(Y)))    stop("Y must be finite.",               call. = FALSE)

      # ---- CONSISTENCY CHECKS REQUIRING DATA ------------------------
      if (nrow(D) != nrow(X) || length(qx) != nrow(X) ||
          nrow(mD) != nrow(X) || length(Y) != nrow(X)) {
        stop("X, D, Y, q_tau(X), and m_D(X) must have matching n.", call. = FALSE)
      }
      if (ncol(mD) != ncol(D)) {
        stop("m_D(X) must have the same number of columns as D.", call. = FALSE)
      }

      p_d <- ncol(D)
      if (length(beta) != p_d) {
        stop(sprintf("`beta` length (%d) must equal ncol(D) (%d).", length(beta), p_d),
             call. = FALSE)
      }
      if (!(j >= 1L && j <= p_d)) {
        stop(sprintf("`j` must be in 1..ncol(D) (ncol(D)=%d).", p_d), call. = FALSE)
      }

      ptr_ctx <- if (exists("rustat_plqr_ctx_create", mode = "function")) {
        rustat_plqr_ctx_create(
          X, Y, D, qx, mD, as.numeric(x0),
          localizer, as.integer(k), as.numeric(if (is.null(h)) -1 else h),
          as.numeric(tau), as.integer(j), as.numeric(beta)
        )
      } else NULL

      list(
        X = X, D = D, qx = qx, mD = mD, Y = Y,
        x0 = as.numeric(x0), tau = as.numeric(tau),
        j = as.integer(j), beta = as.numeric(beta),
        localizer = localizer, k = as.integer(k), h = h,
        s = 32L, ptr_ctx = ptr_ctx,
        .wfun = .plqr_weights
      )
    },

    eval_batch = function(subs, data, ctx, rng) {
      if (!is.null(ctx$ptr_ctx) && exists("rustat_plqr_eval_batch", mode = "function")) {
        return(as.numeric(rustat_plqr_eval_batch(subs, ctx$ptr_ctx)))
      }
      K <- nrow(subs); out <- numeric(K)
      for (kk in seq_len(K)) {
        idxs <- subs[kk, ]
        Xi <- ctx$X[idxs, , drop = FALSE]
        dx <- sqrt(rowSums((Xi - matrix(ctx$x0, nrow(Xi), ncol(Xi), TRUE))^2))
        w  <- ctx$.wfun(dx, ctx$localizer, ctx$k, ctx$h)

        Dsub  <- ctx$D[idxs, , drop = FALSE]
        rDj   <- Dsub[, ctx$j] - ctx$mD[idxs, ctx$j]
        thr   <- as.vector(Dsub %*% ctx$beta) + ctx$qx[idxs]
        resid <- (ctx$tau - as.numeric(ctx$Y[idxs] <= thr))

        out[kk] <- sum(w * resid * rDj)
      }
      out
    },

    eval = function(idxs, data, ctx, rng) {
      Xi <- ctx$X[idxs, , drop = FALSE]
      dx <- sqrt(rowSums((Xi - matrix(ctx$x0, nrow(Xi), ncol(Xi), TRUE))^2))
      w  <- ctx$.wfun(dx, ctx$localizer, ctx$k, ctx$h)

      Dsub  <- ctx$D[idxs, , drop = FALSE]
      rDj   <- Dsub[, ctx$j] - ctx$mD[idxs, ctx$j]
      thr   <- as.vector(Dsub %*% ctx$beta) + ctx$qx[idxs]
      resid <- (ctx$tau - as.numeric(ctx$Y[idxs] <= thr))

      as.numeric(sum(w * resid * rDj))
    },

    randomized = FALSE, symmetric = TRUE, ptr = TRUE
  )
}
