#' CATE score from precomputed pseudo-outcome rho (no nuisance fitting)
#'
#' @param x0 numeric length-p target point
#' @param x_cols columns of X in `data`
#' @param rho_col optional name of a column in `data` holding rho
#' @param rho optional numeric vector (length n) of rho; used if rho_col is NULL
#' @param localizer "dnn1","dnnk","kernel"
#' @param k integer for dnnk
#' @param h numeric bandwidth for "kernel" (if NULL, auto inside evaluator)
#' @export
score_cate_rho <- function(
    x0, x_cols, rho_col = NULL, rho = NULL,
    localizer = c("dnn1","dnnk","kernel"),
    k = 5L, h = NULL
) {
  localizer <- match.arg(localizer)

  # XOR at construction time
  if (is.null(rho_col) == is.null(rho)) {
    stop("Provide exactly one of 'rho_col' or 'rho' (mutually exclusive).")
  }

  new_score(
    name  = paste0("cate_rho_", localizer),
    order = function(n, data, ctx) ctx$s,

    setup = function(data, ...) {
      X <- as.matrix(data[, x_cols, drop = FALSE])
      r <- if (!is.null(rho_col)) data[[rho_col]] else rho

      # validate dimensions early (these will show as clear errors in tests)
      if (length(x0) != ncol(X)) stop("x0 length must equal ncol(X).")
      if (length(r)  != nrow(X)) stop("rho length must equal nrow(data).")

      # compiled ctx, if available
      ptr_ctx <- if (exists("rustat_cate_ctx_create", mode = "function")) {
        rustat_cate_ctx_create(
          X, as.numeric(r), as.numeric(x0),
          localizer, as.integer(k),
          as.numeric(if (is.null(h)) -1 else h)
        )
      } else NULL

      list(
        X = X, rho = as.numeric(r), x0 = as.numeric(x0),
        k = as.integer(k), h = h, localizer = localizer,
        s = 32L, ptr_ctx = ptr_ctx
      )
    },

    eval_batch = function(subs, data, ctx, rng) {
      # Prefer compiled path if both symbol and ctx are present
      if (!is.null(ctx$ptr_ctx) &&
          exists("rustat_cate_eval_batch", mode = "function")) {
        res <- rustat_cate_eval_batch(subs, ctx$ptr_ctx)
        # be explicit: ensure proper type/length to satisfy validator
        return(as.numeric(res))
      }
      # ---- R fallback (defensive, always finite numeric) ----
      K <- nrow(subs)
      out <- numeric(K)
      for (kk in seq_len(K)) {
        idxs <- subs[kk, ]
        Xi   <- ctx$X[idxs, , drop = FALSE]
        dx   <- sqrt(rowSums((Xi - matrix(ctx$x0, nrow(Xi), ncol(Xi), TRUE))^2))
        w    <- .cate_weights(dx, ctx$localizer, ctx$k, ctx$h)
        rsub <- ctx$rho[idxs]
        # guard against accidental NA in rho
        if (anyNA(rsub) || any(!is.finite(rsub))) {
          stop("rho contains non-finite values in the selected subset.")
        }
        out[kk] <- sum(w * rsub)
      }
      return(out)
    },

    eval = function(idxs, data, ctx, rng) {
      Xi <- ctx$X[idxs, , drop = FALSE]
      dx <- sqrt(rowSums((Xi - matrix(ctx$x0, nrow(Xi), ncol(Xi), TRUE))^2))
      w  <- .cate_weights(dx, ctx$localizer, ctx$k, ctx$h)
      rsub <- ctx$rho[idxs]
      if (anyNA(rsub) || any(!is.finite(rsub))) {
        stop("rho contains non-finite values in the selected subset.")
      }
      val <- sum(w * rsub)
      # explicit scalar numeric return
      return(as.numeric(val))
    },

    ptr = TRUE, randomized = FALSE, symmetric = TRUE
  )
}

# helpers used inside the fallbacks
.cate_weights <- function(dx, localizer, k, h) {
  s <- length(dx)
  if (localizer == "dnn1") {
    w <- numeric(s); w[which.min(dx)] <- 1; return(w)
  }
  if (localizer == "dnnk") {
    ord <- order(dx)                  # stable tie-break by index (R default)
    kk  <- min(if (is.null(k)) 1L else as.integer(k), s)
    w <- numeric(s); w[ord[seq_len(kk)]] <- 1/kk; return(w)
  }
  # kernel: Gaussian, with median(auto) bandwidth; guard all edge cases
  hh <- h
  if (is.null(hh) || !(is.finite(hh) && hh > 0)) {
    pos <- dx[is.finite(dx) & dx > 0]
    if (length(pos) == 0L) {
      hh <- 1.0
    } else {
      pos <- sort(pos)
      m   <- length(pos)
      hh  <- if (m %% 2L == 1L) pos[(m+1L)%/%2L] else 0.5*(pos[m/2L] + pos[m/2L + 1L])
    }
  }
  w <- exp(-0.5 * (dx/hh)^2)
  sw <- sum(w)
  if (!is.finite(sw) || sw <= 0) w[] <- 1/s else w <- w/sw
  w
}

