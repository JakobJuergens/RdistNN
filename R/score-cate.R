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

  new_score(
    name  = paste0("cate_rho_", localizer),
    order = function(n, data, ctx) ctx$s,        # engine may override s
    setup = function(data, ...) {
      stopifnot(is.null(rho_col) != is.null(rho))  # exactly one provided
      X <- as.matrix(data[, x_cols, drop = FALSE])
      r <- if (!is.null(rho_col)) data[[rho_col]] else rho
      stopifnot(length(r) == nrow(X))

      # compiled context (if shared lib is present) — no fitting here
      ptr_ctx <- if (exists("rustat_cate_ctx_create")) {
        rustat_cate_ctx_create(X, as.numeric(r), as.numeric(x0),
                               localizer, as.integer(k),
                               as.numeric(if (is.null(h)) -1 else h))
      } else NULL

      list(
        X = X, rho = as.numeric(r), x0 = as.numeric(x0),
        k = as.integer(k), h = h, localizer = localizer,
        s = 32L,                        # sensible default; engine can override
        ptr_ctx = ptr_ctx               # compiled ctx if available
      )
    },
    # R batch fallback (vectorized over rows of subs)
    eval_batch = function(subs, data, ctx, rng) {
      if (!is.null(ctx$ptr_ctx) && exists("rustat_cate_eval_batch")) {
        return(rustat_cate_eval_batch(subs, ctx$ptr_ctx))  # compiled fast-path
      }
      K <- nrow(subs); out <- numeric(K)
      for (kk in seq_len(K)) {
        idxs <- subs[kk, ]
        Xi <- ctx$X[idxs, , drop = FALSE]
        dx <- sqrt(rowSums((Xi - matrix(ctx$x0, nrow(Xi), ncol(Xi), TRUE))^2))
        if (ctx$localizer == "dnn1") {
          w <- rep(0, length(dx)); w[ which.min(dx) ] <- 1
        } else if (ctx$localizer == "dnnk") {
          ord <- order(dx); kk2 <- min(ctx$k, length(dx))
          w <- rep(0, length(dx)); w[ ord[seq_len(kk2)] ] <- 1/kk2
        } else {
          hh <- ctx$h
          if (is.null(hh) || !(hh > 0)) {
            md <- stats::median(dx[dx > 0])
            hh <- ifelse(is.finite(md) && md > 0, md, 1)
          }
          w <- exp(-0.5 * (dx/hh)^2); sw <- sum(w); if (sw > 0) w <- w/sw else w[] <- 1/length(w)
        }
        out[kk] <- sum(w * ctx$rho[idxs])
      }
      out
    },
    # R single-subset fallback (rarely used; for completeness)
    eval = function(idxs, data, ctx, rng) {
      Xi <- ctx$X[idxs, , drop = FALSE]
      dx <- sqrt(rowSums((Xi - matrix(ctx$x0, nrow(Xi), ncol(Xi), TRUE))^2))
      if (ctx$localizer == "dnn1") {
        w <- rep(0, length(dx)); w[ which.min(dx) ] <- 1
      } else if (ctx$localizer == "dnnk") {
        ord <- order(dx); kk2 <- min(ctx$k, length(dx))
        w <- rep(0, length(dx)); w[ ord[seq_len(kk2)] ] <- 1/kk2
      } else {
        hh <- ctx$h
        if (is.null(hh) || !(hh > 0)) {
          md <- stats::median(dx[dx > 0]); hh <- ifelse(is.finite(md) && md > 0, md, 1)
        }
        w <- exp(-0.5 * (dx/hh)^2); sw <- sum(w); if (sw > 0) w <- w/sw else w[] <- 1/length(w)
      }
      sum(w * ctx$rho[idxs])
    },
    # Signal that a compiled fast-path may exist (engine will still check ctx$ptr_ctx)
    ptr = TRUE,
    randomized = FALSE, symmetric = TRUE
  )
}
