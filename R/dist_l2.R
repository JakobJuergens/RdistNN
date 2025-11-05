# R/dist_l2.R
# -----------------------------------------------------------------------------

.dist_l2_R <- function(X, x0, squared = FALSE) {
  X <- as.matrix(X)
  x0 <- as.numeric(x0)
  if (length(x0) != ncol(X)) stop("x0 length must equal ncol(X).")
  Z <- X - matrix(x0, nrow(X), ncol(X), byrow = TRUE)
  if (squared) rowSums(Z * Z) else sqrt(rowSums(Z * Z))
}

#' Fast L2 distance (or squared L2) to a point (C++ with R fallback)
#' @export
dist_l2 <- function(X, x0, squared = FALSE) {
  .assert_no_na(X, x0)
  if (exists("rustat_dist_l2", mode = "function")) {
    return(rustat_dist_l2(as.matrix(X), as.numeric(x0), isTRUE(squared)))
  }
  .dist_l2_R(X, x0, squared)
}

#' L2 distance spec (single-point + batched)
#' @export
distance_l2 <- function(squared = FALSE) {
  list(
    name = if (squared) "l2_sq" else "l2",
    eval = function(X, x0, ctx = NULL) dist_l2(X, x0, squared),
    eval_batch = function(X, X0, ctx = NULL) {
      X <- as.matrix(X)
      X0 <- as.matrix(X0)
      .assert_no_na(X)
      .assert_no_na(X0)
      if (ncol(X0) != ncol(X)) stop("ncol(X0) must match ncol(X).")
      if (exists("rustat_dist_l2_batch", mode = "function")) {
        return(rustat_dist_l2_batch(X, X0, isTRUE(squared)))
      }
      m <- nrow(X0)
      out <- matrix(NA_real_, nrow(X), m)
      for (j in seq_len(m)) out[, j] <- dist_l2(X, X0[j, ], squared)
      out
    }
  )
}

#' Order rows by L2 distance to a point (deterministic ties)
#' @export
order_by_l2 <- function(X, x0, squared = FALSE) {
  dx <- dist_l2(X, x0, squared)
  ord <- order(dx, seq_along(dx))
  list(dx = dx, ord = as.integer(ord))
}

#' Scaled L2 distance (or squared) — C++ fast path if present, else reuse unscaled
#' @export
dist_l2_scaled <- function(X, x0, inv_s, squared = FALSE) {
  .assert_no_na(X, x0)
  X <- as.matrix(X)
  x0 <- as.numeric(x0)
  inv_s <- as.numeric(inv_s)
  if (length(x0) != ncol(X) || length(inv_s) != ncol(X)) {
    stop("x0 and inv_s must be length p == ncol(X).")
  }

  if (exists("rustat_dist_l2_scaled", mode = "function")) {
    return(rustat_dist_l2_scaled(X, x0, inv_s, isTRUE(squared)))
  }
  Xs <- .scale_cols(X, inv_s)
  x0s <- x0 * inv_s
  dist_l2(Xs, x0s, squared = squared)
}

#' Scaled L2 spec (squared optional): reuse unscaled spec after scaling
#' @export
distance_l2_scaled <- function(rule = c("sd", "mad", "iqr", "range", "custom"),
                               scale_vec = NULL, eps = 1e-12,
                               squared = FALSE) {
  rule <- match.arg(rule)
  base <- distance_l2(squared = squared)

  list(
    name = paste0(if (squared) "l2_scaled_sq[" else "l2_scaled[", rule, "]"),
    setup = function(X) {
      .assert_no_na(X)
      inv_s <- .compute_inv_scale(X, rule, scale_vec, eps)
      list(
        inv_s = inv_s,
        Xs = .scale_cols(as.matrix(X), inv_s)
      )
    },
    eval = function(X, x0, ctx) {
      .assert_no_na(X, x0)
      if (exists("rustat_dist_l2_scaled", mode = "function")) {
        return(dist_l2_scaled(X, x0, ctx$inv_s %||% .compute_inv_scale(X, rule, scale_vec, eps), squared))
      }
      if (is.null(ctx)) {
        inv_s <- .compute_inv_scale(X, rule, scale_vec, eps)
        return(base$eval(.scale_cols(as.matrix(X), inv_s), x0 * inv_s, NULL))
      }
      base$eval(ctx$Xs, x0 * ctx$inv_s, NULL)
    },
    eval_batch = function(X, X0, ctx) {
      .assert_no_na(X)
      .assert_no_na(X0)
      if (is.null(ctx)) {
        inv_s <- .compute_inv_scale(X, rule, scale_vec, eps)
        ctx <- list(inv_s = inv_s)
      }
      if (exists("rustat_dist_l2_scaled_batch", mode = "function")) {
        return(rustat_dist_l2_scaled_batch(
          as.matrix(X), as.matrix(X0), as.numeric(ctx$inv_s),
          isTRUE(squared)
        ))
      }
      Xs <- if (!is.null(ctx$Xs)) ctx$Xs else .scale_cols(as.matrix(X), ctx$inv_s)
      X0s <- .scale_cols(as.matrix(X0), ctx$inv_s)
      base$eval_batch(Xs, X0s, NULL)
    }
  )
}

#' Order rows by scaled L2 (deterministic ties) via unscaled orderer after scaling
#' @export
order_by_l2_scaled <- function(X, x0, rule = c("sd", "mad", "iqr", "range", "custom"),
                               scale_vec = NULL, eps = 1e-12, squared = FALSE) {
  dspec <- distance_l2_scaled(
    rule = rule, scale_vec = scale_vec, eps = eps,
    squared = squared
  )
  ctx <- dspec$setup(X)

  if (exists("rustat_dist_l2_scaled", mode = "function")) {
    dx <- dist_l2_scaled(X, x0, ctx$inv_s, squared)
    ord <- order(dx, seq_along(dx))
    return(list(dx = dx, ord = as.integer(ord)))
  }

  order_by_l2(ctx$Xs, x0 * ctx$inv_s, squared = squared)
}
