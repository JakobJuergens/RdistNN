# R/dist_l1.R
# -----------------------------------------------------------------------------

# Pure-R fallback (vectorized)
.dist_l1_R <- function(X, x0) {
  X <- as.matrix(X)
  x0 <- as.numeric(x0)
  if (length(x0) != ncol(X)) stop("x0 length must equal ncol(X).")
  rowSums(abs(X - matrix(x0, nrow(X), ncol(X), byrow = TRUE)))
}

#' Fast L1 (Manhattan) distance to a point (C++ with R fallback)
#' @export
dist_l1 <- function(X, x0) {
  .assert_no_na(X, x0)
  if (exists("rustat_dist_l1", mode = "function")) {
    return(rustat_dist_l1(as.matrix(X), as.numeric(x0)))
  }
  .dist_l1_R(X, x0)
}

#' L1 distance spec (single-point + batched)
#' @export
distance_l1 <- function() {
  list(
    name = "l1",
    eval = function(X, x0, ctx = NULL) dist_l1(X, x0),
    eval_batch = function(X, X0, ctx = NULL) {
      X <- as.matrix(X)
      X0 <- as.matrix(X0)
      .assert_no_na(X)
      .assert_no_na(X0)
      if (ncol(X0) != ncol(X)) stop("ncol(X0) must match ncol(X).")
      if (exists("rustat_dist_l1_batch", mode = "function")) {
        return(rustat_dist_l1_batch(X, X0))
      }
      # R fallback
      m <- nrow(X0)
      out <- matrix(NA_real_, nrow(X), m)
      for (j in seq_len(m)) out[, j] <- dist_l1(X, X0[j, ])
      out
    }
  )
}

#' Order rows by L1 distance to a point (deterministic ties)
#' @export
order_by_l1 <- function(X, x0) {
  dx <- dist_l1(X, x0)
  ord <- order(dx, seq_along(dx)) # stable tie-break on row index
  list(dx = dx, ord = as.integer(ord))
}

#' Scaled L1 distance (C++ fast path if present, else reuse unscaled after scaling)
#' @export
dist_l1_scaled <- function(X, x0, inv_s) {
  .assert_no_na(X, x0)
  X <- as.matrix(X)
  x0 <- as.numeric(x0)
  inv_s <- as.numeric(inv_s)
  if (length(x0) != ncol(X) || length(inv_s) != ncol(X)) {
    stop("x0 and inv_s must be length p == ncol(X).")
  }

  if (exists("rustat_dist_l1_scaled", mode = "function")) {
    return(rustat_dist_l1_scaled(X, x0, inv_s))
  }
  Xs <- .scale_cols(X, inv_s)
  x0s <- x0 * inv_s
  dist_l1(Xs, x0s)
}

#' Scaled L1 spec: compute inv_s once; reuse unscaled eval/batch after scaling
#' @export
distance_l1_scaled <- function(rule = c("sd", "mad", "iqr", "range", "custom"),
                               scale_vec = NULL, eps = 1e-12) {
  rule <- match.arg(rule)
  base <- distance_l1() # unscaled spec

  list(
    name = paste0("l1_scaled[", rule, "]"),
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
      if (exists("rustat_dist_l1_scaled", mode = "function")) {
        return(dist_l1_scaled(X, x0, ctx$inv_s %||% .compute_inv_scale(X, rule, scale_vec, eps)))
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
      # If no context provided, build a minimal one (inv_s; Xs only if needed)
      if (is.null(ctx)) {
        inv_s <- .compute_inv_scale(X, rule, scale_vec, eps)
        ctx <- list(inv_s = inv_s) # Xs only needed for the unscaled fallback
      }
      # Prefer dedicated C++ batch
      if (exists("rustat_dist_l1_scaled_batch", mode = "function")) {
        return(rustat_dist_l1_scaled_batch(as.matrix(X), as.matrix(X0), as.numeric(ctx$inv_s)))
      }
      # Fallback: rescale & reuse unscaled batch
      Xs <- if (!is.null(ctx$Xs)) ctx$Xs else .scale_cols(as.matrix(X), ctx$inv_s)
      X0s <- .scale_cols(as.matrix(X0), ctx$inv_s)
      base$eval_batch(Xs, X0s, NULL)
    }
  )
}

#' Order rows by scaled L1 (deterministic ties) via unscaled orderer after scaling
#' @export
order_by_l1_scaled <- function(X, x0, rule = c("sd", "mad", "iqr", "range", "custom"),
                               scale_vec = NULL, eps = 1e-12) {
  dspec <- distance_l1_scaled(rule = rule, scale_vec = scale_vec, eps = eps)
  ctx <- dspec$setup(X)

  # Prefer dedicated C++ single if present
  if (exists("rustat_dist_l1_scaled", mode = "function")) {
    dx <- dist_l1_scaled(X, x0, ctx$inv_s)
    ord <- order(dx, seq_along(dx))
    return(list(dx = dx, ord = as.integer(ord)))
  }

  # Otherwise: reuse unscaled orderer on pre-scaled inputs
  order_by_l1(ctx$Xs, x0 * ctx$inv_s)
}
