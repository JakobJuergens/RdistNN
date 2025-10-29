# R/distance_l2_scaled.R
# -----------------------------------------------------------------------------
# L2 distance with per-feature scaling (z-score/MAD/IQR/range/custom)
# -----------------------------------------------------------------------------

# small helpers (no NA handling by design)
.col_sd  <- function(X) apply(X, 2L, stats::sd)
.col_mad <- function(X) apply(X, 2L, stats::mad, constant = 1.4826)
.col_iqr <- function(X) apply(X, 2L, stats::IQR)
.col_rng <- function(X) apply(X, 2L, function(z) diff(range(z)))

# pure-R core (vectorized) used as fallback
.dist_l2_scaled_R <- function(X, x0, inv_s, squared = TRUE) {
  X  <- as.matrix(X)
  x0 <- as.numeric(x0)
  inv_s <- as.numeric(inv_s)
  if (length(x0) != ncol(X)) stop("x0 length must equal ncol(X).")
  if (length(inv_s) != ncol(X)) stop("inv_s length must equal ncol(X).")
  dif <- sweep(X,  2L, x0,    "-")
  dif <- sweep(dif,2L, inv_s, "*")
  dx2 <- rowSums(dif * dif)
  if (squared) dx2 else sqrt(dx2)
}

#' Scaled L2 distance spec (per-feature standardization)
#'
#' Compute L2 (or squared L2) distances with per-feature scaling. Supported
#' rules: \code{"sd"}, \code{"mad"}, \code{"iqr"}, \code{"range"}, or
#' \code{"custom"} (provide \code{scale_vec} of length \code{p}).
#'
#' @param rule One of \code{c("sd","mad","iqr","range","custom")}.
#' @param scale_vec Numeric length-\code{p} used when \code{rule="custom"}.
#' @param squared Logical; if TRUE return squared L2 (default).
#' @param eps Small positive floor to avoid zero scales in ill-conditioned cols.
#' @return A `dnn_distance` object.
#' @export
distance_l2_scaled <- function(rule = c("sd","mad","iqr","range","custom"),
                               scale_vec = NULL,
                               squared = TRUE,
                               eps = 1e-12) {
  rule <- match.arg(rule)
  nm <- sprintf("l2_scaled[%s]%s", rule, if (squared) "_sq" else "")

  new_distance(
    name = nm,
    setup = function(data, ...) {
      X <- as.matrix(data)
      p <- ncol(X)
      s <- switch(rule,
                  sd    = .col_sd(X),
                  mad   = .col_mad(X),
                  iqr   = .col_iqr(X),
                  range = .col_rng(X),
                  custom = {
                    if (is.null(scale_vec)) stop("Provide `scale_vec` when rule='custom'.")
                    v <- as.numeric(scale_vec)
                    if (length(v) != p) stop("`scale_vec` must have length p = ncol(data).")
                    v
                  }
      )
      s[!is.finite(s) | s < eps] <- eps
      list(inv_s = 1 / s, squared = isTRUE(squared))
    },
    eval = function(X, x0, ctx = NULL) {
      inv_s <- ctx$inv_s
      if (exists("rustat_dist_l2_scaled", mode = "function")) {
        return(rustat_dist_l2_scaled(as.matrix(X), as.numeric(x0), as.numeric(inv_s), isTRUE(ctx$squared)))
      }
      .dist_l2_scaled_R(X, x0, inv_s, squared = isTRUE(ctx$squared))
    },
    eval_batch = function(X, X0, ctx = NULL) {
      inv_s <- ctx$inv_s
      if (exists("rustat_dist_l2_scaled_batch", mode = "function")) {
        return(rustat_dist_l2_scaled_batch(as.matrix(X), as.matrix(X0), as.numeric(inv_s), isTRUE(ctx$squared)))
      }
      # R fallback loop over targets
      X  <- as.matrix(X); X0 <- as.matrix(X0)
      if (ncol(X0) != ncol(X)) stop("ncol(X0) must match ncol(X).")
      m <- nrow(X0)
      out <- matrix(NA_real_, nrow(X), m)
      for (j in seq_len(m)) {
        out[, j] <- .dist_l2_scaled_R(X, X0[j, ], inv_s, squared = isTRUE(ctx$squared))
      }
      out
    }
  )
}
