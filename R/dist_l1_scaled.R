# R/distance_l1_scaled.R
# -----------------------------------------------------------------------------
# L1 (Manhattan) distance with per-feature scaling (sd/MAD/IQR/range/custom)
# -----------------------------------------------------------------------------

# Reuse small helpers if not already defined (avoid duplicate defs across files)
if (!exists(".col_sd", mode = "function")) {
  .col_sd <- function(X) apply(X, 2L, stats::sd)
}
if (!exists(".col_mad", mode = "function")) {
  .col_mad <- function(X) apply(X, 2L, stats::mad, constant = 1.4826)
}
if (!exists(".col_iqr", mode = "function")) {
  .col_iqr <- function(X) apply(X, 2L, stats::IQR)
}
if (!exists(".col_rng", mode = "function")) {
  .col_rng <- function(X) apply(X, 2L, function(z) diff(range(z)))
}

# Pure-R fallback (vectorized)
.dist_l1_scaled_R <- function(X, x0, inv_s) {
  X <- as.matrix(X)
  x0 <- as.numeric(x0)
  inv_s <- as.numeric(inv_s)
  if (length(x0) != ncol(X)) stop("x0 length must equal ncol(X).")
  if (length(inv_s) != ncol(X)) stop("inv_s length must equal ncol(X).")
  dif <- sweep(X, 2L, x0, "-")
  dif <- sweep(dif, 2L, inv_s, "*")
  rowSums(abs(dif))
}

#' Scaled L1 distance spec (per-feature standardization)
#'
#' Compute L1 (Manhattan) distances with per-feature scaling.
#' Supported rules: \code{"sd"}, \code{"mad"}, \code{"iqr"}, \code{"range"},
#' or \code{"custom"} (provide \code{scale_vec} of length \code{p}).
#'
#' @param rule One of \code{c("sd","mad","iqr","range","custom")}.
#' @param scale_vec Numeric length-\code{p} used when \code{rule="custom"}.
#' @param eps Small positive floor to avoid zero scales in ill-conditioned cols.
#' @return A \code{dnn_distance} object usable with \code{distance_eval} and
#'   \code{distance_eval_batch}.
#' @examples
#' # Example:
#' # dspec <- distance_l1_scaled("mad")
#' # ctx   <- distance_setup(dspec, X)
#' # d     <- distance_eval(dspec, X, x0, ctx)
#' @export
distance_l1_scaled <- function(rule = c("sd", "mad", "iqr", "range", "custom"),
                               scale_vec = NULL,
                               eps = 1e-12) {
  rule <- match.arg(rule)
  nm <- sprintf("l1_scaled[%s]", rule)

  new_distance(
    name = nm,
    setup = function(data, ...) {
      X <- as.matrix(data)
      p <- ncol(X)
      s <- switch(rule,
        sd = .col_sd(X),
        mad = .col_mad(X),
        iqr = .col_iqr(X),
        range = .col_rng(X),
        custom = {
          if (is.null(scale_vec)) stop("Provide `scale_vec` when rule='custom'.")
          v <- as.numeric(scale_vec)
          if (length(v) != p) stop("`scale_vec` must have length p = ncol(data).")
          v
        }
      )
      s[!is.finite(s) | s < eps] <- eps
      list(inv_s = 1 / s)
    },
    eval = function(X, x0, ctx = NULL) {
      inv_s <- ctx$inv_s
      if (exists("rustat_dist_l1_scaled", mode = "function")) {
        return(rustat_dist_l1_scaled(as.matrix(X), as.numeric(x0), as.numeric(inv_s)))
      }
      .dist_l1_scaled_R(X, x0, inv_s)
    },
    eval_batch = function(X, X0, ctx = NULL) {
      inv_s <- ctx$inv_s
      if (exists("rustat_dist_l1_scaled_batch", mode = "function")) {
        return(rustat_dist_l1_scaled_batch(as.matrix(X), as.matrix(X0), as.numeric(inv_s)))
      }
      X <- as.matrix(X)
      X0 <- as.matrix(X0)
      if (ncol(X0) != ncol(X)) stop("ncol(X0) must match ncol(X).")
      m <- nrow(X0)
      out <- matrix(NA_real_, nrow(X), m)
      for (j in seq_len(m)) {
        out[, j] <- .dist_l1_scaled_R(X, X0[j, ], inv_s)
      }
      out
    }
  )
}
