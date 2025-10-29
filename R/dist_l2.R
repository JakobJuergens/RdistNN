# R/dist_l2.R

# Pure-R fallback (vectorized, stable)
.dist_l2_R <- function(X, x0, squared = TRUE) {
  X  <- as.matrix(X)
  x0 <- as.numeric(x0)
  if (length(x0) != ncol(X)) stop("x0 length must equal ncol(X).")
  dx2 <- rowSums((X - matrix(x0, nrow(X), ncol(X), byrow = TRUE))^2)
  if (squared) dx2 else sqrt(dx2)
}

#' Fast L2 distance to a point (C++ with R fallback)
#' @param X numeric matrix (n x p)
#' @param x0 numeric length-p point
#' @param squared logical; return squared L2 if TRUE
#' @export
dist_l2 <- function(X, x0, squared = TRUE) {
  if (exists("rustat_dist_l2", mode = "function")) {
    return(rustat_dist_l2(as.matrix(X), as.numeric(x0), isTRUE(squared)))
  }
  .dist_l2_R(X, x0, squared)
}

#' Order rows by L2 distance to a point (C++ with R fallback)
#' @param X numeric matrix (n x p)
#' @param x0 numeric length-p point
#' @param squared logical; use squared L2 for ordering if TRUE (same order as L2)
#' @return list(dx = numeric n, ord = integer n) with deterministic tie-breaking
#' @export
order_by_l2 <- function(X, x0, squared = TRUE) {
  X  <- as.matrix(X); x0 <- as.numeric(x0)
  if (exists("rustat_order_by_l2", mode = "function")) {
    return(rustat_order_by_l2(X, x0, isTRUE(squared)))
  }
  dx <- .dist_l2_R(X, x0, squared)
  # Deterministic index tie-break
  ord <- order(dx, seq_along(dx))
  list(dx = dx, ord = as.integer(ord))
}

# ----------------- dnn_distance integration -----------------

#' L2 distance spec (plug-in) for the DNN engine
#'
#' @param squared logical; if TRUE use squared L2 (default). Ordering is identical,
#'   and squared avoids a sqrt, so it is a good default for k-NN and kernels.
#' @return a `dnn_distance` object
#' @export
distance_l2 <- function(squared = TRUE) {
  nm <- "l2"
  new_distance(
    name = nm,
    # no setup needed for plain L2
    eval = function(X, x0, ctx = NULL) {
      dist_l2(X, x0, squared = squared)
    },
    # simple R fast-path for multiple targets; loops over columns of X0
    eval_batch = function(X, X0, ctx = NULL) {
      X  <- as.matrix(X); X0 <- as.matrix(X0)
      if (ncol(X0) != ncol(X)) stop("ncol(X0) must match ncol(X).")
      m <- nrow(X0)
      out <- matrix(NA_real_, nrow(X), m)
      for (j in seq_len(m)) out[, j] <- dist_l2(X, X0[j, ], squared = squared)
      out
    }
  )
}
