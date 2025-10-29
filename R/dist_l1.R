# R/dist_l1.R

# Pure-R fallback (vectorized)
.dist_l1_R <- function(X, x0) {
  X  <- as.matrix(X)
  x0 <- as.numeric(x0)
  if (length(x0) != ncol(X)) stop("x0 length must equal ncol(X).")
  rowSums(abs(X - matrix(x0, nrow(X), ncol(X), byrow = TRUE)))
}

#' Fast L1 (Manhattan) distance to a point (C++ with R fallback)
#' @param X numeric matrix (n x p)
#' @param x0 numeric length-p point
#' @return numeric vector of length n with L1 distances
#' @export
dist_l1 <- function(X, x0) {
  if (exists("rustat_dist_l1", mode = "function")) {
    return(rustat_dist_l1(as.matrix(X), as.numeric(x0)))
  }
  .dist_l1_R(X, x0)
}

#' Order rows by L1 distance (C++ with R fallback)
#' @param X numeric matrix (n x p)
#' @param x0 numeric length-p point
#' @return list(dx = numeric n, ord = integer n) with deterministic tie-breaking
#' @export
order_by_l1 <- function(X, x0) {
  X  <- as.matrix(X); x0 <- as.numeric(x0)
  if (exists("rustat_order_by_l1", mode = "function")) {
    return(rustat_order_by_l1(X, x0))
  }
  dx  <- .dist_l1_R(X, x0)
  ord <- order(dx, seq_along(dx))  # tie-break by index
  list(dx = dx, ord = as.integer(ord))
}

# ----------------- dnn_distance integration -----------------

#' L1 distance spec (plug-in) for the DNN engine
#'
#' Returns a `dnn_distance` that computes row-wise L1 (Manhattan)
#' distances from `X` to a target `x0`. Also provides a simple
#' R `eval_batch` for multiple targets `X0`.
#'
#' @return a `dnn_distance` object
#' @export
distance_l1 <- function() {
  new_distance(
    name = "l1",
    eval = function(X, x0, ctx = NULL) {
      dist_l1(X, x0)  # compiled if available; R otherwise
    },
    eval_batch = function(X, X0, ctx = NULL) {
      X  <- as.matrix(X); X0 <- as.matrix(X0)
      if (ncol(X0) != ncol(X)) stop("ncol(X0) must match ncol(X).")
      m <- nrow(X0)
      out <- matrix(NA_real_, nrow(X), m)
      for (j in seq_len(m)) out[, j] <- dist_l1(X, X0[j, ])
      out
    }
  )
}
