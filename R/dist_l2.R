# R/dist_l2.R

# Pure-R fallback (vectorized, stable)
.dist_l2_R <- function(X, x0, squared = FALSE) {
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
dist_l2 <- function(X, x0, squared = FALSE) {
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
order_by_l2 <- function(X, x0, squared = FALSE) {
  X  <- as.matrix(X); x0 <- as.numeric(x0)
  if (exists("rustat_order_by_l2", mode = "function")) {
    return(rustat_order_by_l2(X, x0, isTRUE(squared)))
  }
  dx <- .dist_l2_R(X, x0, squared)
  # Deterministic index tie-break (like C++)
  ord <- order(dx, seq_along(dx))
  list(dx = dx, ord = as.integer(ord))
}
