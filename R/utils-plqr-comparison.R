# PLQR threshold policy over the shared exact sum-of-products arithmetic.
.plqr_exact_comparison <- function(D, theta, g, y) {
  a <- c(D, g, y)
  b <- c(theta, 1, -1)
  if (any(!is.finite(a)) || any(!is.finite(b))) {
    .dnn_abort("internal_invariant",
      "PLQR threshold comparisons require finite inputs.",
      hint = "Report a reproducible example of this internal inconsistency."
    )
  }
  .exact_sum_products_nonnegative(a, b)
}

# Keep the ordinary matrix-product comparison only when a conservative
# forward-error envelope separates its margin from zero. Range-risk and
# uncertain rows use the exact sign, including inclusive exact equality.
.plqr_below_threshold <- function(D, theta, g, y) {
  threshold <- as.vector(D %*% theta) + g
  bound <- abs(g) + abs(y)
  for (j in seq_len(ncol(D))) {
    product <- D[, j] * theta[j]
    bound <- bound + abs(product)
  }
  L <- ncol(D) + 2
  envelope <- 8 * (L + 1) * .Machine$double.eps * bound + 32 * L * 2^-1074
  margin <- threshold - y
  exceptional <- !is.finite(bound) | bound >= 2^1022 | !is.finite(threshold) |
    !is.finite(margin) | abs(margin) <= envelope
  below <- y <= threshold
  for (i in which(exceptional)) {
    below[i] <- .plqr_exact_comparison(D[i, ], theta, g[i], y[i])
  }
  below
}
