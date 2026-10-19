# Internal Beta interval masses for validated rank-grid endpoints and shapes.
# Choose the direct tail using the common Beta-mean crossover. Subtracting two
# lower CDF values near one would discard represented positive right-tail mass.
.beta_interval_mass <- function(lower, upper, a, b) {
  n <- max(length(lower), length(upper), length(a), length(b))
  lower <- rep_len(lower, n)
  upper <- rep_len(upper, n)
  a <- rep_len(a, n)
  b <- rep_len(b, n)
  use_lower <- lower <= a / (a + b)
  mass <- numeric(n)
  if (any(use_lower)) {
    idx <- which(use_lower)
    mass[idx] <- pbeta(upper[idx], a[idx], b[idx]) -
      pbeta(lower[idx], a[idx], b[idx])
  }
  if (any(!use_lower)) {
    idx <- which(!use_lower)
    mass[idx] <- pbeta(lower[idx], a[idx], b[idx], lower.tail = FALSE) -
      pbeta(upper[idx], a[idx], b[idx], lower.tail = FALSE)
  }
  mass
}
