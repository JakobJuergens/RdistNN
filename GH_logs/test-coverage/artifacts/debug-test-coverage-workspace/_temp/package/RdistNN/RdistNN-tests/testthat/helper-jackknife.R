# Helpers for jackknife tests
set.seed(123)

# Squared Euclidean distances (no external deps)
.sqdist <- function(A, B) {
  A <- as.matrix(A)
  B <- as.matrix(B)
  A2 <- rowSums(A^2)
  B2 <- rowSums(B^2)
  -2 * (A %*% t(B)) + matrix(A2, nrow(A), nrow(B)) + matrix(B2, nrow(A), nrow(B), byrow = TRUE)
}

# Deterministic top-k neighbor indices (1-based), with stable tie-breaking
knn_ranks <- function(X, X0 = NULL, k = 5L, exclude_self = TRUE) {
  if (is.null(X0)) X0 <- X
  D <- .sqdist(X0, X)
  Tn <- nrow(D)
  n <- ncol(D)
  ranks <- matrix(NA_integer_, Tn, k)
  for (t in seq_len(Tn)) {
    d <- D[t, ]
    ord <- order(d, seq_len(n)) # stable: break ties by index
    if (exclude_self && identical(X0, X)) ord <- ord[ord != t]
    ranks[t, ] <- ord[seq_len(k)]
  }
  ranks
}
