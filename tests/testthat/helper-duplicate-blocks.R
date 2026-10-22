# Independent explicit-deletion oracle: exact equality of original row values,
# ordinary canonical distance ordering and direct averaged rank-table weights.
.with_duplicate_policy_warnings <- function(code) {
  # Numerical and retained-state fixtures intentionally contain exact-X ties.
  # Dedicated warning tests check the public message and its call scope.
  withCallingHandlers(code, warning = function(w) {
    if (startsWith(conditionMessage(w), "Identical covariate rows detected;")) {
      invokeRestart("muffleWarning")
    }
  })
}

.block_oracle_average_weights <- function(W, X) {
  result <- W
  for (i in seq_len(nrow(X))) {
    same <- vapply(seq_len(nrow(X)), function(h) all(X[h, ] == X[i, ]), logical(1L))
    result[i, ] <- colMeans(W[same, , drop = FALSE])
  }
  result
}

.block_oracle_weights <- function(X, X0, table, ranks = NULL) {
  n <- nrow(X)
  out <- matrix(0, n, nrow(X0))
  for (j in seq_len(nrow(X0))) {
    ord <- if (is.null(ranks)) {
      dx <- sqrt(rowSums(sweep(X, 2L, X0[j, ], "-")^2))
      do.call(order, c(list(dx), unname(as.list(as.data.frame(X))), list(seq_len(n))))
    } else order(ranks[, j])
    raw <- numeric(n)
    raw[ord] <- table
    for (i in seq_len(n)) {
      same <- vapply(seq_len(n), function(h) all(X[h, ] == X[i, ]), logical(1L))
      out[i, j] <- mean(raw[same])
    }
  }
  out
}

.block_exact_table <- function(n, k, s) {
  vapply(seq_len(n), function(r) {
    sum(vapply(seq_len(k), function(j) {
      choose(r - 1, j - 1) * choose(n - r, s - j)
    }, numeric(1L))) / (k * choose(n, s))
  }, numeric(1L))
}

.block_oracle_loo <- function(X, X0, signal, table, ranks = NULL) {
  signal <- as.matrix(signal)
  n <- nrow(X)
  result <- array(0, c(n, ncol(signal), nrow(X0)))
  for (i in seq_len(n)) {
    deleted <- if (is.null(ranks)) NULL else {
      r <- ranks[-i, , drop = FALSE]
      r - (r > matrix(rep(ranks[i, ], each = n - 1L), n - 1L))
    }
    W <- .block_oracle_weights(X[-i, , drop = FALSE], X0, table, deleted)
    result[i, , ] <- crossprod(signal[-i, , drop = FALSE], W)
  }
  result
}
