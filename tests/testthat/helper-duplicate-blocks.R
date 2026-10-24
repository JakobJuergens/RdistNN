# Independent explicit-deletion oracle: exact equality of computed distances,
# ordinary canonical distance ordering and direct averaged rank-table weights.
.with_test_warning_families <- function(code, families) {
  classes <- paste0("rdistnn_", families, "_warning")
  withCallingHandlers(code, warning = function(w) {
    if (inherits(w, classes)) invokeRestart("muffleWarning")
  })
}

.with_localization_policy_warnings <- function(code, near_uniform = TRUE) {
  # Only owned exact-distance and near-uniform condition classes are muffled.
  # Numerical, optimization, rescaling and other warnings remain visible.
  # Dedicated warning tests check the public message and its call scope.
  .with_test_warning_families(code, c("distance_ties",
    if (near_uniform) "weak_localization"))
}

.with_distance_tie_warning <- function(code) {
  .with_localization_policy_warnings(code, near_uniform = FALSE)
}

# Status-bearing fixtures deliberately contain unavailable outputs. Their
# numerical/status suites accept this one owned notice; multiplicity, contents
# and warning promotion are tested by test-reporting-recovery.R.
.with_partial_result_warning <- function(code) {
  .with_test_warning_families(code, "partial_result")
}

.block_oracle_average_weights <- function(W, dx) {
  result <- W
  for (j in seq_len(ncol(W))) for (i in seq_len(nrow(W))) {
    same <- dx[, j] == dx[i, j]
    result[i, j] <- mean(W[same, j])
  }
  result
}

.block_oracle_weights <- function(X, X0, table, ranks = NULL) {
  n <- nrow(X)
  dx <- vapply(seq_len(nrow(X0)), function(j) {
    sqrt(rowSums(sweep(X, 2L, X0[j, ], "-")^2))
  }, numeric(n))
  raw <- vapply(seq_len(nrow(X0)), function(j) {
    ord <- if (is.null(ranks)) order(dx[, j], seq_len(n)) else order(ranks[, j])
    w <- numeric(n)
    w[ord] <- table
    w
  }, numeric(n))
  .block_oracle_average_weights(raw, dx)
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
