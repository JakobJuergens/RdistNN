# Numerical references use the union of original observations with nonzero
# weight in any deletion. Adjacent ranks matter because deletion shifts ranks.
.jackknife_active_ranks <- function(table) {
  c(table != 0, FALSE) | c(FALSE, table != 0)
}

.jackknife_signal_level <- function(signal) {
  if (!nrow(signal)) return(rep(0, ncol(signal)))
  lo <- apply(signal, 2L, min)
  hi <- apply(signal, 2L, max)
  ifelse(lo > 0, lo, ifelse(hi < 0, hi, 0))
}

# Prepare one target at a time, retaining O(n*q) working storage. The same
# prepared inputs enter C++, R prefix sums or explicit frozen-rank deletions.
# Offsets are retained only to reconstruct requested raw deletion estimates.
.dnn_loo_centered <- function(sig, ranks, tab_jk, tab_full = NULL, rw = NULL) {
  ranks <- .rank_check_mat(ranks)
  n <- nrow(ranks)
  sig <- .dnn_signal_matrix(sig, n = n, arg = "sig")
  q <- ncol(sig)
  m <- ncol(ranks)
  active <- .jackknife_active_ranks(tab_jk)
  loo <- array(0, c(n, q, m))
  offset <- matrix(0, q, m)
  for (j in seq_len(m)) {
    used <- active[ranks[, j]]
    level <- .jackknife_signal_level(sig[used, , drop = FALSE])
    z <- matrix(0, n, q)
    z[used, ] <- sweep(sig[used, , drop = FALSE], 2L, level, "-")
    r <- ranks[, j, drop = FALSE]
    value <- if (is.null(rw)) {
      .dnn_loo_prefix(z, r, tab_jk, tab_full)
    } else {
      .dnn_loo_loop(z, r, rw)
    }
    loo[, , j] <- value
    offset[, j] <- level * sum(tab_jk)
  }
  if (q == 1L) loo <- matrix(loo, n, m)
  list(loo = loo, offset = offset)
}

.dnn_restore_loo <- function(loo, offset) {
  n <- dim(loo)[1L]
  q <- nrow(offset)
  m <- ncol(offset)
  if (q == 1L) return(sweep(loo, 2L, offset[1L, ], "+"))
  for (j in seq_len(m)) {
    loo[, , j] <- sweep(matrix(loo[, , j], n, q), 2L, offset[, j], "+")
  }
  loo
}

# Reference paths may recompute sample-dependent geometry. Find their actual
# union of support first; retain no O(n^2) collection of deletion weights.
.jackknife_recomputed_centered <- function(sig, weight_fun, restore = FALSE) {
  n <- nrow(sig)
  q <- ncol(sig)
  first <- weight_fun(1L)
  used <- matrix(FALSE, n, ncol(first))
  m <- ncol(used)
  for (i in seq_len(n)) {
    keep <- seq_len(n) != i
    W <- if (i == 1L) first else weight_fun(i)
    used[keep, ] <- used[keep, ] | (W != 0)
  }
  levels <- vapply(seq_len(m), function(j) {
    .jackknife_signal_level(sig[used[, j], , drop = FALSE])
  }, numeric(q))
  levels <- matrix(levels, q, m)
  loo <- array(0, c(n, q, m))
  for (i in seq_len(n)) {
    keep <- seq_len(n) != i
    W <- weight_fun(i)
    z <- sig[keep, , drop = FALSE]
    for (j in seq_len(m)) {
      active <- W[, j] != 0
      loo[i, , j] <- crossprod(W[active, j],
        sweep(z[active, , drop = FALSE], 2L, levels[, j], "-"))
    }
  }
  if (q == 1L) loo <- matrix(loo, n, m)
  if (restore) return(.dnn_restore_loo(loo, sweep(levels, 2L, colSums(first), "*")))
  loo
}
