# Test-only reference implementation of the current fast delete-1 DNN
# jackknife loop. It intentionally keeps the deletion loop explicit so later
# optimized paths can be checked against the scalar contract without relying
# on stored numeric constants.

.dnn_reference_as_signal_matrix <- function(signal, n) {
  if (is.null(dim(signal))) {
    return(matrix(as.numeric(signal), ncol = 1L))
  }
  sig <- as.matrix(signal)
  if (nrow(sig) != n) {
    stop("Layer-0 reference signal has incompatible row count.", call. = FALSE)
  }
  storage.mode(sig) <- "double"
  sig
}

.dnn_reference_jackknife_se_from_loo <- function(loo, q, n, m) {
  if (q == 1L) {
    loo_bar <- colMeans(loo)
    centered <- sweep(loo, 2L, loo_bar, FUN = "-")
    var_jk <- ((n - 1) / n) * colSums(centered^2)
    return(matrix(sqrt(pmax(var_jk, 0)), nrow = 1L, ncol = m))
  }

  se <- matrix(NA_real_, nrow = q, ncol = m)
  for (j in seq_len(q)) {
    mj <- matrix(loo[, j, ], nrow = n, ncol = m)
    mj_bar <- colMeans(mj)
    centered <- sweep(mj, 2L, mj_bar, FUN = "-")
    var_jk <- ((n - 1) / n) * colSums(centered^2)
    se[j, ] <- sqrt(pmax(var_jk, 0))
  }
  se
}

.dnn_reference_eval_rank_weight_columns <- function(spec, rank_matrix) {
  n <- nrow(rank_matrix)
  m <- ncol(rank_matrix)
  W <- matrix(NA_real_, nrow = n, ncol = m)

  for (j in seq_len(m)) {
    if (is.function(spec$eval)) {
      w <- spec$eval(rank_matrix[, j], ctx = NULL)
    } else {
      w <- spec$eval_batch(matrix(rank_matrix[, j], ncol = 1L), ctx = NULL)[, 1L]
    }
    if (!is.numeric(w) || length(w) != n || any(!is.finite(w))) {
      stop("Layer-0 reference rank-weight output is invalid.", call. = FALSE)
    }
    W[, j] <- as.numeric(w)
  }

  W
}

.dnn_reference_fast_jackknife_reference <- function(
  x_mat, x0_mat, signal, distance, method, k, s, normalize, ctx,
  precomputed_ranks = NULL
) {
  n <- nrow(x_mat)
  m <- nrow(x0_mat)
  sig <- .dnn_reference_as_signal_matrix(signal, n = n)
  q <- ncol(sig)

  ranks_full <- precomputed_ranks
  if (is.null(ranks_full)) {
    ranks_full <- weights_dnn(
      X = x_mat,
      X0 = x0_mat,
      distance = distance,
      k = k,
      s = s,
      method = method,
      normalize = normalize,
      ctx = ctx,
      return = "ranks"
    )$ranks
  }

  dx_full <- distance_eval_batch(.resolve_distance(distance), x_mat, x0_mat, ctx)
  ranks_full <- as.matrix(ranks_full)
  if (nrow(ranks_full) != n || ncol(ranks_full) != m) {
    stop("Layer-0 reference ranks must be an n x m matrix.", call. = FALSE)
  }
  storage.mode(ranks_full) <- "integer"

  rw <- if (identical(method, "exact")) {
    rank_weight_dnn_exact(k = as.integer(k), s = as.integer(s))
  } else {
    rank_weight_dnn_asymp(k = as.integer(k), s = as.integer(s), normalize = normalize)
  }

  if (q == 1L) {
    loo <- matrix(NA_real_, nrow = n, ncol = m)
  } else {
    loo <- array(NA_real_, dim = c(n, q, m))
  }

  idx_all <- seq_len(n)
  for (i in idx_all) {
    keep <- idx_all != i

    r_i <- matrix(ranks_full[i, ], nrow = n - 1L, ncol = m, byrow = TRUE)
    rank_matrix_i <- ranks_full[keep, , drop = FALSE]
    rank_matrix_i <- rank_matrix_i - (rank_matrix_i > r_i)
    storage.mode(rank_matrix_i) <- "integer"

    W_i <- .dnn_reference_eval_rank_weight_columns(rw, rank_matrix_i)
    W_i <- .block_oracle_average_weights(W_i, dx_full[keep, , drop = FALSE])
    est_i <- crossprod(sig[keep, , drop = FALSE], W_i)
    if (q == 1L) {
      loo[i, ] <- as.numeric(est_i[1L, ])
    } else {
      loo[i, , ] <- est_i
    }
  }

  .dnn_reference_jackknife_se_from_loo(loo = loo, q = q, n = n, m = m)
}
