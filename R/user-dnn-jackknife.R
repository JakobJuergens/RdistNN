# Shared deletion orchestration; numeric covariance policy lives in utils-covariance.

.dnn_jackknife_se_from_loo <- function(loo, q, n, m) {
  if (q == 1L) {
    var_jk <- .covariance_matrix(loo, (n - 1) / n, diagonal = TRUE)
    return(matrix(sqrt(pmax(var_jk, 0)), nrow = 1L, ncol = m))
  }

  se <- matrix(NA_real_, nrow = q, ncol = m)
  for (j in seq_len(q)) {
    mj <- matrix(loo[, j, ], nrow = n, ncol = m)
    var_jk <- .covariance_matrix(mj, (n - 1) / n, diagonal = TRUE)
    se[j, ] <- sqrt(pmax(var_jk, 0))
  }
  se
}

.dnn_loo_prefix <- function(sig, ranks, tab_jk, tab_full = NULL) {
  ranks <- .rank_check_mat(ranks)
  n <- nrow(ranks)
  m <- ncol(ranks)
  if (n < 2L) {
    .abort_bad_arg("ranks", "must have at least two rows for delete-1 jackknife")
  }

  sig <- .dnn_signal_matrix(sig, n = n, arg = "sig")
  q <- ncol(sig)

  if (!is.numeric(tab_jk)) {
    .abort_bad_arg("tab_jk", "must be numeric", quoted = FALSE)
  }
  if (length(tab_jk) != n - 1L) {
    .abort_bad_arg(
      "tab_jk",
      "length must equal nrow(ranks) - 1",
      expected = sprintf("length %d", n - 1L),
      actual = sprintf("length %d", length(tab_jk)),
      quoted = FALSE
    )
  }
  .abort_if_na_or_nonfinite(tab_jk, arg = "tab_jk", quoted = FALSE)
  tab_jk <- as.numeric(tab_jk)

  if (is.null(tab_full)) tab_full <- numeric(n)
  if (!is.numeric(tab_full)) .abort_bad_arg("tab_full", "must be numeric", quoted = FALSE)
  if (length(tab_full) != n) {
    .abort_bad_arg("tab_full", "length must equal nrow(ranks)",
      expected = sprintf("length %d", n), actual = sprintf("length %d", length(tab_full)),
      quoted = FALSE)
  }
  .abort_if_na_or_nonfinite(tab_full, arg = "tab_full", quoted = FALSE)
  # One engine owns prefix arithmetic even when the baseline is not requested.
  loo <- jackknife_rank_full(ranks, sig, as.numeric(tab_full), tab_jk)$theta_loo
  if (q == 1L) loo <- matrix(loo, nrow = n, ncol = m)
  loo
}

.dnn_loo_loop <- function(sig, ranks, rw, center = FALSE, groups = NULL) {
  if (!is_rank_weight(rw)) {
    .abort_bad_arg("rw", "must be a rdistnn_rank_weight", quoted = FALSE)
  }
  ranks <- .rank_check_mat(ranks)
  n <- nrow(ranks)
  m <- ncol(ranks)
  if (n < 2L) {
    .abort_bad_arg("ranks", "must have at least two rows for delete-1 jackknife")
  }

  sig <- .dnn_signal_matrix(sig, n = n, arg = "sig")
  q <- ncol(sig)

  idx_all <- seq_len(n)
  weight_fun <- function(i) {
    keep <- (idx_all != i)

    r_i <- matrix(ranks[i, ], nrow = n - 1L, ncol = m, byrow = TRUE)
    rank_matrix_i <- ranks[keep, , drop = FALSE]
    rank_matrix_i <- rank_matrix_i - (rank_matrix_i > r_i)
    if (!is.integer(rank_matrix_i)) storage.mode(rank_matrix_i) <- "integer"

    groups_i <- .dnn_group_subset(groups, keep)
    .dnn_average_distance_weights(rank_weight_eval_batch(rw, rank_matrix_i), groups_i,
      facts = .dnn_compute_group_facts(groups_i))
  }
  .jackknife_recomputed_centered(sig, weight_fun, restore = !center)
}

# Full delete-1 jackknife covariance of a vector estimate per target.
# `loo` uses the same layout as .dnn_jackknife_se_from_loo (n x m for q = 1,
# n x q x m otherwise); the scaling convention ((n-1)/n) matches the
# componentwise SEs, so sqrt(diag(vcov[[j]])) equals the SE matrix column.
.dnn_jackknife_vcov_from_loo <- function(loo, q, n, m, allow_failure = FALSE) {
  vcov <- vector("list", m)
  se <- matrix(NA_real_, nrow = q, ncol = m)
  for (j in seq_len(m)) {
    mj <- if (q == 1L) {
      matrix(loo[, j], nrow = n, ncol = 1L)
    } else {
      matrix(loo[, , j], nrow = n, ncol = q)
    }
    v <- .covariance_matrix(mj, (n - 1) / n, allow_failure = allow_failure)
    vcov[[j]] <- v
    se[, j] <- sqrt(pmax(diag(v), 0))
  }
  list(vcov = vcov, se = se)
}

# Numerically shifted delete-1 estimates for covariance calculations on a
# (possibly multi-component)
# signal at every target: n x m matrix for one component, n x q x m array
# otherwise. Shared by the componentwise-SE and full-covariance routes.
.dnn_jackknife_loo_array <- function(
  x_mat, x0_mat, signal, distance, method, k, s, normalize, ctx,
  signal_arg = "signal",
  jackknife = c("fast", "loo"),
  precomputed_ranks = NULL, precomputed_groups = NULL
) {
  jackknife <- .dnn_match_arg(jackknife)
  n <- nrow(x_mat)
  if (as.integer(s) >= n) {
    .dnn_abort("argument_invalid", sprintf(
      "Delete-1 jackknife requires s < n (got s=%d, n=%d). Use `se = FALSE` or reduce `s`.",
      as.integer(s), n
    ), argument = "s", expected = sprintf("s < n = %d", n), actual = as.integer(s))
  }

  sig <- .dnn_signal_matrix(signal, n = n, arg = signal_arg)
  q <- ncol(sig)
  m <- nrow(x0_mat)

  if (identical(jackknife, "fast")) {
    ranks_full <- precomputed_ranks
    if (is.null(ranks_full) || is.null(precomputed_groups)) {
      ranking_full <- weights_dnn(
        X = x_mat,
        X0 = x0_mat,
        distance = distance,
        k = k,
        s = s,
        method = method,
        normalize = normalize,
        ctx = ctx,
        return = "ranks"
      )
      if (is.null(ranks_full)) ranks_full <- ranking_full$ranks
      precomputed_groups <- ranking_full$groups
    }

    ranks_full <- as.matrix(ranks_full)
    if (nrow(ranks_full) != n || ncol(ranks_full) != m) {
      .dnn_abort("argument_invalid",
        "`precomputed_ranks` must be an n x m ranks matrix.",
        argument = "precomputed_ranks"
      )
    }
    if (!is.integer(ranks_full)) storage.mode(ranks_full) <- "integer"

    rw <- .dnn_rank_weight_spec(method, k, s, normalize)

    tab_jk <- .rank_weight_table(rw, n - 1L)
    loo <- if (is.null(tab_jk)) {
      .dnn_loo_loop(sig, ranks_full, rw, center = TRUE, groups = precomputed_groups)
    } else {
      tab_full <- .rank_weight_table(rw, n)
      .dnn_loo_centered(sig, ranks_full, tab_jk, tab_full = tab_full,
                        groups = precomputed_groups)$loo
    }
    return(loo)
  }

  weight_fun <- function(i) {
    keep <- (seq_len(n) != i)
    weights_dnn(
      X = x_mat[keep, , drop = FALSE],
      X0 = x0_mat,
      distance = distance,
      k = k,
      s = s,
      method = method,
      normalize = normalize,
      ctx = ctx,
      return = "W"
    )$W
  }
  .jackknife_recomputed_centered(sig, weight_fun)
}

.dnn_jackknife_se_loo_matrix <- function(
  x_mat, x0_mat, signal, distance, method, k, s, normalize, ctx,
  signal_arg = "signal",
  jackknife = c("fast", "loo"),
  precomputed_ranks = NULL, precomputed_groups = NULL
) {
  jackknife <- .dnn_match_arg(jackknife)
  n <- nrow(x_mat)
  sig <- .dnn_signal_matrix(signal, n = n, arg = signal_arg)
  loo <- .dnn_jackknife_loo_array(
    x_mat = x_mat,
    x0_mat = x0_mat,
    signal = sig,
    distance = distance,
    method = method,
    k = k,
    s = s,
    normalize = normalize,
    ctx = ctx,
    signal_arg = signal_arg,
    jackknife = jackknife,
    precomputed_ranks = precomputed_ranks,
    precomputed_groups = precomputed_groups
  )
  .dnn_jackknife_se_from_loo(loo = loo, q = ncol(sig), n = n, m = nrow(x0_mat))
}

# Shared rank-weight construction for NP, CATE, paired ratios and grid views.
