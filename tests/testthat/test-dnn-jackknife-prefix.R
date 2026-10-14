.make_dnn_prefix_data <- function(n = 26L, m = 4L) {
  if (n < 12L || m < 1L) {
    stop("Layer-2 prefix data requires n >= 12 and m >= 1.", call. = FALSE)
  }

  base <- matrix(
    c(
      0, 0,
      0, 0,
      1, 0,
      1, 0,
      0, 1,
      0, 1,
      2, 0,
      0, 2,
      1, 1,
      1, 1
    ),
    ncol = 2L,
    byrow = TRUE
  )
  extra_n <- n - nrow(base)
  extra <- cbind(seq_len(extra_n) / n, rev(seq_len(extra_n)) / (n + 3L))
  X <- rbind(base, extra)

  targets <- rbind(
    c(0, 0),
    c(2, 0),
    c(0.4, 0.6),
    c(1, 1),
    c(0.9, 0.1)
  )
  X0 <- targets[seq_len(m), , drop = FALSE]
  y <- 0.25 + 0.7 * X[, 1L] - 0.35 * X[, 2L] + seq_len(n) / (9 * n)

  list(X = X, X0 = X0, y = y)
}

.dnn_prefix_signal <- function(y, q = 1L) {
  if (q == 1L) {
    return(matrix(y, ncol = 1L))
  }
  cbind(
    y,
    y^2 - mean(y^2),
    sin(y + seq_along(y) / length(y))
  )[, seq_len(q), drop = FALSE]
}

.dnn_prefix_loop_loo <- function(sig, ranks_full, rw) {
  sig <- as.matrix(sig)
  ranks_full <- as.matrix(ranks_full)
  storage.mode(ranks_full) <- "integer"
  n <- nrow(sig)
  q <- ncol(sig)
  m <- ncol(ranks_full)

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

    W_i <- rank_weight_eval_batch(rw, rank_matrix_i)
    est_i <- crossprod(sig[keep, , drop = FALSE], W_i)
    if (q == 1L) {
      loo[i, ] <- as.numeric(est_i[1L, ])
    } else {
      loo[i, , ] <- est_i
    }
  }

  loo
}

.dnn_prefix_rank_weight <- function(method, k, s, normalize = TRUE) {
  if (identical(method, "exact")) {
    rank_weight_dnn_exact(k = k, s = s)
  } else {
    rank_weight_dnn_asymp(k = k, s = s, normalize = normalize)
  }
}

.expect_dnn_prefix_matches_loop <- function(method, k, s, n, m, q, normalize = TRUE) {
  d <- .make_dnn_prefix_data(n = n, m = m)
  sig <- .dnn_prefix_signal(d$y, q = q)
  ranks <- weights_dnn(
    X = d$X,
    X0 = d$X0,
    method = method,
    k = k,
    s = s,
    normalize = normalize,
    return = "ranks"
  )$ranks
  rw <- .dnn_prefix_rank_weight(method = method, k = k, s = s, normalize = normalize)
  tab_full <- RdistNN:::.rank_weight_table(rw, n)
  tab_jk <- RdistNN:::.rank_weight_table(rw, n - 1L)

  loo_prefix <- RdistNN:::.dnn_loo_prefix(sig, ranks, tab_jk)
  loo_dispatch <- RdistNN:::.dnn_loo_prefix(sig, ranks, tab_jk, tab_full = tab_full)
  loo_loop <- .dnn_prefix_loop_loo(sig, ranks, rw)

  expect_equal(loo_prefix, loo_loop, tolerance = 1e-12)
  expect_equal(loo_dispatch, loo_loop, tolerance = 1e-12)
  expect_equal(
    RdistNN:::.dnn_jackknife_se_from_loo(loo_dispatch, q = q, n = n, m = m),
    RdistNN:::.dnn_jackknife_se_from_loo(loo_loop, q = q, n = n, m = m),
    tolerance = 1e-12
  )
}

test_that("prefix LOO helper matches deletion loop for scalar signals", {
  .expect_dnn_prefix_matches_loop(
    method = "exact", k = 1L, s = 2L, n = 24L, m = 4L, q = 1L
  )
  .expect_dnn_prefix_matches_loop(
    method = "asymp", k = 2L, s = 12L, n = 30L, m = 4L, q = 1L
  )
})

test_that("prefix LOO helper matches deletion loop for vector signals", {
  .expect_dnn_prefix_matches_loop(
    method = "exact", k = 2L, s = 5L, n = 25L, m = 1L, q = 3L
  )
  .expect_dnn_prefix_matches_loop(
    method = "asymp", k = 1L, s = 10L, n = 28L, m = 3L, q = 2L
  )
})

test_that("fast scalar wrappers keep the layer-0 SE contract", {
  d <- .make_dnn_prefix_data(n = 28L, m = 3L)

  np <- np_dnn(
    X = d$X,
    X0 = d$X0,
    Y = d$y,
    method = "exact",
    k = 2L,
    s = 6L,
    se = TRUE,
    jackknife = "fast",
    return = "ranks"
  )
  np_ref <- .dnn_reference_fast_jackknife_reference(
    x_mat = d$X,
    x0_mat = d$X0,
    signal = d$y,
    distance = "l2",
    method = "exact",
    k = 2L,
    s = 6L,
    normalize = TRUE,
    ctx = NULL,
    precomputed_ranks = np$ranks
  )
  expect_equal(np$se, as.numeric(np_ref[1L, ]), tolerance = 1e-12)

  psi <- d$y - mean(d$y)
  cate <- cate_dnn(
    X = d$X,
    X0 = d$X0,
    pseudo_outcome = psi,
    method = "asymp",
    k = 2L,
    s = 8L,
    se = TRUE,
    jackknife = "fast",
    return = "ranks"
  )
  cate_ref <- .dnn_reference_fast_jackknife_reference(
    x_mat = d$X,
    x0_mat = d$X0,
    signal = psi,
    distance = "l2",
    method = "asymp",
    k = 2L,
    s = 8L,
    normalize = TRUE,
    ctx = NULL,
    precomputed_ranks = cate$ranks
  )
  expect_equal(cate$se, as.numeric(cate_ref[1L, ]), tolerance = 1e-12)

  D <- cbind(1 + d$X[, 1L], -0.25 + d$X[, 2L])
  mD <- cbind(rep(0.5, nrow(d$X)), rep(-0.15, nrow(d$X)))
  theta <- c(0.35, -0.1)
  g <- 0.05 + 0.1 * d$X[, 1L]
  Y <- as.numeric(D %*% theta + g + seq_len(nrow(d$X)) / nrow(d$X) > median(D %*% theta + g))
  plqr <- plqr_dnn(
    X = d$X,
    D = D,
    g = g,
    mD = mD,
    Y = Y,
    theta_eval = theta,
    tau = 0.5,
    X0 = d$X0,
    method = "exact",
    k = 2L,
    s = 6L,
    se = TRUE,
    jackknife = "fast",
    return = "ranks"
  )
  score_mat <- RdistNN:::.plqr_dnn_score_matrix(D, mD, g, Y, theta_eval = theta, tau = 0.5)
  plqr_ref <- .dnn_reference_fast_jackknife_reference(
    x_mat = d$X,
    x0_mat = d$X0,
    signal = score_mat,
    distance = "l2",
    method = "exact",
    k = 2L,
    s = 6L,
    normalize = TRUE,
    ctx = NULL,
    precomputed_ranks = plqr$ranks
  )
  expect_equal(unname(plqr$moment_se), unname(plqr_ref), tolerance = 1e-12)
})

test_that("loop fallback remains available for rank weights without tables", {
  d <- .make_dnn_prefix_data(n = 18L, m = 2L)
  sig <- .dnn_prefix_signal(d$y, q = 2L)
  ranks <- weights_dnn(
    X = d$X,
    X0 = d$X0,
    method = "exact",
    k = 1L,
    s = 4L,
    return = "ranks"
  )$ranks
  rw <- as_rank_weight(function(ranks) {
    raw <- (length(ranks) + 1L - as.numeric(ranks))^2
    raw / sum(raw)
  })

  expect_null(RdistNN:::.rank_weight_table(rw, nrow(ranks) - 1L))
  expect_equal(
    RdistNN:::.dnn_loo_loop(sig, ranks, rw),
    .dnn_prefix_loop_loo(sig, ranks, rw),
    tolerance = 1e-12
  )
})
