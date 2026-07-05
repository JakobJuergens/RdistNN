skip_if_not(
  all(vapply(
    c("np_dnn", "cate_dnn", "plqr_dnn", "weights_dnn"),
    exists, logical(1),
    mode = "function"
  )),
  "DNN wrappers or required helpers not found; skipping file"
)

.make_layer0_boundary_data <- function(n = 18L, p = 2L, m = 3L) {
  stopifnot(n >= 8L, p == 2L, m >= 2L)
  base <- matrix(
    c(
      0, 0,
      0, 0,
      1, 0,
      1, 0,
      0, 1,
      0, 1,
      2, 0,
      0, 2
    ),
    ncol = 2L,
    byrow = TRUE
  )
  extra <- cbind(
    seq_len(n - nrow(base)) / n,
    rev(seq_len(n - nrow(base))) / n
  )
  X <- rbind(base, extra)
  X0 <- rbind(c(0, 0), c(2, 0))
  if (m > 2L) {
    X0 <- rbind(
      X0,
      matrix(c(0.5, 0.5), nrow = m - 2L, ncol = 2L, byrow = TRUE)
    )
  }
  y <- 1 + 0.7 * X[, 1] - 0.4 * X[, 2] + seq_len(n) / (10 * n)
  list(X = X, X0 = X0, y = y)
}

.expect_layer0_ties_and_rank_boundaries <- function(ranks, dx) {
  expect_equal(sort(unique(as.integer(ranks[, 1L]))), seq_len(nrow(ranks)))
  expect_true(any(ranks == 1L))
  expect_true(any(ranks == nrow(ranks)))
  expect_equal(dx[1L, 1L], dx[2L, 1L], tolerance = 1e-12)
  expect_equal(ranks[1:2, 1L], c(1L, 2L))
}

test_that("layer 0 reference locks np_dnn scalar estimates, SEs, and extras", {
  d <- .make_layer0_boundary_data(n = 20L, m = 3L)

  out <- NULL
  expect_warning(
    out <- np_dnn(
      X = d$X,
      X0 = d$X0,
      Y = d$y,
      method = "exact",
      k = 1L,
      s = 1L,
      se = TRUE,
      jackknife = "fast",
      return = c("W", "ranks", "ord", "dx", "y", "ess", "neighbors")
    ),
    "near-\\)uniform|near-uniform",
    ignore.case = TRUE
  )

  ref_se <- .dnn_layer0_fast_jackknife_reference(
    x_mat = d$X,
    x0_mat = d$X0,
    signal = d$y,
    distance = "l2",
    method = "exact",
    k = 1L,
    s = 1L,
    normalize = TRUE,
    ctx = NULL,
    precomputed_ranks = out$ranks
  )

  expect_s3_class(out, "np_result")
  expect_equal(out$mu, as.numeric(crossprod(d$y, out$W)), tolerance = 1e-12)
  expect_equal(out$se, as.numeric(ref_se[1L, ]), tolerance = 1e-12)
  expect_equal(dim(out$W), c(nrow(d$X), nrow(d$X0)))
  expect_equal(dim(out$ranks), dim(out$W))
  expect_equal(dim(out$ord), dim(out$W))
  expect_equal(dim(out$dx), dim(out$W))
  expect_equal(out$y, d$y, tolerance = 1e-12)
  expect_equal(out$ess, 1 / colSums(out$W^2), tolerance = 1e-12)
  expect_equal(dim(out$neighbors), c(1L, nrow(d$X0)))
  .expect_layer0_ties_and_rank_boundaries(out$ranks, out$dx)
})

test_that("layer 0 reference locks cate_dnn scalar estimates, SEs, and extras", {
  d <- .make_layer0_boundary_data(n = 22L, m = 3L)
  psi <- d$y - mean(d$y)

  out <- cate_dnn(
    X = d$X,
    X0 = d$X0,
    pseudo_outcome = psi,
    method = "asymp",
    k = 2L,
    s = 4L,
    se = TRUE,
    jackknife = "fast",
    return = c("W", "ranks", "ord", "dx", "psi", "ess", "neighbors")
  )

  ref_se <- .dnn_layer0_fast_jackknife_reference(
    x_mat = d$X,
    x0_mat = d$X0,
    signal = psi,
    distance = "l2",
    method = "asymp",
    k = 2L,
    s = 4L,
    normalize = TRUE,
    ctx = NULL,
    precomputed_ranks = out$ranks
  )

  expect_s3_class(out, "cate_result")
  expect_equal(out$tau, as.numeric(crossprod(psi, out$W)), tolerance = 1e-12)
  expect_equal(out$se, as.numeric(ref_se[1L, ]), tolerance = 1e-12)
  expect_equal(dim(out$W), c(nrow(d$X), nrow(d$X0)))
  expect_equal(dim(out$ranks), dim(out$W))
  expect_equal(dim(out$ord), dim(out$W))
  expect_equal(dim(out$dx), dim(out$W))
  expect_equal(out$psi, psi, tolerance = 1e-12)
  expect_equal(out$ess, 1 / colSums(out$W^2), tolerance = 1e-12)
  expect_equal(dim(out$neighbors), c(2L, nrow(d$X0)))
  .expect_layer0_ties_and_rank_boundaries(out$ranks, out$dx)
})

test_that("layer 0 reference locks plqr_dnn scalar moment, SE matrix, and extras", {
  set.seed(2105)
  n <- 28L
  p <- 2L
  p_d <- 2L
  m <- 3L
  d <- .make_layer0_boundary_data(n = n, p = p, m = m)
  D <- cbind(1 + d$X[, 1], -0.5 + d$X[, 2])
  mD <- cbind(rep(0.75, n), rep(-0.25, n))
  beta <- c(0.4, -0.2)
  q <- 0.1 + 0.2 * d$X[, 1]
  Y <- as.numeric(D %*% beta + q + seq_len(n) / n > median(D %*% beta + q))

  out <- plqr_dnn(
    X = d$X,
    D = D,
    q = q,
    mD = mD,
    Y = Y,
    beta = beta,
    tau = 0.5,
    X0 = d$X0,
    method = "exact",
    k = 2L,
    s = 5L,
    se = TRUE,
    jackknife = "fast",
    return = c("W", "ranks", "ord", "dx", "g", "ess", "neighbors")
  )

  g_mat <- RdistNN:::.plqr_dnn_g_matrix(D, mD, q, Y, beta = beta, tau = 0.5)
  ref_se <- .dnn_layer0_fast_jackknife_reference(
    x_mat = d$X,
    x0_mat = d$X0,
    signal = g_mat,
    distance = "l2",
    method = "exact",
    k = 2L,
    s = 5L,
    normalize = TRUE,
    ctx = NULL,
    precomputed_ranks = out$ranks
  )

  expect_s3_class(out, "plqr_result")
  expect_equal(unname(out$moment), unname(crossprod(g_mat, out$W)), tolerance = 1e-12)
  expect_equal(unname(out$se), unname(ref_se), tolerance = 1e-12)
  expect_equal(dim(out$moment), c(p_d, m))
  expect_equal(dim(out$se), c(p_d, m))
  expect_equal(dim(out$W), c(nrow(d$X), nrow(d$X0)))
  expect_equal(dim(out$ranks), dim(out$W))
  expect_equal(dim(out$ord), dim(out$W))
  expect_equal(dim(out$dx), dim(out$W))
  expect_equal(out$g, g_mat, tolerance = 1e-12)
  expect_equal(out$ess, 1 / colSums(out$W^2), tolerance = 1e-12)
  expect_equal(dim(out$neighbors), c(2L, nrow(d$X0)))
  .expect_layer0_ties_and_rank_boundaries(out$ranks, out$dx)
})

test_that("layer 0 important scalar error paths remain stable across wrappers", {
  d <- .make_layer0_boundary_data(n = 12L, m = 2L)
  D <- cbind(1 + d$X[, 1], -0.5 + d$X[, 2])
  mD <- matrix(0, nrow = nrow(d$X), ncol = 2L)
  beta <- c(0.4, -0.2)
  q <- rep(0, nrow(d$X))
  Y <- d$y

  expect_error(
    np_dnn(X = d$X, X0 = d$X0, Y = d$y, k = 2L, s = 1L, se = FALSE),
    "k <= s"
  )
  expect_error(
    cate_dnn(X = d$X, X0 = d$X0, pseudo_outcome = d$y, k = 1L, s = nrow(d$X), se = TRUE),
    "s < n"
  )
  expect_error(
    plqr_dnn(
      X = d$X, D = D, q = q, mD = mD, Y = Y,
      beta = beta, X0 = d$X0, x0 = d$X0[1L, ],
      k = 1L, s = 2L, se = FALSE
    ),
    "either `X0` or `x0`",
    ignore.case = TRUE
  )
})
