# tests/testthat/test-user-cate-dnn-custom.R

test_that("custom distance works in cate_dnn; custom rank weights evaluate in batches", {

  set.seed(123)

  nm <- paste0("linf_", as.integer(runif(1, 1e6, 2e6)))
  register_distance(nm, function(X, x0, ctx) {
    apply(abs(X - matrix(x0, nrow(X), ncol(X), TRUE)), 1L, max)
  })
  spec <- get_distance(nm)

  n <- 25
  p <- 3
  m <- 4
  X <- matrix(rnorm(n * p), n, p)
  X0 <- matrix(rnorm(m * p), m, p)
  psi <- rnorm(n)

  df <- data.frame(
    Y = rnorm(n),
    T = rbinom(n, 1, 0.5),
    X1 = X[, 1], X2 = X[, 2], X3 = X[, 3]
  )
  x_cols <- c("X1", "X2", "X3")

  wt <- weights_dnn(
    X = X, X0 = X0, distance = nm,
    k = 3L, s = 4L, method = "exact",
    return = "W"
  )
  tau_manual <- as.numeric(crossprod(psi, wt$W))

  out1 <- cate_dnn(
    data = df, x_cols = x_cols, y_col = "Y", t_col = "T",
    pseudo_outcome = psi, X0 = X0,
    distance = nm, k = 3L, s = 4L, method = "exact",
    se = FALSE, return = "W"
  )
  expect_equal(out1$tau, tau_manual, tolerance = 1e-12)

  DX <- distance_eval_batch(spec, X, X0)
  ORD <- apply(DX, 2L, function(col) order(col, seq_along(col)))
  RKS <- apply(ORD, 2L, function(ord) {
    ranks <- integer(length(ord))
    ranks[ord] <- seq_along(ord)
    ranks
  })

  rw_power <- function(alpha) {
    as_rank_weight(function(r, ctx = NULL) {
      w <- 1 / (as.numeric(r)^alpha)
      w / sum(w)
    }, name = paste0("inv_rank^", alpha))
  }
  rw1 <- rw_power(1)
  rw2 <- rw_power(2)

  W1 <- rank_weight_eval_batch(rw1, RKS)
  W2 <- rank_weight_eval_batch(rw2, RKS)
  tau_rw1 <- as.numeric(crossprod(psi, W1))
  tau_rw2 <- as.numeric(crossprod(psi, W2))

  expect_true(all(abs(colSums(W1) - 1) < 1e-12))
  expect_true(all(abs(colSums(W2) - 1) < 1e-12))
  expect_false(isTRUE(all.equal(tau_rw1, tau_rw2)))
})
