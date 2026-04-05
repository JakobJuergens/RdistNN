# tests/testthat/test-user-cate-dnn-ties.R

skip_if_not(
  all(vapply(
    c("cate_dnn", "get_distance", "distance_eval_batch"),
    exists, logical(1),
    mode = "function"
  )),
  "required APIs not available"
)

test_that("deterministic ties in distances lead to deterministic ordering and stable tau", {
  set.seed(42)
  n <- 10L
  p <- 2L
  m <- 3L
  X <- cbind(seq_len(n), rep(0, n))
  X0 <- matrix(
    c(
      5, 0,
      5, 1,
      5, -1
    ),
    m, p,
    byrow = TRUE
  )

  psi <- rnorm(n)
  df <- data.frame(
    Y = rnorm(n), T = rbinom(n, 1, 0.5),
    X1 = X[, 1], X2 = X[, 2]
  )

  out <- cate_dnn(
    data = df, x_cols = c("X1", "X2"), y_col = "Y", t_col = "T",
    pseudo_outcome = psi, X0 = X0,
    distance = "l2", k = 3L, s = 3L, method = "exact",
    se = FALSE, return = c("W", "ord")
  )
  expect_true("ord" %in% names(out))

  l2 <- get_distance("l2")
  DX <- distance_eval_batch(l2, X, X0)
  exp_ord1 <- order(DX[, 1], seq_along(DX[, 1]))
  expect_equal(out$ord[, 1], exp_ord1)

  out2 <- cate_dnn(
    data = df, x_cols = c("X1", "X2"), y_col = "Y", t_col = "T",
    pseudo_outcome = psi, X0 = X0,
    distance = "l2", k = 3L, s = 3L, method = "exact",
    se = FALSE
  )
  expect_equal(out$tau, out2$tau, tolerance = 1e-12)
})
