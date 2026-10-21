# tests/testthat/test-user-cate-dnn-ties.R

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

  # k = s is intentional here (uniform weights make tie ordering the only
  # source of variation), so the degenerate-localization warning is expected.
  out <- suppressWarnings(cate_dnn(
    data = df, x_cols = c("X1", "X2"), y_col = "Y", t_col = "T",
    pseudo_outcome = psi, X0 = X0,
    distance = "l2", k = 3L, s = 3L, method = "exact",
    se = FALSE, return = c("W", "ord")
  ))
  expect_true("ord" %in% names(out))

  l2 <- get_distance("l2")
  DX <- distance_eval_batch(l2, X, X0)
  exp_ord1 <- .tie_reference_order(DX[, 1], X)
  expect_equal(out$ord[, 1], exp_ord1)

  out2 <- suppressWarnings(cate_dnn(
    data = df, x_cols = c("X1", "X2"), y_col = "Y", t_col = "T",
    pseudo_outcome = psi, X0 = X0,
    distance = "l2", k = 3L, s = 3L, method = "exact",
    se = FALSE
  ))
  expect_equal(out$tau, out2$tau, tolerance = 1e-12)
})

test_that("identical covariate rows retain index allocation in nonuniform weights", {
  X <- matrix(c(0, 0, 1, 1), ncol = 1L)
  y <- c(1, 3, 10, 20)
  reference <- function(X, y) {
    subsets <- utils::combn(length(y), 2L)
    mean(vapply(seq_len(ncol(subsets)), function(j) {
      ids <- subsets[, j]
      y[ids[order(X[ids, 1L]^2, ids)[1L]]]
    }, numeric(1L)))
  }
  loo <- vapply(1:4, function(i) reference(X[-i, , drop = FALSE], y[-i]), numeric(1L))
  V <- 3 / 4 * sum((loo - mean(loo))^2)
  f <- cate_dnn(X = X, x0 = 0, pseudo_outcome = y, k = 1, s = 2)
  expect_equal(f$tau, reference(X, y), tolerance = 1e-12)
  expect_equal(unname(vcov(f)), matrix(V), tolerance = 1e-12)
  expect_equal(V, 89 / 12, tolerance = 1e-12)
  perm <- c(2, 1, 3, 4)
  changed <- cate_dnn(X = X[perm, , drop = FALSE], x0 = 0,
    pseudo_outcome = y[perm], k = 1, s = 2)
  expect_equal(changed$tau - f$tau, 1 / 3, tolerance = 1e-12)
})
