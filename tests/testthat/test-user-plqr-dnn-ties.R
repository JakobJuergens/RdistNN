skip_if_not(
  all(vapply(
    c("plqr_dnn", "get_distance", "distance_eval_batch"),
    exists, logical(1),
    mode = "function"
  )),
  "required APIs not available"
)

test_that("deterministic ties in distances lead to deterministic ordering and stable moment", {
  set.seed(43)
  n <- 10L
  p <- 2L
  p_d <- 1L
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

  D <- matrix(rnorm(n * p_d), n, p_d)
  mD <- matrix(rnorm(n * p_d), n, p_d)
  q <- rnorm(n)
  Y <- rnorm(n)
  beta <- rep(0.3, p_d)

  # k = s is intentional here (uniform weights make tie ordering the only
  # source of variation), so the degenerate-localization warning is expected.
  out <- suppressWarnings(plqr_dnn(
    X = X, D = D, q = q, mD = mD, Y = Y,
    beta = beta, tau = 0.5, X0 = X0,
    distance = "l2", k = 3L, s = 3L, method = "exact",
    se = FALSE, return = c("W", "ord")
  ))
  expect_true("ord" %in% names(out))

  l2 <- get_distance("l2")
  DX <- distance_eval_batch(l2, X, X0)
  exp_ord1 <- order(DX[, 1], seq_along(DX[, 1]))
  expect_equal(out$ord[, 1], exp_ord1)

  out2 <- suppressWarnings(plqr_dnn(
    X = X, D = D, q = q, mD = mD, Y = Y,
    beta = beta, tau = 0.5, X0 = X0,
    distance = "l2", k = 3L, s = 3L, method = "exact",
    se = FALSE
  ))
  expect_equal(out$moment, out2$moment, tolerance = 1e-12)
})
