test_that("dist_l1 matches pure R implementation", {
  set.seed(42)
  X  <- matrix(rnorm(300), 100, 3)
  x0 <- c(0.3, -0.7, 0.2)

  dx_ref <- rowSums(abs(X - matrix(x0, nrow(X), 3, TRUE)))
  dx     <- dist_l1(X, x0)

  expect_equal(dx, dx_ref, tolerance = 1e-12)
})

test_that("order_by_l1 is deterministic on ties", {
  # Rows 1 and 2 at equal L1 distance (1), rest far away
  X  <- rbind(c(0, 1), c(0, -1), matrix(10, 5, 2))
  x0 <- c(0, 0)
  ob <- order_by_l1(X, x0)
  expect_identical(ob$ord[1:2], c(1L, 2L))  # index tie-break
})

test_that("dimension checks for L1", {
  X  <- matrix(0, 10, 3)
  x0 <- 1:2
  expect_error(dist_l1(X, x0), "x0 length must equal ncol\\(X\\)")
  expect_error(order_by_l1(X, x0), "x0 length must equal ncol\\(X\\)")
})
