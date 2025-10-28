# tests/testthat/test-dist-l2.R

test_that("dist_l2 matches pure R and handles squared correctly", {
  set.seed(1)
  X  <- matrix(rnorm(200), 100, 2)
  x0 <- c(0.1, -0.2)

  dx_r  <- sqrt(rowSums((X - matrix(x0, nrow(X), 2, TRUE))^2))
  dx_r2 <- rowSums((X - matrix(x0, nrow(X), 2, TRUE))^2)

  dx    <- dist_l2(X, x0, squared = FALSE)
  dx2   <- dist_l2(X, x0, squared = TRUE)

  expect_equal(dx,  dx_r,  tolerance = 1e-12)
  expect_equal(dx2, dx_r2, tolerance = 1e-12)
})

test_that("order_by_l2 returns deterministic ties", {
  # Construct equal distances for first two rows
  X  <- rbind(c(0, 1), c(0, -1), matrix(10, 5, 2))
  x0 <- c(0, 0)
  ob <- order_by_l2(X, x0)
  # rows 1 and 2 are both distance 1 -> row 1 must come before row 2
  expect_identical(ob$ord[1:2], c(1L, 2L))
})

test_that("dimension checks", {
  X  <- matrix(0, 10, 3)
  x0 <- 1:2
  expect_error(dist_l2(X, x0), "x0 length must equal ncol\\(X\\)")
  expect_error(order_by_l2(X, x0), "x0 length must equal ncol\\(X\\)")
})
