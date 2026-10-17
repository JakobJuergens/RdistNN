test_that("CLATE matrix data-frame and tibble interfaces agree", {
  a <- .clate_fit_arguments()
  matrix_fit <- do.call(clate_dnn, a)
  d <- .clate_test_frame(a)
  args <- .clate_frame_arguments(a, d)
  frame_fit <- do.call(clate_dnn, args)
  expect_equal(frame_fit$theta, matrix_fit$theta)
  expect_equal(frame_fit$theta_vcov, matrix_fit$theta_vcov)
  expect_identical(colnames(frame_fit$x0), c("x1", "x2"))
  if (requireNamespace("tibble", quietly = TRUE)) {
    args$data <- tibble::as_tibble(d)
    tibble_fit <- do.call(clate_dnn, args)
    expect_equal(tibble_fit$theta_vcov, matrix_fit$theta_vcov)
  }
})

test_that("CLATE input routes and missing columns fail explicitly", {
  a <- .clate_fit_arguments()
  d <- .clate_test_frame(a)
  args <- .clate_frame_arguments(a, d)
  expect_error(do.call(clate_dnn, c(args, list(X = a$X))), "data/X")
  expect_error(do.call(clate_dnn, c(args, list(Y = a$Y))), "data")
  for (nm in c("y_col", "d_col", "a_col", "e_col", "m_y0_col", "m_y1_col", "m_d0_col", "m_d1_col")) {
    bad <- args
    bad[[nm]] <- "absent"
    expect_error(do.call(clate_dnn, bad), "absent")
  }
  a$Y <- a$Y[-1L]
  expect_error(do.call(clate_dnn, a), "Y")
})

test_that("tabular CLATE preserves reordered selectors and singleton shapes", {
  a <- .clate_fit_arguments()
  d <- .clate_test_frame(a)
  args <- .clate_frame_arguments(a, d)
  original <- do.call(clate_dnn, args)
  args$x_cols <- c("x2", "x1")
  args$X0 <- args$X0[, 2:1, drop = FALSE]
  reordered <- do.call(clate_dnn, args)
  expect_equal(reordered$theta, original$theta)
  expect_equal(reordered$theta_vcov, original$theta_vcov)
  expect_identical(colnames(reordered$x0), c("x2", "x1"))
  args$x_cols <- "x1"
  args$X0 <- a$X0[1L, 1L, drop = FALSE]
  single <- do.call(clate_dnn, args)
  a$X <- a$X[, 1L, drop = FALSE]
  a$X0 <- args$X0
  expect_equal(single$theta, do.call(clate_dnn, a)$theta)
  expect_identical(dim(single$x0), c(1L, 1L))
  expect_length(single$component_vcov, 1L)
  args$data <- data.frame(d, Y = d$Y, check.names = FALSE)
  expect_error(do.call(clate_dnn, args), "unique column names")
})
