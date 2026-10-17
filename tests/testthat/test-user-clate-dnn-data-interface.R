.clate_test_frame <- function(a) {
  data.frame(x1 = a$X[, 1L], x2 = a$X[, 2L], a[names(.clate_test_predictions())])
}

.clate_frame_arguments <- function(a, data) {
  list(data = data, x_cols = c("x1", "x2"), X0 = a$X0, s = a$s,
       y_col = "Y", d_col = "D", a_col = "A", e_col = "e",
       m_y0_col = "m_y0", m_y1_col = "m_y1", m_d0_col = "m_d0", m_d1_col = "m_d1")
}

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
