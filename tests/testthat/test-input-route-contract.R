test_that("core input routes reject ignored direct inputs and selectors", {
  d <- data.frame(x = 1:8, y = 1:8)
  for (fun in list(np_dnn, np_tdnn)) {
    order_args <- if (identical(fun, np_dnn)) list(s = 3L) else list(s1 = 2L, s2 = 4L)
    args <- c(list(data = d, x_cols = "x", y_col = "y", x0 = 0, se = FALSE), order_args)
    expect_error(do.call(fun, c(args, list(X = matrix(1:8)))), "X")
    expect_error(do.call(fun, c(args, list(Y = 1:8))), "Y")
    expect_error(do.call(fun, c(list(X = matrix(1:8), Y = 1:8, x_cols = "x", x0 = 0,
      se = FALSE), order_args)), "data")
    args$y_col <- NULL; args$y <- d$y
    expect_s3_class(do.call(fun, args), "np_result")
  }
  expect_error(cate_dnn(data = d, x_cols = "x", X = matrix(1:8),
    pseudo_outcome = 1:8, x0 = 0, s = 3, se = FALSE), "X")
  expect_s3_class(cate_dnn(data = d, x_cols = "x", pseudo_outcome = 1:8,
    x0 = 0, s = 3, se = FALSE), "cate_result")
})
