test_that("retained selectors reproduce duplicated and matrix-valued data columns", {
  d <- data.frame(x = 1:8, x.1 = 8:1, y = 1:8)
  for (family in c("np", "cate", "tdnn")) {
    fit <- switch(family,
      np = np_dnn(data = d, x_cols = c("x", "x"), y_col = "y", x0 = c(0, 0), s = 3, se = FALSE),
      cate = cate_dnn(data = d, x_cols = c("x", "x"), pseudo_outcome = d$y, x0 = c(0, 0), s = 3, se = FALSE),
      tdnn = np_tdnn(data = d, x_cols = c("x", "x"), y_col = "y", x0 = c(0, 0), s1 = 2, s2 = 3, se = FALSE))
    expect_equal(unname(predict(fit, newdata = data.frame(x = 0, x.1 = 100))),
      unname(predict(fit, x0 = c(0, 0))), tolerance = 0)
    expect_no_error(predict(fit, newdata = data.frame(x = 0)))
  }
  d$pair <- I(cbind(x = d$x, z = d$x.1))
  for (family in c("np", "cate", "tdnn")) {
    fit <- switch(family,
      np = np_dnn(data = d, x_cols = "pair", y_col = "y", x0 = c(0, 0), s = 3, se = FALSE),
      cate = cate_dnn(data = d, x_cols = "pair", pseudo_outcome = d$y, x0 = c(0, 0), s = 3, se = FALSE),
      tdnn = np_tdnn(data = d, x_cols = "pair", y_col = "y", x0 = c(0, 0), s1 = 2, s2 = 3, se = FALSE))
    expect_equal(unname(predict(fit, newdata = d)), unname(predict(fit, X0 = as.matrix(d$pair))))
  }
})
