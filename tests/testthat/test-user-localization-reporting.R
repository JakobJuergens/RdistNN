test_that("fit and summary printing report localization fractions without changing fits", {
  X <- matrix(1:8, ncol = 1L)
  y <- c(2, 1, 4, 3, 6, 5, 7, 9)
  fits <- list(np_dnn(X = X, x0 = 2, Y = y, s = 3),
    cate_dnn(X = X, x0 = 2, pseudo_outcome = y, s = 3),
    affine_ratio_dnn(X = X, x0 = 2, numerator = y, denominator = rep(1, 8), s = 3),
    np_tdnn(X = X, x0 = 2, Y = y, s1 = 2, s2 = 4, dim = 1))
  for (f in fits) {
    before <- serialize(f, NULL)
    expect_match(paste(capture.output(print(f)), collapse = " "), "s[12]?/n")
    expect_match(paste(capture.output(print(summary(f))), collapse = " "), "s[12]?/n")
    expect_identical(serialize(f, NULL), before)
  }
  expect_match(paste(capture.output(print(fits[[1L]])), collapse = " "), "0.375")
})

test_that("causal diagnostics expose available arm and contribution summaries", {
  d <- cate_demo_data(n = 40, seed = 5)$data
  f <- cate_dnn(data = d, x_cols = c("X1", "X2"), x0 = c(0, 0),
    y_col = "Y", t_col = "D", e_col = "e_oracle",
    m0_col = "m0_oracle", m1_col = "m1_oracle", s = 5)
  diag <- diagnostics(f)
  before <- serialize(diag, NULL)
  output <- paste(capture.output(print(diag)), collapse = " ")
  expect_match(output, "treated ESS")
  expect_match(output, "control ESS")
  expect_match(output, "propensity")
  expect_match(output, "absolute score")
  expect_identical(serialize(diag, NULL), before)
  f <- cate_dnn(X = matrix(1:8, ncol = 1L), x0 = 2, pseudo_outcome = 1:8, s = 3)
  expect_match(paste(capture.output(print(diagnostics(f))), collapse = " "), "unavailable")
})
