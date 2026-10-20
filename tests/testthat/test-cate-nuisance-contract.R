test_that("propensity limits and partial nuisance routes fail early", {
  d <- data.frame(x = 1:8, y = 1:8, d = rep(0:1, 4), e = .5, m = 0)
  expect_error(cate_dnn(data = d, x_cols = "x", y_col = "y", t_col = "d",
    e_col = "e", m0_col = "m", x0 = 0, s = 3), "m0_col.*m1_col")
  expect_error(cate_dnn(X = matrix(1:8), pseudo_outcome = 1:8, x0 = 0,
    s = 3, clip = 1e-20), "clip")
  expect_error(iv_augmented_signals(Y = 1:2, D = 0:1, A = 0:1, e = 0:1,
    m_y0 = 1:2, m_y1 = 1:2, m_d0 = c(0, 0), m_d1 = c(1, 1), clip = 1e-20), "clip")
  expect_no_error(cate_dnn(data = d, x_cols = "x", pseudo_outcome = d$y,
    m0_col = "absent", x0 = 0, s = 3, se = FALSE))
})

test_that("AIPW retains its population identity with correct arm means", {
  d <- y <- c(0, 0, 0, 1)
  for (e in c(.24, .26, .30)) {
    expect_equal(mean(cate_rho_aipw(y, d, rep(1, 4), rep(0, 4), rep(e, 4))), 1)
  }
})
