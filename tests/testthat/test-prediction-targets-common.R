.common_prediction_fits <- function() {
  d <- data.frame(junk = 101:107, a = c(-2, -1, 0, 0, 1, 2, 4),
    b = c(3, 1, 0, 0, 3, 2, -1), y = c(2, -1, 3, 7, -2, 4, 1),
    den = 1, arm = c(0, 1, 0, 1, 0, 1, 0), propensity = .5,
    my0 = 0, my1 = 0, md0 = 0, md1 = 1)
  X0 <- rbind(c(-.2, .1), c(1.7, 1), c(0, 0))
  list(
    np = np_dnn(data = d, x_cols = 2:3, y_col = "y", X0 = X0, s = 3, se = FALSE),
    cate = cate_dnn(data = d, x_cols = 2:3, pseudo_outcome = d$y,
      X0 = X0, s = 3, se = FALSE),
    ratio = affine_ratio_dnn(data = d, x_cols = 2:3, numerator_col = "y",
      denominator_col = "den", X0 = X0, s = 3, se = FALSE),
    clate = clate_dnn(data = d, x_cols = 2:3, y_col = "y", d_col = "arm", a_col = "arm",
      e_col = "propensity", m_y0_col = "my0", m_y1_col = "my1",
      m_d0_col = "md0", m_d1_col = "md1", X0 = X0, s = 3, se = FALSE)
  )
}

test_that("prediction adapters share canonical reporting-target identities", {
  .with_localization_policy_warnings({
  fits <- .common_prediction_fits()
  for (f in fits) {
    at <- f$x0
    nd <- data.frame(b = at[, 2], extra = NA_real_, a = at[, 1])
    estimate <- unname(coef(f))
    expect_equal(predict(f, newdata = nd), estimate)
    expect_equal(predict(f, X0 = at), estimate)
    expect_equal(predict(f, x0 = at[1, ]), estimate[1])
    expect_equal(predict(unserialize(serialize(f, NULL)), newdata = nd), estimate)
  }
  })
})

test_that("prediction adapters validate the selected target interface consistently", {
  .with_localization_policy_warnings({
  fits <- .common_prediction_fits()
  for (f in fits) {
    at <- f$x0
    expect_error(predict(f, newdata = at, X0 = at), "newdata/X0/x0")
    expect_error(predict(f, X0 = at, x0 = at[1, ]), "X0/x0")
    expect_error(predict(f, newdata = matrix(1, 1, 1)), "newdata")
    expect_error(predict(f, newdata = matrix(numeric(), 0, 2)), "newdata.*at least 1 row")
    bad <- at
    bad[1, 1] <- Inf
    expect_error(predict(f, newdata = bad), "newdata.*non-finite")
    expect_error(predict(f, x0 = c(NA_real_, 0)), "x0.*NA")
    expect_error(predict(f, newdata = data.frame(a = 0)), "b.*x_cols")
  }
  })
})
