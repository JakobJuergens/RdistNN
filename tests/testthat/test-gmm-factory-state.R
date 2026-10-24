test_that("DR-IV model factories capture column selectors before setup", {
  dat <- data.frame(Y = 1:3, D = 4:6, Z = c(0, 1, 0), pi = .5,
                    muY = 0, muD = 1)
  y <- "Y"
  d <- "D"
  z <- "Z"
  pi <- "pi"
  mu_y <- "muY"
  mu_d <- "muD"
  model <- gmm_model_dr_iv(y_col = y, d_col = d, z_cols = z,
    pi_cols = pi, muY_col = mu_y, muD_col = mu_d)
  y <- d <- z <- pi <- mu_y <- mu_d <- "missing"
  ctx <- model$setup(dat, n = nrow(dat))
  expect_identical(ctx$y, as.numeric(dat$Y))
  expect_identical(ctx$d, as.numeric(dat$D))
  expect_equal(ctx$z_resid, as.matrix(dat["Z"]) - as.matrix(dat["pi"]),
    ignore_attr = TRUE)
  expect_identical(ctx$muY, dat$muY)
  expect_identical(ctx$muD, dat$muD)
})
