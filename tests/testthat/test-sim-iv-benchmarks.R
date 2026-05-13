test_that(".sim_metric_iv handles overidentification diagnostics", {
  metric <- getFromNamespace(".sim_metric_iv", "RdistNN")

  just <- metric(theta_hat = 1, se_hat = 0.2, theta0 = 1, alpha = 0.05,
                 neff = 20, J = NA_real_, df = NA_integer_)
  expect_true(is.na(just[["J"]]))
  expect_true(is.na(just[["df"]]))
  expect_true(is.na(just[["pval_J"]]))
  expect_true(is.na(just[["reject_J"]]))

  over <- metric(theta_hat = 1, se_hat = 0.2, theta0 = 1, alpha = 0.05,
                 neff = 20, J = 1.5, df = 2L)
  expect_equal(over[["df"]], 2)
  expect_true(is.finite(over[["pval_J"]]))
  expect_true(over[["reject_J"]] %in% c(0, 1))
})

test_that(".method_iv_oracle matches manual weighted efficient GMM", {
  meth <- getFromNamespace(".method_iv_oracle", "RdistNN")
  dgp <- sim_iv_data("dgp2_iv_homog_overid", n = 160L, kx = 2L, q = 2L, seed = 11L)
  x_cols <- c("X1", "X2")
  k <- 2L
  s <- 12L
  alpha <- 0.05

  res <- meth(dgp = dgp, x_cols = x_cols, k = k, s = s, alpha = alpha)

  wt <- weights_dnn(
    X = as.matrix(dgp$data[, x_cols]),
    x0 = dgp$x0,
    k = k,
    s = s,
    method = "exact",
    return = "W"
  )
  w <- as.numeric(wt$W[, 1L])
  Z_resid <- as.matrix(dgp$data[, c("Z1", "Z2")]) -
    as.matrix(dgp$data[, c("piZ1_oracle", "piZ2_oracle")])
  A_obs <- Z_resid * (dgp$data$Y - dgp$data$muY_oracle)
  B_obs <- Z_resid * (dgp$data$D - dgp$data$muD_oracle)
  a_bar <- as.numeric(crossprod(A_obs, w))
  b_bar <- as.numeric(crossprod(B_obs, w))
  theta_1 <- sum(b_bar * a_bar) / sum(b_bar^2)
  g1 <- A_obs - B_obs * theta_1
  Omega <- crossprod(g1, g1 * w)
  Omega <- (Omega + t(Omega)) / 2
  Omega_inv <- solve(Omega)
  info <- as.numeric(t(b_bar) %*% Omega_inv %*% b_bar)
  theta_hat <- as.numeric(t(b_bar) %*% Omega_inv %*% a_bar / info)
  neff <- 1 / sum(w^2)
  se_hat <- sqrt(1 / (neff * info))
  gbar <- a_bar - b_bar * theta_hat
  J <- as.numeric(neff * t(gbar) %*% Omega_inv %*% gbar)

  expect_equal(res[["theta_hat"]], theta_hat, tolerance = 1e-12)
  expect_equal(res[["se"]], se_hat, tolerance = 1e-12)
  expect_equal(res[["neff"]], neff, tolerance = 1e-12)
  expect_equal(res[["J"]], J, tolerance = 1e-12)
  expect_equal(res[["df"]], 1)
})

test_that(".method_iv_oracle returns finite metrics and NA J for q = 1", {
  meth <- getFromNamespace(".method_iv_oracle", "RdistNN")
  dgp <- sim_iv_data("dgp1_iv_justid_linear", n = 140L, kx = 2L, q = 3L, seed = 12L)
  res <- meth(dgp = dgp, x_cols = c("X1", "X2"), k = 2L, s = 10L, alpha = 0.05)

  expect_true(is.finite(res[["theta_hat"]]))
  expect_true(is.finite(res[["se"]]))
  expect_true(res[["se"]] > 0)
  expect_true(is.na(res[["J"]]))
  expect_true(is.na(res[["df"]]))
})

test_that("classical IV benchmark methods return full metric schema", {
  methods <- c(".method_2sls_global", ".method_2sls_interact", ".method_local_2sls")
  dgp <- sim_iv_data("dgp1_iv_justid_linear", n = 180L, kx = 2L, q = 1L, seed = 13L)
  expected <- c(
    "theta_hat", "theta0", "bias", "sq_err", "se", "cover", "ci_len",
    "neff", "J", "df", "pval_J", "reject_J"
  )

  for (name in methods) {
    fn <- getFromNamespace(name, "RdistNN")
    res <- fn(dgp = dgp, x_cols = c("X1", "X2"), alpha = 0.05)
    expect_named(res, expected)
    expect_true(is.finite(res[["theta_hat"]]), info = name)
    expect_true(is.finite(res[["se"]]), info = name)
    expect_true(res[["cover"]] %in% c(0, 1), info = name)
    expect_true(is.na(res[["J"]]), info = name)
  }
})

test_that("DNN IV benchmark methods return full metric schema", {
  expected <- c(
    "theta_hat", "theta0", "bias", "sq_err", "se", "cover", "ci_len",
    "neff", "J", "df", "pval_J", "reject_J"
  )
  dgp <- sim_iv_data("dgp2_iv_homog_overid", n = 120L, kx = 2L, q = 2L, seed = 14L)
  methods <- c(".method_dnn_iv_oracle", ".method_dnn_iv_fullsample",
               ".method_dnn_iv_crossfit")

  for (name in methods) {
    fn <- getFromNamespace(name, "RdistNN")
    res <- fn(
      dgp = dgp,
      x_cols = c("X1", "X2"),
      k = 2L,
      s = 10L,
      alpha = 0.05,
      n_folds = 3L
    )
    expect_named(res, expected)
    expect_true(is.finite(res[["theta_hat"]]), info = name)
    expect_true(is.finite(res[["se"]]), info = name)
    expect_equal(res[["df"]], 1, info = name)
    expect_true(is.finite(res[["J"]]), info = name)
  }
})

test_that("Mahalanobis DNN IV benchmark works on correlated instruments", {
  fn <- getFromNamespace(".method_dnn_iv_oracle_mahal", "RdistNN")
  dgp <- sim_iv_data("dgp5_iv_correlated_instruments", n = 140L, kx = 2L,
                     q = 2L, seed = 15L)
  res <- fn(dgp = dgp, x_cols = c("X1", "X2"), k = 2L, s = 11L, alpha = 0.05)

  expect_true(is.finite(res[["theta_hat"]]))
  expect_true(is.finite(res[["se"]]))
  expect_equal(res[["df"]], 1)
  expect_true(is.finite(res[["pval_J"]]))
})

test_that(".get_sim_iv_method dispatches registered methods", {
  getter <- getFromNamespace(".get_sim_iv_method", "RdistNN")
  expect_true(is.function(getter("iv_oracle")))
  expect_true(is.function(getter("2sls_global")))
  expect_error(getter("bad_method"), "Unknown IV simulation method")
})
