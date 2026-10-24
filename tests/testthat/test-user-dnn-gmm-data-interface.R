# GMM matrix/data-frame equivalence and column selectors.
test_that("gmm_dnn DR-IV data-frame and matrix interfaces agree", {
  dgp <- sim_iv_data("dgp2_iv_homog_overid", n = 120L, kx = 2L, q = 2L, seed = 2L)
  model <- .test_dr_iv_model(dgp$q)
  x_cols <- c("X1", "X2")

  fit_df <- gmm_dnn(
    data = dgp$data,
    x_cols = x_cols,
    x0 = dgp$x0,
    model = model,
    k = 2L,
    s = 10L,
    return = c("ess", "omega", "jacobian")
  )
  fit_mat <- gmm_dnn(
    X = as.matrix(dgp$data[, x_cols]),
    model_data = dgp$data,
    x0 = dgp$x0,
    model = model,
    k = 2L,
    s = 10L,
    return = "ess"
  )

  expect_s3_class(fit_df, "rdistnn_gmm_result")
  expect_equal(fit_df$theta, fit_mat$theta, tolerance = 1e-8)
  expect_equal(fit_df$se, fit_mat$se, tolerance = 1e-8)
  expect_true(is.finite(fit_df$neff[1L]))
  expect_true(is.matrix(fit_df$omega[[1L]]))
  expect_true(is.matrix(fit_df$jacobian[[1L]]))
})

test_that("gmm_dnn supports numeric columns and tibble data like data frames", {
  dgp <- sim_iv_data("dgp2_iv_homog_overid", n = 110L, kx = 2L, q = 2L, seed = 24L)
  nm <- names(dgp$data)
  model_names <- .test_dr_iv_model(dgp$q)
  model_idx <- gmm_model_dr_iv(
    y_col = match("Y", nm),
    d_col = match("D", nm),
    z_cols = match(c("Z1", "Z2"), nm),
    pi_cols = match(c("piZ1_oracle", "piZ2_oracle"), nm),
    muY_col = match("muY_oracle", nm),
    muD_col = match("muD_oracle", nm)
  )

  fit_names <- gmm_dnn(
    data = dgp$data,
    x_cols = c("X1", "X2"),
    x0 = dgp$x0,
    model = model_names,
    k = 2L,
    s = 9L
  )
  fit_idx <- gmm_dnn(
    data = dgp$data,
    x_cols = match(c("X1", "X2"), nm),
    x0 = dgp$x0,
    model = model_idx,
    k = 2L,
    s = 9L
  )

  expect_equal(fit_idx$theta, fit_names$theta, tolerance = 1e-8)
  expect_equal(fit_idx$se, fit_names$se, tolerance = 1e-8)

  testthat::skip_if_not_installed("tibble")
  fit_tbl <- gmm_dnn(
    data = tibble::as_tibble(dgp$data),
    x_cols = c("X1", "X2"),
    x0 = dgp$x0,
    model = model_names,
    k = 2L,
    s = 9L
  )
  expect_equal(fit_tbl$theta, fit_names$theta, tolerance = 1e-8)
  expect_equal(fit_tbl$se, fit_names$se, tolerance = 1e-8)
})

test_that("gmm_dnn rejects ambiguous data-frame column-name schemas", {
  dgp <- sim_iv_data("dgp2_iv_homog_overid", n = 70L, kx = 2L, q = 2L, seed = 31L)
  nms <- names(dgp$data)
  model <- .test_dr_iv_model(dgp$q)
  df_dup <- dgp$data
  names(df_dup)[match("X2", nms)] <- "X1"

  expect_error(
    gmm_dnn(
      data = df_dup,
      x_cols = match(c("X1", "X2"), nms),
      x0 = dgp$x0,
      model = model,
      k = 1L, s = 2L, se = FALSE
    ),
    "unique column names",
    ignore.case = TRUE
  )
})
