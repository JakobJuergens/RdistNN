# Shared diagnostics validation and cross-estimator reporting.

test_that("general diagnostics reject unsupported objects and malformed retained matrices", {
  .test_rng_scope({
    n <- 35L
    x_mat <- matrix(rnorm(n * 2L), n, 2L)
    y_vec <- x_mat[, 1L] + rnorm(n, sd = 0.2)
    fit <- np_dnn(
      X = x_mat, x0 = c(0, 0), Y = y_vec,
      method = "exact", k = 2L, s = 3L, se = FALSE,
      return = c("W", "dx")
    )

    expect_error(diagnostics(list()), "not implemented for class 'list'")

    bad_w <- fit
    bad_w$W[1L, 1L] <- NA_real_
    expect_error(diagnostics(bad_w), "object\\$W.*[Rr]efit")

    bad_dx_dim <- fit
    bad_dx_dim$dx <- fit$dx[-1L, , drop = FALSE]
    expect_error(diagnostics(bad_dx_dim), "same dimensions as `W`")

    bad_dx_value <- fit
    bad_dx_value$dx[1L, 1L] <- Inf
    expect_error(diagnostics(bad_dx_value), "dx cannot contain non-finite values")
  }, seed = 777L)
})

test_that("print.dnn_diagnostics runs", {
  .test_rng_scope({
    n <- 30L
    x_mat <- matrix(rnorm(n * 2L), n, 2L)
    y_vec <- x_mat[, 1L] + rnorm(n, sd = 0.2)
    fit <- np_dnn(
      X = x_mat, x0 = c(0, 0), Y = y_vec,
      method = "exact", k = 2L, s = 3L, se = FALSE, return = "W"
    )
    d <- diagnostics(fit)
    txt <- capture.output(print(d))
    expect_true(any(grepl("<dnn_diagnostics>", txt, fixed = TRUE)))
    expect_true(any(grepl("outcome:", txt, fixed = TRUE)))

    demo_cate <- cate_demo_data(n = 45L, p = 2L, seed = 26L)
    fit_cate <- cate_dnn(
      data = demo_cate$data,
      x_cols = c("X1", "X2"),
      y_col = "Y",
      t_col = "D",
      m0_col = "m0_oracle", m1_col = "m1_oracle",
      e_col = "e_oracle",
      x0 = demo_cate$x0,
      method = "exact",
      k = 2L,
      s = 4L,
      se = FALSE,
      return = "W"
    )
    txt_cate <- capture.output(print(diagnostics(fit_cate)))
    expect_true(any(grepl("overlap:", txt_cate, fixed = TRUE)))

    demo_plqr <- plqr_demo_data(n = 55L, p = 2L, p_d = 1L, seed = 27L)
    fit_plqr <- plqr_dnn(
      data = demo_plqr$data,
      x_cols = c("X1", "X2"),
      D_cols = "D1",
      g_col = "g_oracle",
      mD_cols = "mD1",
      y_col = "Y",
      theta_eval = demo_plqr$truth$theta_true,
      x0 = demo_plqr$x0,
      method = "exact",
      k = 2L,
      s = 4L,
      se = FALSE,
      return = "W"
    )
    txt_plqr <- capture.output(print(diagnostics(fit_plqr)))
    expect_true(any(grepl("gram   :", txt_plqr, fixed = TRUE)))
  }, seed = 777L)
})
