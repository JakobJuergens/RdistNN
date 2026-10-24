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
    expect_error(diagnostics(bad_dx_dim), class = "rdistnn_state_error")

    bad_dx_value <- fit
    bad_dx_value$dx[1L, 1L] <- Inf
    expect_error(diagnostics(bad_dx_value), class = "rdistnn_state_error")
  }, seed = 777L)
})

.retained_diagnostic_extra_fits <- function() {
  common <- list(X = matrix(seq_len(8L)), X0 = matrix(c(.25, 9.25)),
    s = 3L, se = FALSE, return = c("W", "dx"))
  list(
    np = do.call(np_dnn, c(common, list(Y = seq_len(8L)^2))),
    cate = do.call(cate_dnn, c(common, list(pseudo_outcome = seq_len(8L)^2))),
    plqr = do.call(plqr_dnn, c(common, list(Y = seq_len(8L)^2,
      D = matrix(seq_len(8L)), mD = matrix(0, 8L), g = rep(0, 8L), theta_eval = 0)))
  )
}

test_that("retained diagnostic matrices reject malformed extras with refit guidance", {
  invalid <- list(matrix(1, 7L, 2L), matrix(NA_real_, 8L, 2L),
    matrix(Inf, 8L, 2L), matrix(TRUE, 8L, 2L), matrix(1 + 1i, 8L, 2L),
    rep(1, 16L), as.data.frame(matrix(1, 8L, 2L)))
  for (fit in .retained_diagnostic_extra_fits()) {
    for (field in c("W", "dx")) for (value in invalid) {
      bad <- fit
      bad[[field]] <- value
      err <- tryCatch(diagnostics(bad), error = identity)
      expect_s3_class(err, "rdistnn_state_error")
      if (inherits(err, "condition")) {
        expect_identical(err$code, "state_invalid")
        expect_identical(err$argument, paste0("object$", field))
        expect_match(conditionMessage(err), "Refit", fixed = TRUE)
      }
    }
  }
})

test_that("optional diagnostic extras use exact names and survive current serialization", {
  for (fit in .retained_diagnostic_extra_fits()) {
    reference <- diagnostics(fit)
    restored <- unserialize(serialize(fit, NULL))
    expect_equal(diagnostics(restored), reference)
    absent_dx <- restored
    absent_dx$dx <- NULL
    expect_true(all(!diagnostics(absent_dx)$distance_diagnostics_available))
    renamed_dx <- restored
    names(renamed_dx)[names(renamed_dx) == "dx"] <- "dx_old"
    expect_equal(diagnostics(renamed_dx), diagnostics(absent_dx))
    for (remove in c(FALSE, TRUE)) {
      missing_W <- restored
      if (remove) missing_W$W <- NULL else names(missing_W)[names(missing_W) == "W"] <- "W_old"
      err <- tryCatch(diagnostics(missing_W), error = identity)
      expect_s3_class(err, "rdistnn_extra_unavailable")
      if (inherits(err, "condition")) {
        expect_identical(err$argument, "object$W")
        expect_match(conditionMessage(err), 'return = "W"', fixed = TRUE)
      }
    }
  }
  # The numeric helper's own input boundary remains distinct from saved state.
  expect_error(RdistNN:::.dnn_weight_diagnostics(matrix(1 / 8, 8L, 2L),
    dx = matrix(NA_real_, 8L, 2L)), class = "rdistnn_input_error")
})

test_that("print.rdistnn_diagnostics runs", {
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
    expect_true(any(grepl("<rdistnn_diagnostics>", txt, fixed = TRUE)))
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
