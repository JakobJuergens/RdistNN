test_that("diagnostics(rdistnn_plqr_result) reports residual-regressor Gram diagnostics", {
  .test_rng_scope({
    demo <- plqr_demo_data(n = 75L, p = 2L, p_d = 1L, tau = 0.5, seed = 21L)
    fit <- plqr_dnn(
      data = demo$data,
      x_cols = c("X1", "X2"),
      D_cols = "D1",
      g_col = "g_oracle",
      mD_cols = "mD1",
      y_col = "Y",
      theta_eval = demo$truth$theta_true,
      tau = 0.5,
      x0 = demo$x0,
      method = "exact",
      k = 2L,
      s = 4L,
      se = FALSE,
      return = "W"
    )
    d <- diagnostics(fit)

    d_resid <- demo$data$D1 - demo$data$mD1
    expected_gram <- as.numeric(crossprod(d_resid^2, fit$W))

    expect_s3_class(d, "rdistnn_diagnostics")
    expect_true(all(d$gram_available))
    expect_equal(d$gram_min_eig, expected_gram, tolerance = 1e-12)
    expect_equal(d$gram_max_eig, expected_gram, tolerance = 1e-12)
    expect_equal(d$gram_condition_number, rep(1, nrow(d)), tolerance = 1e-12)

    bad_state <- fit
    bad_state$.model <- NULL
    err <- tryCatch(diagnostics(bad_state), error = identity)
    expect_s3_class(err, "rdistnn_state_error")
    expect_identical(err$code, "state_invalid")
    expect_match(err$argument, "object$.model", fixed = TRUE)
    expect_match(conditionMessage(err), "Refit", fixed = TRUE)
  }, seed = 777L)
})

test_that("diagnostics(rdistnn_plqr_result) reports multi-column and singular Gram diagnostics", {
  .test_rng_scope({
    demo <- plqr_demo_data(n = 85L, p = 2L, p_d = 2L, tau = 0.5, seed = 24L)
    fit <- plqr_dnn(
      data = demo$data,
      x_cols = c("X1", "X2"),
      D_cols = c("D1", "D2"),
      g_col = "g_oracle",
      mD_cols = c("mD1", "mD2"),
      y_col = "Y",
      theta_eval = demo$truth$theta_true,
      tau = 0.5,
      x0 = demo$x0,
      method = "exact",
      k = 2L,
      s = 4L,
      se = FALSE,
      return = "W"
    )
    d <- diagnostics(fit)

    d_resid <- fit$.model$d_resid_mat
    expected <- lapply(seq_len(ncol(fit$W)), function(j) {
      G <- crossprod(d_resid, d_resid * as.numeric(fit$W[, j]))
      G <- (G + t(G)) / 2
      eigen(G, symmetric = TRUE, only.values = TRUE)$values
    })
    expect_equal(d$gram_min_eig, vapply(expected, min, numeric(1L)), tolerance = 1e-12)
    expect_equal(d$gram_max_eig, vapply(expected, max, numeric(1L)), tolerance = 1e-12)
    expect_true(all(is.finite(d$gram_condition_number)))

    singular <- plqr_dnn(
      X = fit$.model$x_mat, X0 = fit$x0, Y = demo$data$Y,
      D = fit$.model$d_mat[, c(1L, 1L)],
      mD = fit$.model$mD_mat[, c(1L, 1L)], g = demo$data$g_oracle,
      theta_eval = demo$truth$theta_true, tau = .5, k = 2L, s = 4L,
      se = FALSE, return = "W")
    d_singular <- diagnostics(singular)
    expect_true(all(d_singular$gram_min_eig <= d_singular$gram_rank_tol))
    expect_true(all(is.infinite(d_singular$gram_condition_number)))
  }, seed = 777L)
})

test_that("diagnostics(rdistnn_plqr_result) validates retained residual-regressor state", {
  .test_rng_scope({
    demo <- plqr_demo_data(n = 70L, p = 2L, p_d = 2L, tau = 0.5, seed = 25L)
    fit <- plqr_dnn(
      data = demo$data,
      x_cols = c("X1", "X2"),
      D_cols = c("D1", "D2"),
      g_col = "g_oracle",
      mD_cols = c("mD1", "mD2"),
      y_col = "Y",
      theta_eval = demo$truth$theta_true,
      tau = 0.5,
      x0 = demo$x0,
      method = "exact",
      k = 2L,
      s = 4L,
      se = FALSE,
      return = "W"
    )

    bad_rows <- fit
    bad_rows$.model$d_resid_mat <- bad_rows$.model$d_resid_mat[-1L, , drop = FALSE]
    error <- tryCatch(diagnostics(bad_rows), error = identity)
    expect_s3_class(error, "rdistnn_state_error")
    expect_identical(error$code, "state_invalid")
    expect_identical(error$argument, "object$.model$d_resid_mat")
    expect_match(conditionMessage(error), "70 by 2.*[Rr]efit")

    bad_value <- fit
    bad_value$.model$d_resid_mat[1L, 1L] <- NaN
    error <- tryCatch(diagnostics(bad_value), error = identity)
    expect_s3_class(error, "rdistnn_state_error")
    expect_identical(error$argument, "object$.model$d_resid_mat")
    expect_match(conditionMessage(error), "nonfinite.*[Rr]efit")
  }, seed = 777L)
})

test_that("diagnostics report target-level Jacobian conditioning for direct-sandwich fits", {
  .test_rng_scope({
    d <- .make_plqr_inference_data(n = 70L, m = 2L)
    fit <- plqr_dnn(
      X = d$X, D = d$D, g = d$g, mD = d$mD, Y = d$Y,
      theta_eval = d$theta_eval, tau = d$tau, density = d$density,
      X0 = d$X0, method = "exact", k = 2L, s = 8L,
      se = TRUE, return = "W"
    )
    dd <- diagnostics(fit)

    cols <- c(
      "jacobian_available", "jacobian_min_singular_value",
      "jacobian_condition_number", "jacobian_rank",
      "density_min_support", "density_mean_weighted",
      "density_median_weighted", "density_max_support",
      "theta_vcov_available"
    )
    expect_true(all(cols %in% names(dd)))
    expect_true(all(dd$jacobian_available))
    expect_true(all(is.finite(dd$jacobian_min_singular_value)))
    expect_true(all(dd$jacobian_min_singular_value > 0))
    expect_true(all(is.finite(dd$jacobian_condition_number)))
    expect_true(all(dd$jacobian_rank == length(d$theta_eval)))
    expect_true(all(dd$theta_vcov_available))
    expect_true(all(dd$density_min_support >= min(d$density)))
    expect_true(all(dd$density_max_support <= max(d$density)))
    expect_equal(
      dd$density_mean_weighted,
      as.numeric(crossprod(fit$W, d$density)),
      tolerance = 1e-12
    )
    expect_true(all(dd$density_mean_weighted >= dd$density_min_support))
    expect_true(all(dd$density_mean_weighted <= dd$density_max_support))
    expect_true(all(is.finite(dd$density_median_weighted)))
  }, seed = 818L)
})

test_that("diagnostics flag singular PLQR Jacobians and unavailable coefficient inference", {
  .test_rng_scope({
    d <- .make_plqr_inference_data(n = 60L, m = 1L)
    # D == mD makes the residualized design identically zero, so J is singular.
    d_singular <- d
    d_singular$mD <- d_singular$D
    fit <- suppressWarnings(plqr_dnn(
      X = d_singular$X, D = d_singular$D, g = d_singular$g,
      mD = d_singular$mD, Y = d_singular$Y,
      theta_eval = d_singular$theta_eval, tau = d_singular$tau,
      density = d_singular$density,
      X0 = d_singular$X0, method = "exact", k = 2L, s = 8L,
      se = TRUE, return = "W"
    ))
    dd <- diagnostics(fit)

    expect_true(all(dd$jacobian_available))
    expect_true(all(dd$jacobian_condition_number == Inf))
    expect_true(all(dd$jacobian_rank == 0L))
    expect_false(any(dd$theta_vcov_available))

    # Fits without density report the Jacobian block as unavailable.
    fit_no_density <- plqr_dnn(
      X = d$X, D = d$D, g = d$g, mD = d$mD, Y = d$Y,
      theta_eval = d$theta_eval, tau = d$tau,
      X0 = d$X0, method = "exact", k = 2L, s = 8L,
      se = FALSE, return = "W"
    )
    dd_no <- diagnostics(fit_no_density)
    expect_false(any(dd_no$jacobian_available))
    expect_true(all(is.na(dd_no$jacobian_min_singular_value)))
    expect_false(any(dd_no$theta_vcov_available))
  }, seed = 818L)
})

test_that("nearly singular PLQR Jacobians report large condition numbers", {
  .test_rng_scope({
    d <- .make_plqr_inference_data(n = 60L, m = 1L)
    # Make the second residualized column an almost-exact copy of the first so
    # J is finite-rank but ill-conditioned rather than exactly singular.
    eps <- 1e-7
    d$D[, 2L] <- d$D[, 1L] + eps * d$D[, 2L]
    d$mD[, 2L] <- d$mD[, 1L] + eps * d$mD[, 2L]
    fit <- plqr_dnn(
      X = d$X, D = d$D, g = d$g, mD = d$mD, Y = d$Y,
      theta_eval = d$theta_eval, tau = d$tau, density = d$density,
      X0 = d$X0, method = "exact", k = 2L, s = 8L,
      se = TRUE, return = "W"
    )
    dd <- diagnostics(fit)
    expect_true(all(dd$jacobian_available))
    expect_true(all(dd$jacobian_condition_number > 1e6))
  }, seed = 818L)
})
test_that("diagnostic support includes every represented positive localization weight", {
  n <- 80L
  X <- matrix(seq_len(n))
  density <- rep(1, n)
  density[20L] <- 1000
  fits <- list(
    np = np_dnn(X = X, Y = seq_len(n), x0 = 0, s = 50L, se = FALSE),
    cate = cate_dnn(X = X, pseudo_outcome = seq_len(n), x0 = 0, s = 50L, se = FALSE),
    ratio = affine_ratio_dnn(X = X, numerator = seq_len(n), denominator = rep(1, n),
      x0 = 0, s = 50L, se = FALSE),
    plqr = plqr_dnn(X = X, Y = rep(1, n), D = matrix(1, n), mD = matrix(0, n),
      g = rep(0, n), theta_eval = 0, density = density, x0 = 0, s = 50L)
  )
  for (fit in fits) {
    expect_identical(sum(fit$W[, 1L] != 0), 31L)
    expect_true(fit$W[20L, 1L] > 0 && fit$W[20L, 1L] < sqrt(.Machine$double.eps))
    expect_identical(diagnostics(fit)$n_nonzero_weight, 31L)
  }
  dd <- diagnostics(fits$plqr)
  expect_identical(dd$density_min_support, 1)
  expect_identical(dd$density_max_support, 1000)
  expect_identical(dd$density_median_weighted,
    RdistNN:::.dnn_weighted_quantile(density, fits$plqr$W[, 1L], .5))
})
