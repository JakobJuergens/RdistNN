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

test_that("PLQR diagnostics preserve a large finite scalar Gram", {
  fit_at_scale <- function(scale) plqr_dnn(
    X = matrix(seq_len(8L)), x0 = 0, D = matrix(scale, 8L, 1L),
    mD = matrix(0, 8L, 1L), Y = rep(c(-1, 1), 4L), g = rep(0, 8L),
    theta_eval = 0, density = rep(.5, 8L), s = 3L, se = TRUE, return = "W"
  )
  ordinary <- fit_at_scale(1)
  large <- fit_at_scale(1e154)
  expect_true(all(is.finite(vcov(large))) && all(vcov(large) > 0))
  expect_equal(unname(vcov(large)) * 1e308, unname(vcov(ordinary)), tolerance = 1e-12)
  ordinary_diagnostics <- diagnostics(ordinary)
  expect_no_warning(large_diagnostics <- diagnostics(large))
  expect_true(all(is.finite(large_diagnostics$gram_min_eig)))
  expect_true(all(is.finite(large_diagnostics$gram_max_eig)))
  expect_equal(large_diagnostics$gram_min_eig / 1e308,
    ordinary_diagnostics$gram_min_eig, tolerance = 1e-14)
  expect_equal(large_diagnostics$gram_max_eig / 1e308,
    ordinary_diagnostics$gram_max_eig, tolerance = 1e-14)
  for (field in c("gram_condition_number", "jacobian_condition_number", "jacobian_rank")) {
    expect_identical(large_diagnostics[[field]], ordinary_diagnostics[[field]])
  }
  restored <- unserialize(serialize(large, NULL))
  expect_identical(diagnostics(restored)$gram_min_eig, large_diagnostics$gram_min_eig)
})

test_that("PLQR diagnostics mark an out-of-range Gram unavailable", {
  fit_at_scale <- function(regressor, density) plqr_dnn(
    X = matrix(seq_len(8L)), x0 = 0, D = matrix(regressor, 8L, 1L),
    mD = matrix(0, 8L, 1L), Y = rep(c(-1, 1), 4L), g = rep(0, 8L),
    theta_eval = 0, density = rep(density, 8L), s = 3L, se = TRUE, return = "W"
  )
  ordinary <- fit_at_scale(1, 1)
  large <- fit_at_scale(2e154, 1e-100)
  expect_true(all(is.finite(vcov(large))) && all(vcov(large) > 0))
  expect_equal(unname(vcov(large)) * 4e108, unname(vcov(ordinary)), tolerance = 1e-12)
  expect_no_warning(d <- diagnostics(large))
  expect_identical(d$gram_available, FALSE)
  for (field in c("gram_min_eig", "gram_max_eig", "gram_condition_number")) {
    expect_true(all(is.na(d[[field]])))
  }
  expect_identical(d$jacobian_available, TRUE)
  expect_identical(d$jacobian_rank, 1L)
  expect_identical(d$jacobian_condition_number, 1)
  expect_equal(d$jacobian_min_singular_value / 4e208, 1, tolerance = 1e-14)
  expect_identical(d$theta_vcov_available, TRUE)
  restored <- unserialize(serialize(large, NULL))
  expect_identical(diagnostics(restored)$gram_available, FALSE)
  expect_identical(diagnostics(restored)$jacobian_available, TRUE)
})

test_that("public PLQR inference retains a finite dense norm-overflow bread", {
  .with_localization_policy_warnings({
    a <- 1e154 * sqrt(1.8)
    b <- 1e154 * sqrt(.2)
    D <- rbind(c(a, a), c(-a, -a), c(b, -b), c(-b, b))
    fit <- plqr_dnn(X = matrix(1:4), x0 = 0, D = D, mD = D * 0,
      Y = c(-1, 1, -1, 1), g = rep(0, 4L), theta_eval = c(0, 0),
      density = rep(1, 4L), s = 1L, return = "W")
    expected <- matrix(c(25, -20, -20, 25), 2L) / 144
    expect_true(all(is.finite(fit$jacobian[[1L]])))
    expect_true(all(is.finite(vcov(fit))))
    expect_equal(unname(vcov(fit)) / 1e-308, expected, tolerance = 1e-12)
    d <- diagnostics(fit)
    expect_identical(d$gram_available, FALSE)
    expect_true(is.na(d$gram_condition_number))
    expect_identical(d$jacobian_available, TRUE)
    expect_identical(d$jacobian_rank, 2L)
    expect_equal(d$jacobian_condition_number, 9, tolerance = 1e-13)
    expect_equal(d$jacobian_min_singular_value / 2e307, 1, tolerance = 1e-13)
    expect_identical(d$theta_vcov_available, TRUE)
    restored <- unserialize(serialize(fit, NULL))
    expect_identical(vcov(restored), vcov(fit))
    expect_identical(diagnostics(restored)$jacobian_rank, 2L)
    expect_s3_class(summary(restored), "rdistnn_summary_plqr_result")
    prediction <- predict(restored, x0 = 0, se.fit = TRUE)
    expect_equal(unname(prediction$se.fit), unname(fit$theta_se), tolerance = 1e-13)
  })
})

test_that("PLQR Jacobian diagnostics preserve conditioning beyond singular-value range", {
  J <- 1.5e308 * matrix(c(1, 1, 1, -1), 2L)
  object <- list(jacobian = list(J), theta_vcov = list(diag(1e-309, 2L)),
    .model = list(density = rep(1, 2L)))
  d <- RdistNN:::.plqr_jacobian_diagnostics(object, matrix(.5, 2L, 1L))
  expect_identical(d$jacobian_available, TRUE)
  expect_identical(d$jacobian_rank, 2L)
  expect_equal(d$jacobian_condition_number, 1, tolerance = 1e-14)
  expect_true(is.na(d$jacobian_min_singular_value))
  expect_identical(d$theta_vcov_available, TRUE)
})
