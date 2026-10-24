test_that("plqr_dnn computes direct PLQR sandwich from supplied density", {
  .test_rng_scope({
    d <- .make_plqr_inference_data()
    fit <- plqr_dnn(
      X = d$X, D = d$D, g = d$g, mD = d$mD, Y = d$Y,
      theta_eval = d$theta_eval, tau = d$tau, density = d$density,
      X0 = d$X0, method = "exact", k = 2L, s = 10L,
      se = TRUE, return = "W"
    )

    expect_s3_class(fit, "rdistnn_plqr_result")
    expect_equal(
      unname(fit$theta),
      matrix(d$theta_eval, nrow = length(d$theta_eval), ncol = nrow(d$X0))
    )
    expect_equal(length(fit$score_covariance), nrow(d$X0))
    expect_equal(length(fit$jacobian), nrow(d$X0))
    expect_equal(length(fit$theta_omega), nrow(d$X0))
    expect_equal(length(fit$theta_vcov), nrow(d$X0))
    expect_equal(dim(fit$theta_se), c(length(d$theta_eval), nrow(d$X0)))

    manual <- .manual_plqr_sandwich(d, fit$W[, 1L])
    expect_equal(unname(fit$score_covariance[[1L]]), unname(manual$S), tolerance = 1e-12)
    expect_equal(unname(fit$jacobian[[1L]]), unname(manual$J), tolerance = 1e-12)
    expect_equal(unname(fit$theta_omega[[1L]]), unname(manual$omega), tolerance = 1e-10)
    expect_equal(unname(fit$theta_vcov[[1L]]), unname(manual$vcov), tolerance = 1e-10)
    expect_equal(unname(fit$theta_se[, 1L]), unname(sqrt(diag(manual$vcov))), tolerance = 1e-10)
    expect_equal(unname(fit$neff), unname(1 / colSums(fit$W^2)), tolerance = 1e-12)
    expect_identical(fit$meta$theta_se_method, "direct_sandwich")
    expect_identical(fit$meta$moment_se_method, "none")
    expect_identical(fit$meta$theta_source, "supplied")

    txt <- capture.output(print(fit))
    expect_true(any(grepl("coefficient SEs available", txt, fixed = TRUE)))
    expect_true(any(grepl("coefficient SEs reported above", txt, fixed = TRUE)))
  }, seed = 818L)
})

test_that("direct sandwich internals reject unnormalized weights", {
  .test_rng_scope({
    d <- .make_plqr_inference_data(n = 40L, m = 1L)
    fit <- plqr_dnn(
      X = d$X, D = d$D, g = d$g, mD = d$mD, Y = d$Y,
      theta_eval = d$theta_eval, tau = d$tau, density = d$density,
      X0 = d$X0, method = "exact", k = 2L, s = 6L,
      se = TRUE, return = "W"
    )

    # Rescaled weight columns must fail clearly rather than silently producing
    # scale-dependent coefficient SEs.
    expect_error(
      RdistNN:::.plqr_direct_sandwich(
        d_mat = d$D, mD_mat = d$mD, density_vec = d$density,
        W = 2 * fit$W, tau = d$tau,
        coef_names = colnames(d$D), target_names = "target_1"
      ),
      "normalized weight columns",
      ignore.case = TRUE
    )

    # The jackknife-meat route shares the same internal guard.
    W_corrupt <- fit$W
    W_corrupt[1L, 1L] <- W_corrupt[1L, 1L] + 0.5
    expect_error(
      RdistNN:::.plqr_jackknife_meat_sandwich(
        x_mat = d$X, x0_mat = d$X0,
        score_mat = (d$tau - as.numeric(d$Y <= as.vector(d$D %*% d$theta_eval) + d$g)) * (d$D - d$mD),
        d_mat = d$D, mD_mat = d$mD, density_vec = d$density,
        W = W_corrupt,
        distance = "l2", method = "exact", k = 2L, s = 6L,
        normalize = TRUE, ctx = NULL, jackknife = "fast",
        precomputed_ranks = NULL,
        coef_names = colnames(d$D), target_names = "target_1"
      ),
      "normalized weight columns",
      ignore.case = TRUE
    )
  }, seed = 818L)
})
# Analytic inference does not delete observations, so the largest feasible
# weight order must not depend on the irrelevant jackknife selector.
test_that("analytic PLQR full-sample order is independent of jackknife selector", {
  .with_distance_tie_warning({
  args <- list(X = matrix(1:5), D = matrix(c(1, 2, 3, 2, 1)),
    mD = matrix(0, 5L, 1L), Y = c(-1, 1, -2, 2, 0), g = rep(0, 5L),
    density = rep(.5, 5L), theta_eval = .2, x0 = 2, s = 5L)
  default <- do.call(plqr_dnn, args)
  none <- do.call(plqr_dnn, c(args, list(jackknife = "none")))
  expect_identical(default$meta, none$meta)
  expect_identical(default$theta_se, none$theta_se)
  expect_equal(unname(default$theta_se), matrix(.5), tolerance = 1e-12)
  expect_equal(unname(vcov(default)), matrix(.25), tolerance = 1e-12)
  expect_null(default$moment_se)

  expect_error(do.call(plqr_dnn, c(args, list(theta_meat = "jackknife"))), "s < n")
  args$density <- NULL
  expect_error(do.call(plqr_dnn, args), "s < n")
  })
})
