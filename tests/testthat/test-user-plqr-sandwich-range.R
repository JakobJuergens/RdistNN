.plqr_range_fit <- function(theta_meat = "analytic", D = 1e154,
                            density = 100, mD = 0, se = TRUE) {
  plqr_dnn(X = matrix(seq_len(8L)), x0 = 4,
    Y = seq_len(8L), D = matrix(rep(D, 8L)),
    g = rep(0, 8L), mD = matrix(rep(mD, 8L)), theta_eval = 0,
    density = rep(density, 8L), s = 3L, se = se,
    theta_meat = theta_meat)
}

test_that("nonfinite PLQR bread cannot become a zero coefficient covariance", {
  .with_distance_tie_warning({
  for (meat in c("analytic", "jackknife")) {
    expect_warning(fit <- .plqr_range_fit(theta_meat = meat),
      "Jacobian.*outside the numeric range.*[Rr]escale")
    expect_true(all(is.na(fit$theta_se)))
    expect_true(all(is.na(fit$theta_vcov[[1L]])))
    expect_true(all(is.na(fit$theta_omega[[1L]])))
    expect_true(all(is.finite(fit$moment)))
    expect_warning(pred <- predict(fit, x0 = 4, se.fit = TRUE),
      "Jacobian.*outside the numeric range.*[Rr]escale")
    expect_true(all(is.na(pred$se.fit)))
  }
  # The analytic coefficient covariance itself is positive and representable;
  # this guard reports the unsupported raw-bread range instead of false zero.
  normalized <- .plqr_range_fit(D = 1, density = 100)
  expected_se <- normalized$theta_se[1L, 1L] / 1e154
  expect_gt(expected_se, 0)
  expect_gt(expected_se^2, 0)
  })
})

test_that("PLQR prediction reports nonfinite bread as unavailable uncertainty", {
  .with_distance_tie_warning({
  fit <- .plqr_range_fit(se = FALSE)
  expect_warning(pred <- predict(fit, x0 = 4, se.fit = TRUE,
    interval = "confidence"), "Jacobian.*outside the numeric range.*[Rr]escale")
  expect_true(all(is.na(pred$se.fit)))
  expect_true(all(is.na(pred$lwr)))
  expect_true(all(is.na(pred$upr)))
  expect_identical(as.numeric(pred$fit), 0)
  })
})

test_that("finite PLQR bread preserves ordinary uncertainty and zero meat", {
  .with_distance_tie_warning({
  ordinary <- .plqr_range_fit(D = 2, density = 1)
  expected <- .5 / (2 * sqrt(ordinary$neff))
  expect_equal(as.numeric(ordinary$theta_se), as.numeric(expected), tolerance = 1e-14)
  zero <- .plqr_range_fit(theta_meat = "jackknife", D = 2, density = 1)
  expect_identical(as.numeric(zero$theta_se), 0)
  expect_identical(as.numeric(zero$theta_vcov[[1L]]), 0)
  expect_warning(singular <- .plqr_range_fit(D = 2, density = 1, mD = 2),
    "Jacobian is singular")
  expect_true(all(is.na(singular$theta_se)))
  })
})

.plqr_product_range_fit <- function(theta_meat = "analytic", density = 1e-200,
                                    mD = 0, Y = rep(c(-1, 1), 4L)) {
  plqr_dnn(X = matrix(seq_len(8L)), x0 = 0, Y = Y,
    D = matrix(1, 8L), mD = matrix(mD, 8L), g = rep(0, 8L),
    theta_eval = 0, density = rep(density, 8L), s = 3L,
    theta_meat = theta_meat, return = character())
}

test_that("PLQR sandwich overflow preserves finite coefficients and unavailable uncertainty", {
  old <- options(RdistNN.use_compiled = FALSE)
  on.exit(options(old), add = TRUE)
  for (meat in c("analytic", "jackknife")) for (mD in c(0, 2)) {
    expect_warning(fit <- .plqr_product_range_fit(theta_meat = meat, mD = mD),
      "sandwich covariance.*outside the numeric range.*[Rr]escale")
    expect_identical(as.numeric(coef(fit)), 0)
    expect_true(all(is.finite(fit$moment)))
    expect_true(all(is.finite(fit$jacobian[[1L]])))
    expect_true(all(is.na(fit$theta_se)))
    expect_true(all(is.na(vcov(fit))))
    expect_true(all(is.na(fit$theta_omega[[1L]])))
    expect_equal(as.numeric(fit$jacobian[[1L]]) / 1e-200, 1 - mD, tolerance = 1e-14)
    # Positive scalar meat and a bread of magnitude 1e-200 imply a
    # coefficient variance above the largest double; it cannot be zero.
    local_meat <- if (meat == "analytic") fit$score_covariance[[1L]] / fit$neff else fit$moment_vcov[[1L]]
    expect_gt(as.numeric(local_meat), 0)
    log_variance <- log(as.numeric(local_meat)) - 2 * log(abs(as.numeric(fit$jacobian[[1L]])))
    expect_gt(log_variance, log(.Machine$double.xmax))
    restored <- unserialize(serialize(fit, NULL))
    expect_identical(coef(restored), coef(fit))
    expect_identical(vcov(restored), vcov(fit))
    expect_s3_class(summary(restored), "rdistnn_summary_plqr_result")
    expect_warning(pred <- predict(restored, x0 = 0, se.fit = TRUE),
      "sandwich covariance.*outside the numeric range.*[Rr]escale")
    expect_identical(as.numeric(pred$fit), 0)
    expect_true(all(is.na(pred$se.fit)))
  }
})

test_that("PLQR product guards preserve finite signed bread and genuine zero meat", {
  for (mD in c(0, 2)) {
    ordinary <- .plqr_product_range_fit(density = 1, mD = mD)
    expect_equal(as.numeric(vcov(ordinary)), .25 / as.numeric(ordinary$neff), tolerance = 1e-14)
    expect_true(all(is.finite(ordinary$theta_se)))
    zero <- .plqr_product_range_fit(theta_meat = "jackknife", mD = mD, Y = rep(1, 8L))
    expect_identical(as.numeric(vcov(zero)), 0)
    expect_identical(as.numeric(zero$theta_se), 0)
    expect_true(all(is.finite(zero$jacobian[[1L]])))
  }
})
