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
