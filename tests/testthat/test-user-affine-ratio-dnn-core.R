.ratio_fixture <- function(...) {
  args <- list(X = matrix(1:4, ncol = 1L), x0 = 0,
               numerator = c(2, 1, 4, 3), denominator = c(1, -1, 2, 1),
               s = 2L, return = character())
  extra <- list(...)
  args[names(extra)] <- extra
  do.call(affine_ratio_dnn, args)
}

test_that("hand fixture returns paired ratio covariance and exact array identities", {
  fit <- .ratio_fixture(return = c("W", "signals", "deletions", "pseudo_values", "ess", "tuning"))
  expect_s3_class(fit, "affine_ratio_result")
  expect_equal(unname(fit$components), matrix(c(2, .5), 1L))
  expect_equal(fit$theta, 4)
  expect_equal(unname(fit$W[, 1L]), c(.5, 1 / 3, 1 / 6, 0))
  expect_equal(unname(fit$component_vcov[[1L]]), matrix(c(.5, .5, .5, .75), 2L))
  expect_equal(unname(fit$theta_vcov[[1L]]), matrix(34, 1L))
  expect_equal(unname(fit$theta_omega[[1L]]), matrix(612 / 7, 1L))
  expect_equal(fit$se, sqrt(34))
  expect_equal(fit$neff, 18 / 7)
  expect_identical(fit$ess, fit$neff)
  expect_identical(dim(fit$deletions), c(4L, 2L, 1L))
  expect_equal(as.numeric(fit$deletions[, 2L, 1L]), c(0, 4 / 3, 1 / 3, 1 / 3))
  expect_equal(as.numeric(fit$pseudo_values[, 1L, 1L]), c(2, 0, 3, 3))
  expect_equal(as.numeric(fit$pseudo_values[, 2L, 1L]), c(2, -2, 1, 1))
  expect_equal(4 * fit$pseudo_values[, 2L, 1L] - fit$pseudo_values[, 1L, 1L], c(6, -8, 1, 1))
  expect_equal(colMeans(fit$pseudo_values[, , 1L]), as.numeric(fit$components[1L, ]), ignore_attr = TRUE)
  expect_true(fit$status$wald_available)
})

test_that("subset enumeration matches full and retained arrays including ties", {
  for (n in 2:8) for (s in seq_len(n - 1L)) {
    X <- matrix(seq_len(n) %% 3L, ncol = 1L)
    signal <- cbind(seq_len(n)^2 / n, 1 + seq_len(n) / n)
    ref <- .clate_subset_reference(X, 0, signal, s)
    run <- function() affine_ratio_dnn(
      X = X, x0 = 0, numerator = signal[, 1L], denominator = signal[, 2L],
      s = s, return = c("deletions", "pseudo_values")
    )
    if (s == 1L) expect_warning(fit <- run(), "uniform") else fit <- run()
    expect_equal(as.numeric(fit$components), ref$components, tolerance = 1e-12)
    expect_equal(unname(fit$deletions[, , 1L]), unname(ref$loo), tolerance = 1e-12)
    expect_equal(unname(fit$component_vcov[[1L]]), unname(ref$vcov), tolerance = 1e-12)
  }
})

test_that("signed and zero denominators retain every target and component covariance", {
  fit <- affine_ratio_dnn(
    X = matrix(1:4, ncol = 1L), X0 = matrix(c(0, 3, 5), ncol = 1L),
    numerator = c(2, 1, 4, 3), denominator = c(0, 0, -1, 2), s = 3L
  )
  expect_equal(fit$status$target_id, 1:3)
  expect_equal(fit$status$ratio_status, c("zero_denominator", "ok", "ok"))
  expect_identical(fit$status$negative_denominator, c(FALSE, TRUE, FALSE))
  expect_true(is.na(fit$theta[1L]))
  expect_true(is.na(fit$se[1L]))
  expect_true(all(is.na(fit$theta_vcov[[1L]])))
  expect_true(all(is.finite(fit$component_vcov[[1L]])))
  expect_true(all(is.finite(fit$theta[-1L])))
  expect_identical(dim(fit$components), c(3L, 2L))
  expect_identical(dim(fit$x0), c(3L, 1L))
})

test_that("affine signal transformations preserve ratio and variance semantics", {
  base <- .ratio_fixture()
  for (a in c(-3, 0, 2)) {
    fit <- .ratio_fixture(numerator = a * c(2, 1, 4, 3) + 4 * c(1, -1, 2, 1))
    expect_equal(fit$theta, a * base$theta + 4)
    expect_equal(fit$theta_vcov[[1L]], a^2 * base$theta_vcov[[1L]], tolerance = 1e-12)
    if (a == 0) {
      expect_equal(fit$se, 0)
      expect_identical(fit$status$covariance_status, "zero_score_variance")
      expect_false(fit$status$wald_available)
    }
  }
  for (scale in c(-1e6, 1e-6, 1e6)) {
    fit <- .ratio_fixture(numerator = scale * c(2, 1, 4, 3),
                          denominator = scale * c(1, -1, 2, 1))
    expect_equal(fit$theta, base$theta, tolerance = 1e-12)
    expect_equal(fit$theta_vcov, base$theta_vcov, tolerance = 1e-10)
  }
})

test_that("singular paired covariance is valid with a deterministic first stage", {
  fit <- .ratio_fixture(denominator = rep(1, 4L))
  expect_true(fit$se > 0)
  expect_equal(unname(fit$component_vcov[[1L]][2L, ]), c(0, 0))
  expect_true(fit$status$wald_available)
})

test_that("matrix and data-frame interfaces preserve target dimensions and names", {
  d <- data.frame(x = 1:4, y = c(2, 1, 4, 3), d = c(1, -1, 2, 1))
  f <- affine_ratio_dnn(data = d, x_cols = "x", numerator_col = "y",
                        denominator_col = "d", X0 = matrix(c(0, 5), ncol = 1L), s = 2L)
  g <- affine_ratio_dnn(X = matrix(d$x, ncol = 1L), X0 = matrix(c(0, 5), ncol = 1L),
                        numerator = d$y, denominator = d$d, s = 2L)
  expect_equal(f$theta, g$theta)
  expect_equal(f$theta_vcov, g$theta_vcov)
  expect_identical(colnames(f$x0), "x")
  expect_length(f$component_vcov, 2L)
  no_cov <- .ratio_fixture(se = FALSE)
  expect_null(no_cov$se)
  expect_null(no_cov$theta_vcov)
  expect_null(no_cov$component_vcov)
  expect_identical(no_cov$status$covariance_status, "not_computed")
  expect_false("W" %in% names(no_cov))
  expect_equal(.ratio_fixture(se = TRUE, jackknife = "none")$theta, no_cov$theta)
  expect_null(.ratio_fixture(se = TRUE, jackknife = "none")$se)
})

test_that("one-observation point estimates and same-order boundary deletions work", {
  fit <- affine_ratio_dnn(X = matrix(2, 1L), x0 = 0,
                          numerator = 6, denominator = 2, s = 1L, se = FALSE)
  expect_equal(fit$theta, 3)
  expect_equal(fit$neff, 1)
  expect_identical(dim(fit$x0), c(1L, 1L))
  expect_equal(.ratio_fixture(s = 4L, se = FALSE)$theta, 2)
  expect_warning(fit2 <- affine_ratio_dnn(X = matrix(1:2, ncol = 1L), x0 = 0,
                           numerator = c(2, 4), denominator = c(1, 1), s = 1L), "uniform")
  expect_equal(fit2$theta, 3)
  expect_equal(fit2$se, 1)
})
