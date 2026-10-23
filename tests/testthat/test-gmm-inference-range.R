.gmm_range_model <- function(scale) {
  force(scale)
  new_gmm_model("information_range",
    setup = function(sample, n, ...) list(y = sample$y),
    moment_obs = function(sample, theta, ctx) {
      matrix(scale * theta[1L] - ctx$y, ncol = 1L)
    },
    jacobian_obs = function(sample, theta, ctx) {
      matrix(scale, length(ctx$y), 1L)
    }, coef_names = "location")
}

.gmm_range_fit <- function(scale) {
  gmm_dnn(X = matrix(seq_len(8L)), x0 = 4,
    model_data = data.frame(y = rep(c(-1, 1), 4L)),
    model = .gmm_range_model(scale), start = 0, s = 1L,
    se = TRUE, omega_ridge = 0, return = c("W", "jacobian"))
}

.gmm_range_capture <- function(expr) {
  seen <- character()
  value <- withCallingHandlers(expr, warning = function(w) {
    seen <<- c(seen, conditionMessage(w))
    invokeRestart("muffleWarning")
  })
  list(value = value, warnings = seen)
}

test_that("nonfinite GMM information cannot become zero coefficient uncertainty", {
  .with_distance_tie_warning({
  captured <- .gmm_range_capture(.gmm_range_fit(1e155))
  fit <- captured$value
  expect_identical(as.numeric(coef(fit)), 0)
  expect_true(all(is.finite(fit$jacobian[[1L]])))
  expect_true(all(is.na(fit$se)))
  expect_null(vcov(fit))
  expect_true(any(grepl("information.*outside the numeric range.*[Rr]escale",
    captured$warnings)))
  # The raw information overflows, while the independent coefficient
  # covariance remains strictly positive and representable.
  expected_vcov <- (1 / 1e155)^2 / 8
  expect_gt(expected_vcov, 0)
  expect_true(is.finite(expected_vcov))
  })
})

test_that("GMM safe solve rejects nonfinite operands before LAPACK", {
  for (invalid in c(Inf, -Inf, NA_real_, NaN)) {
    expect_null(RdistNN:::.gmm_safe_solve(matrix(invalid)))
    expect_null(RdistNN:::.gmm_safe_solve(diag(2L), c(0, invalid)))
  }
  expect_identical(RdistNN:::.gmm_safe_solve(diag(2L)), diag(2L))
  expect_identical(RdistNN:::.gmm_safe_solve(diag(2L), c(0, 0)), c(0, 0))
})

test_that("finite GMM information preserves ordinary and singular controls", {
  .with_distance_tie_warning({
  ordinary <- .gmm_range_capture(.gmm_range_fit(1))
  expect_identical(as.numeric(coef(ordinary$value)), 0)
  expect_equal(as.numeric(vcov(ordinary$value)), 1 / 8, tolerance = 1e-14)
  expect_equal(as.numeric(ordinary$value$se), 1 / sqrt(8), tolerance = 1e-14)
  expect_false(any(grepl("outside the numeric range", ordinary$warnings)))

  singular <- .gmm_range_capture(.gmm_range_fit(0))
  expect_identical(as.numeric(coef(singular$value)), 0)
  expect_true(all(is.na(singular$value$se)))
  expect_null(vcov(singular$value))
  expect_false(any(grepl("outside the numeric range", singular$warnings)))
  })
})
