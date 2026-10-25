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

test_that("GMM moment rescaling preserves covariance and conditioning diagnostics", {
  .with_distance_tie_warning({
    for (scale in c(1, 1e-10, 1e-50)) {
      fit <- .gmm_range_capture(gmm_dnn(
        X = matrix(seq_len(8L)), x0 = 0,
        model_data = data.frame(y = scale * rep(c(-1, 1), 4L)),
        model = .gmm_range_model(scale), start = 0, s = 1L,
        se = TRUE, omega_ridge = 0, jacobian = "analytic", return = "omega"
      ))$value
      expect_identical(as.numeric(coef(fit)), 0)
      expect_equal(as.numeric(vcov(fit)), 1 / 8, tolerance = 1e-14)
      expect_equal(as.numeric(fit$omega[[1L]]) / scale^2, 1, tolerance = 1e-14)
      expect_equal(fit$convergence$omega_condition, 1, tolerance = 1e-14)
      restored <- unserialize(serialize(fit, NULL))
      expect_equal(summary(restored)$convergence$omega_condition, 1, tolerance = 1e-14)
    }
  })
})

test_that("GMM spectral conditioning handles singular and large finite matrices", {
  expect_identical(RdistNN:::.gmm_condition_number(matrix(0)), Inf)
  expect_identical(RdistNN:::.gmm_condition_number(diag(c(1, 0))), Inf)
  # Both entries and the true condition are finite, although the largest
  # unscaled eigenvalue 1.8e308 is outside double precision's range.
  dense <- matrix(c(1e308, 8e307, 8e307, 1e308), nrow = 2L)
  expect_equal(RdistNN:::.gmm_condition_number(dense), 9, tolerance = 1e-14)
})

test_that("finite norm-overflow inverses retain signed and nonsymmetric solutions", {
  cases <- list(
    list(matrix(c(1, .8, .8, 1), 2L), matrix(c(25, -20, -20, 25), 2L) / 9),
    list(-matrix(c(1, .8, .8, 1), 2L), -matrix(c(25, -20, -20, 25), 2L) / 9),
    list(matrix(c(1, -.8, .8, 1), 2L), matrix(c(25, 20, -20, 25), 2L) / 41)
  )
  for (safe in list(RdistNN:::.gmm_safe_solve, RdistNN:::.plqr_safe_solve)) {
    for (case in cases) {
      inverse <- safe(1e308 * case[[1L]])
      expect_true(all(is.finite(inverse)))
      expect_equal(inverse / 1e-308, case[[2L]], tolerance = 1e-13)
    }
    expect_null(safe(matrix(1e308, 2L, 2L)))
    expect_null(safe(2^1023 * matrix(c(1, 1, 1, 1 + .Machine$double.eps), 2L)))
  }
})

test_that("exceptional GMM right-hand sides preserve finite scales", {
  A <- 1e308 * matrix(c(1, .8, .8, 1), 2L)
  expect_equal(RdistNN:::.gmm_safe_solve(A, c(1e308, -1e308)), c(5, -5),
    tolerance = 1e-13)
  rhs <- cbind(c(1e308, -1e308), c(1e-3, -1e-3))
  actual <- RdistNN:::.gmm_safe_solve(A, rhs)
  expect_equal(actual[, 1L], c(5, -5), tolerance = 1e-13)
  expect_equal(actual[, 2L] / 1e-311, c(5, -5), tolerance = 1e-12)
  delta <- 8 * .Machine$double.eps
  A <- 2^1023 * matrix(c(1, 1, 1, 1 + delta), 2L)
  expected <- (1e-20 / (2^1023 * delta)) * c(2 + delta, -2)
  expect_true(all(expected != 0))
  expect_equal(RdistNN:::.gmm_safe_solve(A, c(1e-20, -1e-20)) / 1e-313,
    expected / 1e-313, tolerance = 1e-9)
  mixed <- diag(c(1e308, 1e308, 1e294))
  mixed[1L, 2L] <- mixed[2L, 1L] <- 8e307
  expect_null(RdistNN:::.gmm_safe_solve(mixed, c(1e308, -1e308, 1e-20)))
})

test_that("normalized inverse branches preserve unexpected solver identities", {
  original <- structure(list(message = "opaque normalized solver failure", call = NULL),
    class = c("simulated_resource_error", "error", "condition"))
  local_mocked_bindings(solve = function(...) stop(original), .package = "base")
  A <- matrix(c(1e308, 8e307, 8e307, 1e308), 2L)
  for (safe in list(RdistNN:::.gmm_safe_solve, RdistNN:::.plqr_safe_solve)) {
    expect_identical(tryCatch(safe(A), error = identity), original)
  }
})

test_that("normalized inverse branches preserve unexpected condition-number identities", {
  original <- structure(list(message = "opaque normalized condition failure", call = NULL),
    class = c("simulated_resource_error", "error", "condition"))
  local_mocked_bindings(rcond = function(...) stop(original), .package = "base")
  A <- matrix(c(1e308, 8e307, 8e307, 1e308), 2L)
  for (safe in list(RdistNN:::.gmm_safe_solve, RdistNN:::.plqr_safe_solve)) {
    expect_identical(tryCatch(safe(A), error = identity), original)
  }
})

test_that("public GMM inference inverts a finite dense norm-overflow covariance", {
  .with_localization_policy_warnings({
    L <- 1e154 * rbind(c(sqrt(.9), sqrt(.1)), c(sqrt(.9), -sqrt(.1)))
    Y <- rbind(c(1, 1), c(-1, -1), c(1, -1), c(-1, 1))
    model <- new_gmm_model("dense_covariance",
      setup = function(sample, n, ...) list(y = Y, L = L),
      moment_obs = function(sample, theta, ctx) {
        (matrix(theta, nrow(ctx$y), 2L, byrow = TRUE) - ctx$y) %*% t(ctx$L)
      },
      jacobian_obs = function(sample, theta, ctx) {
        out <- array(0, c(nrow(ctx$y), 2L, 2L))
        for (j in 1:2) out[, , j] <- matrix(rep(ctx$L[, j], each = nrow(ctx$y)),
          nrow(ctx$y))
        out
      }, coef_names = c("one", "two"))
    fit <- gmm_dnn(X = matrix(1:4), x0 = 0, model = model, start = c(0, 0),
      s = 1L, omega_ridge = 0, jacobian = "analytic", return = c("omega", "jacobian"))
    expect_identical(as.numeric(coef(fit)), c(0, 0))
    expect_true(all(is.finite(fit$omega[[1L]])))
    expect_equal(unname(fit$omega[[1L]]) / 1e308,
      matrix(c(1, .8, .8, 1), 2L), tolerance = 1e-13)
    expect_equal(unname(vcov(fit)), diag(.25, 2L), tolerance = 1e-12)
    expect_equal(fit$convergence$omega_condition, 9, tolerance = 1e-13)
    restored <- unserialize(serialize(fit, NULL))
    expect_identical(vcov(restored), vcov(fit))
    expect_s3_class(summary(restored), "rdistnn_summary_gmm_result")
  })
})
