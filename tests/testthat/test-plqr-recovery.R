.plqr_recovery_fit <- function() {
  plqr_dnn(X = matrix(seq_len(8L)), x0 = 1.25, Y = rep(1, 8L),
    D = matrix(seq_len(8L)), mD = matrix(0, 8L), g = rep(0, 8L),
    theta_eval = 0, density = rep(1, 8L), s = 3L)
}

test_that("PLQR inference preserves unexpected solver failure identities", {
  original <- structure(list(message = "opaque solver failure", call = NULL),
    class = c("simulated_resource_error", "error", "condition"))
  local_mocked_bindings(solve = function(...) stop(original), .package = "base")
  observed <- tryCatch(.plqr_recovery_fit(), error = identity)
  expect_identical(observed, original)
})

test_that("PLQR does not relabel caller-promoted solver warnings as singularity", {
  old <- options(warn = 2)
  on.exit(options(old), add = TRUE)
  local_mocked_bindings(solve = function(...) warning("solver warning sentinel"), .package = "base")
  observed <- tryCatch(.plqr_recovery_fit(), error = identity)
  expect_s3_class(observed, "error")
  expect_match(conditionMessage(observed), "solver warning sentinel")
})

test_that("PLQR numerical preflight preserves solve's existing singular cutoff", {
  safe <- RdistNN:::.plqr_safe_solve
  eps <- .Machine$double.eps
  for (a in c(0, eps / 2, eps, 2 * eps, .5)) {
    A <- diag(c(1, a))
    expect_equal(safe(A), tryCatch(solve(A), error = function(e) NULL))
  }
  for (value in c(NA_real_, NaN, Inf, -Inf)) expect_null(safe(matrix(value)))
})
