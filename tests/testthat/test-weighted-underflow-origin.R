test_that("finite-support deletions preserve zero survivors and constants under tiny weights", {
  old <- options(RdistNN.use_compiled = TRUE)
  on.exit(options(old), add = TRUE)
  eta <- .Machine$double.xmin * .Machine$double.eps
  expect_gt(eta, 0)
  ranks <- matrix(1:3, nrow = 1L)
  signal <- cbind(c(1.5, 0, 0), rep(1.5, 3L))
  expected <- cbind(c(0, 1, 1), rep(2, 3L))
  calls <- 0L
  native <- RdistNN:::rustat_jackknife_fast_rank
  testthat::local_mocked_bindings(rustat_jackknife_fast_rank = function(...) {
    calls <<- calls + 1L
    native(...)
  }, .package = "RdistNN")
  for (compiled in c(FALSE, TRUE)) {
    options(RdistNN.use_compiled = compiled)
    expect_identical(is.null(RdistNN:::.rdistnn_compiled_fun("rustat_jackknife_fast_rank")), !compiled)
    fast <- jackknife_fast_rank(ranks, c(eta, eta), signal, normalize = TRUE)
    slow <- jackknife_slow_rank(ranks, c(eta, eta), signal, normalize = TRUE)
    expect_equal(fast$theta0, c(1, 2), tolerance = 0)
    expect_equal(fast$theta_loo, expected, tolerance = 0)
    expect_equal(slow$theta_loo, expected, tolerance = 0)
    expect_identical(RdistNN:::.covariance_matrix(fast$theta_loo[, 2L, drop = FALSE], 2 / 3), matrix(0))
    expect_identical(calls, as.integer(compiled))
  }
})

test_that("nonzero subnormal signed mass cannot create constant-signal variation", {
  old <- options(RdistNN.use_compiled = TRUE)
  on.exit(options(old), add = TRUE)
  eta <- .Machine$double.xmin * .Machine$double.eps
  w <- 1.25 * .Machine$double.xmin
  table <- c(w, -(w - eta))
  signal <- matrix(rep(1.5 + .Machine$double.eps, 3L))
  expected <- matrix(1, 3L, 1L)
  for (compiled in c(FALSE, TRUE)) {
    options(RdistNN.use_compiled = compiled)
    fast <- jackknife_fast_rank(matrix(1:3, nrow = 1L), table, signal, normalize = TRUE)
    slow <- jackknife_slow_rank(matrix(1:3, nrow = 1L), table, signal, normalize = TRUE)
    expect_equal(fast$theta0, 1, tolerance = 0)
    expect_equal(fast$theta_loo, expected, tolerance = 0)
    expect_equal(slow$theta_loo, expected, tolerance = 0)
    expect_identical(RdistNN:::.covariance_matrix(fast$theta_loo, 2 / 3), matrix(0))
    weight_fun <- function(i) matrix(table, ncol = 1L)
    raw <- RdistNN:::.jackknife_recomputed_centered(signal, weight_fun, restore = TRUE)
    expect_equal(as.numeric(raw) / eta, rep(1, 3L), tolerance = 0)
  }
})

test_that("full-rank signed and unnormalized tiny weights use original rounded products", {
  old <- options(RdistNN.use_compiled = TRUE)
  on.exit(options(old), add = TRUE)
  eta <- .Machine$double.xmin * .Machine$double.eps
  ranks <- cbind(1:3, 3:1)
  signal <- cbind(c(1.5, .5, 0), rep(1.5, 3L))
  for (table in list(c(eta, eta), c(eta, -eta), c(2 * eta, -eta))) {
    expected <- array(0, c(3L, 2L, 2L))
    for (target in 1:2) for (i in 1:3) {
      keep <- (1:3) != i
      reduced <- ranks[keep, target] - (ranks[keep, target] > ranks[i, target])
      for (component in 1:2) expected[i, component, target] <- sum(table[reduced] * signal[keep, component])
    }
    for (compiled in c(FALSE, TRUE)) {
      options(RdistNN.use_compiled = compiled)
      expect_identical(is.null(RdistNN:::.rdistnn_compiled_fun("rustat_jackknife_rank_full")), !compiled)
      actual <- jackknife_rank_full(ranks, signal, c(table, 0), table)
      expect_identical(dim(actual$theta_loo), c(3L, 2L, 2L))
      expect_equal(actual$theta_loo / eta, expected / eta, tolerance = 0)
    }
  }
})

test_that("recomputed and grouped paths apply the weighted amplitude origin guard", {
  old <- options(RdistNN.use_compiled = TRUE)
  on.exit(options(old), add = TRUE)
  eta <- .Machine$double.xmin * .Machine$double.eps
  signal <- cbind(c(1.5, 0, 0), rep(1.5, 3L))
  expected <- cbind(c(0, 2, 2), rep(4, 3L))
  weight_fun <- function(i) matrix(c(eta, eta), ncol = 1L)
  for (compiled in c(FALSE, TRUE)) {
    options(RdistNN.use_compiled = compiled)
    reference <- RdistNN:::.jackknife_recomputed_centered(signal, weight_fun, restore = TRUE)
    expect_equal(matrix(reference, 3L, 2L) / eta, expected, tolerance = 0)
    for (groups in list(NULL, c(1L, 1L, 2L))) {
      shifted <- RdistNN:::.dnn_loo_centered(signal, matrix(1:3), c(eta, eta), groups = groups)
      raw <- RdistNN:::.dnn_restore_loo(shifted$loo, shifted$offset)
      expect_equal(matrix(raw, 3L, 2L) / eta, expected, tolerance = 0)
    }
  }
})

test_that("weighted origins preserve ordinary tiny tails and exact zero mass", {
  eta <- .Machine$double.xmin * .Machine$double.eps
  level <- RdistNN:::.numeric_signal_level
  expect_equal(level(matrix(c(1.5, 0, 0)), weights = c(eta, eta)), 0, tolerance = 0)
  expect_equal(level(matrix(c(1.5, 0, 0)), weights = c(1, eta)), 1.5, tolerance = 0)
  expect_equal(level(matrix(rep(2^1000, 3L)), weights = c(eta, eta)), 0, tolerance = 0)
  expect_equal(level(matrix(rep(1.5, 3L))), 1.5, tolerance = 0)
  expect_equal(level(matrix(rep(1.5, 3L)), weights = c(1, -1)), 1.5, tolerance = 0)
})

test_that("subnormal table mass preserves amplified rounded block coefficients", {
  old <- options(RdistNN.use_compiled = TRUE)
  on.exit(options(old), add = TRUE)
  eta <- .Machine$double.xmin * .Machine$double.eps
  signal <- matrix(rep(2^1000, 4L))
  for (units in c(1, 4)) {
    table <- c(1, -1, units * eta)
    coefficient <- if (units == 1) 0 else eta
    expected <- rep(3 * (coefficient * 2^1000), 4L)
    for (compiled in c(FALSE, TRUE)) {
      options(RdistNN.use_compiled = compiled)
      shifted <- RdistNN:::.dnn_loo_centered(signal, matrix(1:4), table, groups = rep(1L, 4L))
      raw <- RdistNN:::.dnn_restore_loo(shifted$loo, shifted$offset)
      expect_equal(as.numeric(raw) / 2^-74, expected / 2^-74, tolerance = 0)
    }
  }
})
