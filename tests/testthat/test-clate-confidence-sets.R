test_that("quadratic inversion covers every geometry and both half-line directions", {
  cases <- list(
    list(c(1, 0, -1), "bounded", matrix(c(-1, 1), 1L)),
    list(c(1, -2, 1), "singleton", matrix(c(1, 1), 1L)),
    list(c(.25, -.5, .25), "singleton", matrix(c(1, 1), 1L)),
    list(c(1, 0, 1), "empty", matrix(numeric(), 0L, 2L)),
    list(c(-1, 0, 1), "two_rays", rbind(c(-Inf, -1), c(1, Inf))),
    list(c(-1, 0, 0), "all_real", matrix(c(-Inf, Inf), 1L)),
    list(c(-1, 0, -1), "all_real", matrix(c(-Inf, Inf), 1L)),
    list(c(0, 2, -4), "half_line", matrix(c(-Inf, 2), 1L)),
    list(c(0, -2, 4), "half_line", matrix(c(2, Inf), 1L)),
    list(c(0, 0, -1), "all_real", matrix(c(-Inf, Inf), 1L)),
    list(c(0, 0, 1), "empty", matrix(numeric(), 0L, 2L)),
    list(c(0, 0, 0), "all_real", matrix(c(-Inf, Inf), 1L))
  )
  for (x in cases) {
    fit <- do.call(.ratio_quadratic_set, as.list(c(x[[1L]], tol = 64 * .Machine$double.eps)))
    expect_identical(fit$type, x[[2L]])
    expect_equal(unname(fit$intervals), x[[3L]], tolerance = 1e-12)
  }
})

test_that("paired confidence sets reduce to Wald geometry with deterministic denominator", {
  comp <- matrix(c(2, 1), 1L)
  V <- list(diag(c(.5, 0)))
  set <- ratio_confidence_set(comp, V)
  critical <- stats::qnorm(.025, lower.tail = FALSE)
  expect_s3_class(set, "ratio_confidence_sets")
  expect_identical(set[[1L]]$type, "bounded")
  expect_equal(as.numeric(set[[1L]]$intervals), 2 + c(-1, 1) * critical * sqrt(.5))
  expect_output(ret <- withVisible(print(set)), "bounded")
  expect_false(ret$visible)
  expect_identical(ret$value, set)
  level <- 1 - .Machine$double.eps / 2
  expect_true(all(is.finite(ratio_confidence_set(comp, V, level = level)[[1L]]$intervals)))
})

test_that("zero denominators and singular covariances retain legitimate set-valued output", {
  comp <- rbind(c(0, 0), c(3, 0), c(0, .1), c(3, 0), c(2, 1))
  cov <- list(matrix(0, 2L, 2L), matrix(0, 2L, 2L), diag(2L), diag(2L), matrix(0, 2L, 2L))
  sets <- ratio_confidence_set(comp, cov)
  expect_equal(vapply(sets, function(x) x$type, character(1L)),
               c("all_real", "empty", "all_real", "two_rays", "singleton"))
  expect_equal(sets[[5L]]$intervals[1L, ], c(lower = 2, upper = 2))
  expect_equal(nrow(sets[[2L]]$intervals), 0L)
  expect_equal(nrow(sets[[4L]]$intervals), 2L)
  path <- tempfile()
  on.exit(unlink(path), add = TRUE)
  saveRDS(sets, path)
  expect_equal(readRDS(path), sets)
  expect_output(print(sets), "Inf")
})

test_that("resolved sets agree with the original moment inequality and contain the fitted ratio", {
  set.seed(803)
  comp <- matrix(stats::rnorm(24L), 12L, 2L)
  comp[, 2L] <- comp[, 2L] + 1
  V <- lapply(seq_len(12L), function(i) {
    B <- matrix(stats::rnorm(4L), 2L)
    crossprod(B) / 10 + diag(.02, 2L)
  })
  sets <- ratio_confidence_set(comp, V)
  contains <- function(set, value) {
    any(vapply(seq_len(nrow(set$intervals)), function(i) {
      value >= set$intervals[i, 1L] && value <= set$intervals[i, 2L]
    }, logical(1L)))
  }
  z <- stats::qnorm(.025, lower.tail = FALSE)
  for (j in seq_len(12L)) {
    expect_false(sets[[j]]$type == "unresolved")
    expect_true(contains(sets[[j]], comp[j, 1L] / comp[j, 2L]))
    for (theta in c(-10, -2, 0, 1, 5, 20)) {
      lhs <- (comp[j, 1L] - theta * comp[j, 2L])^2
      rhs <- z^2 * (V[[j]][1L, 1L] - 2 * theta * V[[j]][1L, 2L] + theta^2 * V[[j]][2L, 2L])
      expect_identical(contains(sets[[j]], theta), lhs <= rhs)
    }
  }
})

test_that("confidence inversion handles scale changes and reports unresolved boundaries", {
  comp <- matrix(c(2, 1), 1L)
  cov <- list(diag(c(.5, .1)))
  original <- ratio_confidence_set(comp, cov)[[1L]]$intervals
  for (b in c(-3, .1, 1e100)) {
    transform <- diag(c(b, 1))
    scaled <- ratio_confidence_set(comp %*% transform, list(transform %*% cov[[1L]] %*% transform))[[1L]]
    expect_equal(as.numeric(scaled$intervals), sort(b * as.numeric(original)), tolerance = 1e-11)
  }
  near <- .ratio_quadratic_set(1, -2, 1 + 1e-15, tol = 1e-12)
  expect_identical(near$type, "unresolved")
  expect_identical(near$numerical_status, "near_boundary")
  small_root <- .ratio_quadratic_set(1, 1e12, 1, tol = 0)
  expect_equal(unname(small_root$intervals[1L, 2L]), -1e-12, tolerance = 1e-12)
})

test_that("confidence sets reject incoherent input without silently repairing covariance", {
  good <- matrix(c(2, 1), 1L)
  for (bad in list(matrix(1, 1L, 1L), matrix(c(NA, 1), 1L), matrix(c(Inf, 1), 1L))) {
    expect_error(ratio_confidence_set(bad, list(diag(2L))), "components")
  }
  expect_error(ratio_confidence_set(good, diag(2L)), "component_vcov")
  expect_error(ratio_confidence_set(good, list(diag(3L))), "component_vcov")
  expect_error(ratio_confidence_set(good, list(matrix(c(1, 3, 0, 1), 2L))), "symmetric")
  expect_error(ratio_confidence_set(good, list(matrix(c(1, 2, 2, 1), 2L))), "positive semidefinite")
  expect_error(ratio_confidence_set(good, list(diag(c(-.1, 1)))), "positive semidefinite")
  for (level in list(0, 1, NA_real_, "x")) {
    expect_error(ratio_confidence_set(good, list(diag(2L)), level = level), "level")
  }
  expect_error(ratio_confidence_set(good, list(diag(2L)), tol = -1), "tol")
})

test_that("power-of-two coordinate scaling preserves extreme but representable sets", {
  comp <- matrix(c(2, 1), 1L)
  V <- diag(c(.5, .1))
  reference <- ratio_confidence_set(comp, list(V))[[1L]]
  for (exponent in c(-500, 500)) {
    for (coordinate in 1:2) {
      scale <- rep(1, 2L)
      scale[coordinate] <- 2^exponent
      out <- ratio_confidence_set(
        sweep(comp, 2L, scale, "*"), list(V * outer(scale, scale))
      )[[1L]]
      expect_identical(out$type, reference$type)
      expect_equal(as.numeric(out$intervals) / (scale[1L] / scale[2L]),
                   as.numeric(reference$intervals), tolerance = 1e-12)
    }
  }
  out <- ratio_confidence_set(comp, list(V))[[1L]]
  z2 <- stats::qnorm(.025, lower.tail = FALSE)^2
  expect_equal(out$quadratic, c(a = 1 - z2 * .1, b = -4, c = 4 - z2 * .5))
  expect_identical(ratio_confidence_set(matrix(c(pi, 3), 1L),
                                      list(matrix(0, 2L, 2L)))[[1L]]$type, "singleton")
})

test_that("unresolvable ranges and nearly indefinite covariance have explicit statuses", {
  comp <- matrix(c(2, 1), 1L)
  almost <- matrix(c(1, 1 + 2 * .Machine$double.eps,
                     1 + 2 * .Machine$double.eps, 1), 2L)
  out <- ratio_confidence_set(comp, list(almost))[[1L]]
  expect_identical(out$type, "unresolved")
  expect_identical(out$numerical_status, "near_psd_boundary")
  extreme <- ratio_confidence_set(matrix(c(1e300, 1e-300), 1L),
                                  list(matrix(0, 2L, 2L)))[[1L]]
  expect_identical(extreme$type, "unresolved")
  expect_identical(extreme$numerical_status, "endpoint_range")
  expect_error(ratio_confidence_set(matrix(numeric(), 0L, 2L), list()), "components")
  expect_error(ratio_confidence_set(comp, list(matrix(c(0, 1e-300, 1e-300, 1), 2L))),
               "positive semidefinite")
})

test_that("subnormal covariance is preserved and underflowed discriminants are labelled", {
  comp <- matrix(c(1, 2), 1L)
  covariance <- matrix(1, 2L, 2L)
  ordinary <- ratio_confidence_set(comp, list(covariance))[[1L]]
  tiny <- ratio_confidence_set(comp * 2^-537, list(covariance * 2^-1074))[[1L]]
  expect_identical(tiny$type, ordinary$type)
  expect_equal(tiny$intervals, ordinary$intervals, tolerance = 1e-12)
  expect_true(is.na(tiny$discriminant))
  expect_true(is.finite(tiny$scaled_discriminant))
})
