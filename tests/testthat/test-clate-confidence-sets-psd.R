test_that("ratio covariance checks detect negative determinants hidden by correlation rounding", {
  components <- matrix(c(2, 1), 1L)
  # The represented sqrt(2) is above the exact square root. Its squared
  # covariance therefore exceeds 1*2, although sqrt(2)/sqrt(2) rounds to 1.
  # The second example has exact determinant -2^-82, while both rounded
  # products in the ordinary determinant are equal to 1-2^-40.
  cases <- list(matrix(c(1, sqrt(2), sqrt(2), 2), 2L),
    matrix(c(1 + 2^-20, 1 - 2^-41, 1 - 2^-41, 1 - 2^-20), 2L))
  for (V in cases) {
    expect_error(ratio_confidence_set(components, list(V), tol = 0), "positive semidefinite")
    out <- ratio_confidence_set(components, list(V))[[1L]]
    expect_identical(out$type, "unresolved")
    expect_identical(out$covariance_status, "near_psd_boundary")
    expect_identical(out$numerical_status, "near_psd_boundary")
    expect_identical(nrow(out$intervals), 0L)
  }
})

test_that("exactly singular and positive covariance survives rounded correlation boundaries", {
  components <- matrix(c(2, 1), 1L)
  # Integer products establish determinant zero without any square roots.
  # For a=3 the quotient of rounded square roots can exceed one.
  for (a in c(1, 2, 3, 5, 7, 11, 31)) {
    V <- matrix(c(a, 2 * a, 2 * a, 4 * a), 2L)
    out <- ratio_confidence_set(components, list(V), tol = 0)[[1L]]
    expect_identical(out$covariance_status, "ok")
    expect_false(out$type == "unresolved")
  }
  V <- matrix(c(1 + 2^-20, 1 - 2^-41 - 2^-53,
    1 - 2^-41 - 2^-53, 1 - 2^-20), 2L)
  out <- ratio_confidence_set(components, list(V), tol = 0)[[1L]]
  expect_identical(out$covariance_status, "ok")
  expect_false(out$type == "unresolved")
})

test_that("ratio PSD classification is preserved under exact binary coordinate scaling", {
  components <- matrix(c(2, 1), 1L)
  invalid <- matrix(c(1, sqrt(2), sqrt(2), 2), 2L)
  singular <- matrix(c(3, 6, 6, 12), 2L)
  for (powers in list(c(-500, 500), c(500, -500), c(-500, -500), c(500, 500))) {
    scale <- 2^powers
    delta <- sweep(components, 2L, scale, "*")
    V <- invalid * outer(scale, scale)
    expect_error(ratio_confidence_set(delta, list(V), tol = 0), "positive semidefinite")
    out <- ratio_confidence_set(delta, list(V))[[1L]]
    expect_identical(out$numerical_status, "near_psd_boundary")
    exact <- ratio_confidence_set(delta, list(singular * outer(scale, scale)), tol = 0)[[1L]]
    expect_identical(exact$covariance_status, "ok")
    expect_false(exact$type == "unresolved")
  }
})
