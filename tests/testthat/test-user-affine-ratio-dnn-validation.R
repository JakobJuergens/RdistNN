test_that("unsupported controls fail before distance evaluation", {
  args <- list(X = matrix(1:4, ncol = 1L), x0 = 0,
               numerator = 1:4, denominator = rep(1, 4L), s = 2L)
  distance <- function(...) stop("distance should not run")
  for (item in list(list(k = 2), list(method = "asymp"), list(s = 5),
                    list(s = 4), list(s = 0), list(s = 1.5),
                    list(s = 1e100), list(denominator_tol = -1))) {
    a <- args
    a[names(item)] <- item
    a$distance <- distance
    expect_error(do.call(affine_ratio_dnn, a), paste(names(item), collapse = "|"))
  }
})

test_that("mixed and incomplete interfaces produce argument-specific errors", {
  d <- data.frame(x = 1:4, y = 1:4, d = rep(1, 4L))
  expect_error(affine_ratio_dnn(data = d, X = matrix(d$x), x_cols = "x",
                               numerator_col = "y", denominator_col = "d", x0 = 0), "data/X")
  expect_error(affine_ratio_dnn(data = d, x_cols = "x", numerator = d$y,
                               denominator_col = "d", x0 = 0), "data")
  expect_error(affine_ratio_dnn(data = d, x_cols = "missing",
                               numerator_col = "y", denominator_col = "d", x0 = 0), "missing")
  expect_error(affine_ratio_dnn(X = matrix(d$x), x0 = 0, numerator = d$y), "denominator")
  expect_error(affine_ratio_dnn(X = matrix(d$x), x0 = 0, X0 = matrix(0),
                               numerator = d$y, denominator = d$d), "X0/x0")
  expect_error(affine_ratio_dnn(X = matrix(d$x), numerator = d$y,
                               denominator = d$d), "X0/x0")
})

test_that("signals targets IDs and optional-return controls are validated", {
  args <- list(X = matrix(1:4, ncol = 1L), x0 = 0,
               numerator = 1:4, denominator = rep(1, 4L), s = 2L)
  for (nm in c("numerator", "denominator")) for (bad in list(1:3, matrix(1:4),
      c(1, 2, 3, NA_real_), c(1, 2, 3, NaN), c(1, 2, 3, Inf))) {
    a <- args
    a[[nm]] <- bad
    expect_error(do.call(affine_ratio_dnn, a), nm)
  }
  expect_error(do.call(affine_ratio_dnn, c(args, list(row_ids = c(1, 1, 2, 3)))), "row_ids")
  expect_warning(f <- do.call(affine_ratio_dnn, c(args, list(return = "unknown"))), "unknown")
  expect_false("unknown" %in% names(f))
  args$x0 <- c(0, 1)
  expect_error(do.call(affine_ratio_dnn, args), "x0")
})

test_that("numerical cancellation tolerance never changes a representable ratio", {
  args <- list(X = matrix(1:4, ncol = 1L), x0 = 0,
               numerator = 1:4, denominator = c(1, -3 + 1e-9, 0, 0), s = 3L)
  a <- do.call(affine_ratio_dnn, c(args, list(denominator_tol = 1e-8)))
  b <- do.call(affine_ratio_dnn, c(args, list(denominator_tol = 0)))
  expect_true(a$status$near_zero_denominator)
  expect_false(b$status$near_zero_denominator)
  expect_identical(a$theta, b$theta)
  expect_identical(a$theta_vcov, b$theta_vcov)
})
