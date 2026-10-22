test_that("fast identity comparisons preserve exact values and structure", {
  expect_true(exists(".expect_identical_fast", mode = "function"))
  if (!exists(".expect_identical_fast", mode = "function")) return()

  x <- list(values = matrix(c(1, NA_real_, Inf, 2^-1074), 2L),
    status = c("ok", "unavailable"))
  returned <- withVisible(.expect_identical_fast(x, x))
  expect_identical(returned$value, x)
  expect_false(returned$visible)
  changed <- x
  changed$values[1L] <- 1 + .Machine$double.eps
  expect_failure(.expect_identical_fast(changed, x), "identical")
  changed <- x
  dim(changed$values) <- c(1L, 4L)
  expect_failure(.expect_identical_fast(changed, x), "identical")
  expect_failure(.expect_identical_fast(1L, 1), "identical")
  expect_failure(.expect_identical_fast(NA_real_, NaN), "identical")
  expect_failure(.expect_identical_fast(c(a = 1), 1), "identical")
  expect_failure(.expect_identical_fast(0, 2^-1074), "identical")
})

test_that("identity failures do not inspect retained closure environments", {
  if (!exists(".expect_identical_fast", mode = "function")) {
    fail("fast identity assertion is required")
    return()
  }
  forced <- FALSE
  env <- new.env(parent = emptyenv())
  delayedAssign("data", { forced <<- TRUE; stop("unrelated binding forced") },
    assign.env = env)
  callback <- function() NULL
  environment(callback) <- env
  actual <- list(values = 1, callback = callback)
  expected <- actual
  environment(expected$callback) <- new.env(parent = emptyenv())
  expect_failure(.expect_identical_fast(actual, expected), "identical")
  expect_false(forced)
})

test_that("identity failure messages stay bounded for long extreme vectors", {
  if (!exists(".expect_identical_fast", mode = "function")) {
    fail("fast identity assertion is required")
    return()
  }
  actual <- rep(2^1000, 2000L)
  expected <- actual
  expected[1000L] <- 0
  failure <- tryCatch(.expect_identical_fast(actual, expected,
    info = "native distance table"), expectation_failure = identity)
  expect_s3_class(failure, "expectation_failure")
  expect_match(conditionMessage(failure), "native distance table")
  expect_lt(nchar(conditionMessage(failure)), 1000L)
})
