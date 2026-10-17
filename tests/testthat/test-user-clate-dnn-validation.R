test_that("CLATE validates every supplied vector at the public boundary", {
  original <- .clate_fit_arguments()
  for (nm in names(.clate_test_predictions())) {
    for (bad in list(NA_real_, NaN, Inf, -Inf)) {
      a <- original
      a[[nm]][1L] <- bad
      expect_error(do.call(clate_dnn, a), nm)
    }
    a <- original
    a[[nm]] <- matrix(a[[nm]], ncol = 1L)
    expect_error(do.call(clate_dnn, a), nm)
  }
  for (nm in c("D", "A")) {
    a <- original
    a[[nm]] <- factor(a[[nm]])
    expect_error(do.call(clate_dnn, a), nm)
  }
  for (nm in c("e", "m_d0", "m_d1")) {
    a <- original
    a[[nm]][1L] <- 1.1
    expect_error(do.call(clate_dnn, a), nm)
  }
})

test_that("one-sided CLATE checks both observed coding and supplied zero-arm means", {
  a <- .clate_fit_arguments()
  expect_error(do.call(clate_dnn, c(a, list(one_sided = TRUE))), "one_sided")
  a$D <- a$D * a$A
  expect_error(do.call(clate_dnn, c(a, list(one_sided = TRUE))), "m_d0")
  a$m_d0 <- NULL
  expect_error(do.call(clate_dnn, a), "m_d0")
  expect_s3_class(do.call(clate_dnn, c(a, list(one_sided = TRUE))), "clate_result")
})

test_that("CLATE unknown returns warn once while optional payloads remain additive", {
  a <- .clate_fit_arguments()
  baseline <- do.call(clate_dnn, a)
  a$return <- c("nuisance", "signals", "unknown")
  warnings <- character()
  z <- withCallingHandlers(do.call(clate_dnn, a), warning = function(w) {
    warnings <<- c(warnings, conditionMessage(w))
    invokeRestart("muffleWarning")
  })
  expect_length(warnings, 1L)
  expect_match(warnings, "unknown")
  expect_equal(z$theta, baseline$theta)
  expect_equal(z$theta_vcov, baseline$theta_vcov)
  expect_s3_class(z$nuisance, "iv_signals")
})

test_that("CLATE default tuning is resolved once and recorded", {
  a <- .clate_fit_arguments()
  a$s <- NULL
  expect_warning(z <- do.call(clate_dnn, a), "uniform")
  expect_equal(z$meta$s, 1L)
  expect_equal(z$meta$tuning$source, "paper_default")
  a$s <- 2L
  a$tuning <- tuning_dnn_paper()
  expect_error(do.call(clate_dnn, a), "s/tuning")
})
