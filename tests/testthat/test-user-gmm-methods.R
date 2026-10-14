# GMM result selectors, inference availability, and summaries.
test_that("gmm_result methods validate availability and selectors", {
  dgp <- sim_iv_data("dgp2_iv_homog_overid", n = 100L, kx = 2L, q = 2L, seed = 27L)
  model <- .test_dr_iv_model(dgp$q)
  fit <- gmm_dnn(
    data = dgp$data,
    x_cols = c("X1", "X2"),
    X0 = rbind(dgp$x0, dgp$x0 + c(0.1, 0.1)),
    model = model,
    k = 2L,
    s = 9L
  )

  expect_equal(coef(fit), fit$theta)
  expect_equal(length(vcov(fit)), 2L)
  expect_equal(vcov(fit, target = 1L), fit$vcov[[1L]])
  expect_error(vcov(fit, target = 3L), "target")

  ci <- confint(fit, parm = "theta", target = 1L)
  expect_equal(nrow(ci), 1L)
  expect_named(ci, c("target", "term", "conf.low", "conf.high"))
  expect_error(confint(fit, parm = 2L), "parm")
  expect_error(confint(fit, parm = 1.5), "integer coefficient indices", ignore.case = TRUE)
  expect_error(confint(fit, parm = "theta_typo"), "valid coefficient names", ignore.case = TRUE)
  expect_error(confint(fit, parm = character()), "valid coefficient names", ignore.case = TRUE)
  expect_error(confint(fit, target = 1.5), "integer target indices", ignore.case = TRUE)
  expect_error(confint(fit, target = Inf), "non-finite", ignore.case = TRUE)
  expect_error(confint(fit, target = 1e20), "valid target indices", ignore.case = TRUE)
  expect_error(confint(fit, target = c(0L, 1L)), "valid target indices", ignore.case = TRUE)
  expect_error(summary(fit, level = 1), "level")

  fit_no_se <- gmm_dnn(
    data = dgp$data,
    x_cols = c("X1", "X2"),
    x0 = dgp$x0,
    model = model,
    k = 2L,
    s = 9L,
    se = FALSE
  )
  expect_error(vcov(fit_no_se), "not available")
  expect_error(confint(fit_no_se), "not available")
})

test_that("gmm_result summaries expose inference and overidentification state", {
  overid_dgp <- sim_iv_data(
    "dgp2_iv_homog_overid", n = 100L, kx = 2L, q = 2L, seed = 34L
  )
  fit_overid <- gmm_dnn(
    data = overid_dgp$data,
    x_cols = c("X1", "X2"),
    x0 = overid_dgp$x0,
    model = .test_dr_iv_model(overid_dgp$q),
    k = 2L,
    s = 9L
  )
  summary_overid <- summary(fit_overid)

  expect_s3_class(summary_overid, "summary.gmm_result")
  expect_named(
    summary_overid$coefficients,
    c(
      "target", "term", "estimate", "std.error", "statistic", "p.value",
      "conf.low", "conf.high"
    )
  )
  expect_named(
    summary_overid$overid,
    c("target", "J", "df", "p.value", "neff")
  )
  expect_true(any(is.finite(summary_overid$overid$J)))

  fit_text <- capture.output(fit_returned <- print(fit_overid))
  summary_text <- capture.output(summary_returned <- print(summary_overid))
  expect_identical(fit_returned, fit_overid)
  expect_identical(summary_returned, summary_overid)
  expect_true(any(grepl("<gmm_result>", fit_text, fixed = TRUE)))
  expect_true(any(grepl("<summary.gmm_result>", summary_text, fixed = TRUE)))
  expect_true(any(grepl(
    "Overidentification diagnostics:", summary_text, fixed = TRUE
  )))

  justid_dgp <- sim_iv_data(
    "dgp1_iv_justid_linear", n = 90L, kx = 2L, q = 1L, seed = 35L
  )
  fit_justid <- gmm_dnn(
    data = justid_dgp$data,
    x_cols = c("X1", "X2"),
    x0 = justid_dgp$x0,
    model = .test_dr_iv_model(justid_dgp$q),
    k = 2L,
    s = 9L,
    se = FALSE
  )
  summary_justid <- summary(fit_justid)
  justid_fit_text <- capture.output(justid_returned <- print(fit_justid))
  justid_summary_text <- capture.output(print(summary_justid))
  inferential <- c(
    "std.error", "statistic", "p.value", "conf.low", "conf.high"
  )

  expect_identical(justid_returned, fit_justid)
  expect_true(any(grepl("SEs not requested", justid_fit_text, fixed = TRUE)))
  expect_true(all(is.na(summary_justid$coefficients[inferential])))
  expect_true(all(is.na(summary_justid$overid$J)))
  expect_false(any(grepl(
    "Overidentification diagnostics:", justid_summary_text, fixed = TRUE
  )))
})
