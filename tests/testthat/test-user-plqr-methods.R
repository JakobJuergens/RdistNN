test_that("plqr_result methods report coefficient inference when direct sandwich is available", {
  .test_rng_scope({
    d <- .make_plqr_inference_data(n = 85L, m = 2L)
    fit <- plqr_dnn(
      X = d$X, D = d$D, g = d$g, mD = d$mD, Y = d$Y,
      theta_eval = d$theta_eval, tau = d$tau, density = d$density,
      X0 = d$X0, method = "exact", k = 2L, s = 9L,
      se = TRUE
    )

    expect_equal(coef(fit), fit$theta, ignore_attr = "theta_source")
    expect_equal(vcov(fit, target = 1L), fit$theta_vcov[[1L]])

    s <- summary(fit)
    expect_s3_class(s, "summary.plqr_result")
    expect_true(is.data.frame(s$coefficients))
    expect_equal(
      colnames(s$coefficients),
      c("target", "term", "estimate", "std.error", "statistic", "p.value", "conf.low", "conf.high")
    )
    expect_equal(nrow(s$coefficients), length(d$theta_eval) * nrow(d$X0))
    expect_true(all(is.finite(s$coefficients$std.error)))

    ci <- confint(fit, parm = 1L, target = 2L)
    expect_equal(nrow(ci), 1L)
    expect_named(ci, c("target", "term", "conf.low", "conf.high"))
    expect_error(vcov(fit, target = 3L), "target")
    expect_error(confint(fit, parm = "not_a_term"), "valid coefficient names", ignore.case = TRUE)
  }, seed = 818L)
})

test_that("plqr_result methods fail clearly when direct sandwich is unavailable", {
  .test_rng_scope({
    d <- .make_plqr_inference_data(n = 70L, m = 1L)
    fit <- plqr_dnn(
      X = d$X, D = d$D, g = d$g, mD = d$mD, Y = d$Y,
      theta_eval = d$theta_eval, tau = d$tau,
      X0 = d$X0, method = "exact", k = 2L, s = 8L,
      se = FALSE
    )

    expect_null(fit$theta_se)
    expect_error(vcov(fit), "not available", ignore.case = TRUE)
    expect_error(confint(fit), "not available", ignore.case = TRUE)

    s <- summary(fit)
    expect_false(isTRUE(s$has_inference))
    expect_true(all(is.na(s$coefficients$std.error)))
  }, seed = 818L)
})

test_that("PLQR S3 methods expose supplied-theta semantics", {
  .test_rng_scope({
    d <- .make_plqr_inference_data(n = 60L, m = 2L)
    fit <- plqr_dnn(
      X = d$X, D = d$D, g = d$g, mD = d$mD, Y = d$Y,
      theta_eval = d$theta_eval, tau = d$tau, density = d$density,
      X0 = d$X0, method = "exact", k = 2L, s = 8L,
      se = TRUE
    )

    s <- summary(fit)
    expect_identical(s$theta_source, "supplied")
    expect_true(all(is.na(s$coefficients$statistic)))
    expect_true(all(is.na(s$coefficients$p.value)))
    expect_true(all(is.finite(s$coefficients$std.error)))
    txt <- capture.output(print(s))
    expect_true(any(grepl(
      "theta_source: supplied; no local PLQR root solving was performed",
      txt,
      fixed = TRUE
    )))

    expect_identical(attr(coef(fit), "theta_source"), "supplied")
    expect_identical(attr(confint(fit), "theta_source"), "supplied")

    pred <- predict(fit, X0 = d$X0)
    expect_identical(attr(pred, "theta_source"), "supplied")
    pred_se <- predict(fit, X0 = d$X0, se.fit = TRUE)
    expect_identical(attr(pred_se$fit, "theta_source"), "supplied")
  }, seed = 818L)
})
