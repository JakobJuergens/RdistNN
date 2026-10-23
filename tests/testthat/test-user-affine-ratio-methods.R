test_that("ratio coefficient covariance and intervals have stable target semantics", {
  f <- .ratio_fixture()
  expect_identical(names(coef(f)), "target_1")
  expect_equal(unname(coef(f)), 4)
  expect_equal(vcov(f), f$theta_vcov[[1L]])
  ci <- confint(f)
  expect_identical(dim(ci), c(1L, 2L))
  expect_identical(colnames(ci), c("conf.low", "conf.high"))
  expect_equal(as.numeric(ci), 4 + c(-1, 1) * stats::qnorm(.975) * sqrt(34))
  batch <- .ratio_fixture(X0 = matrix(c(0, 5), ncol = 1L), x0 = NULL)
  expect_length(vcov(batch), 2L)
  expect_equal(vcov(batch, target = 2L), batch$theta_vcov[[2L]])
  expect_equal(confint(batch, parm = 2L), confint(batch)[2L, , drop = FALSE])
  expect_error(vcov(batch, target = 3L), "target")
  expect_error(confint(batch, parm = 0L), "parm")
  for (level in list(0, 1, NA_real_, "x", c(.9, .95))) {
    expect_error(confint(f, level = level), "level")
    expect_error(summary(f, level = level), "level")
    expect_error(predict(f, x0 = 0, level = level), "level")
  }
})

test_that("confidence quantiles remain finite at the largest representable level below one", {
  f <- .ratio_fixture()
  level <- 1 - .Machine$double.eps / 2
  ci <- confint(f, level = level)
  critical <- stats::qnorm((1 - level) / 2, lower.tail = FALSE)
  expect_true(all(is.finite(ci)))
  expect_equal(as.numeric(ci), 4 + c(-1, 1) * critical * sqrt(34))
})

test_that("ratio summaries retain undefined and degenerate targets without fabricated intervals", {
  .with_distance_tie_warning({
  expect_warning(f <- .ratio_fixture(X0 = matrix(c(0, 3, 5), ncol = 1L), x0 = NULL,
                       denominator = c(0, 0, -1, 2), s = 3L), class = "rdistnn_partial_result_warning")
  expect_warning(tab <- summary(f), class = "rdistnn_partial_result_warning")
  expect_s3_class(tab, "summary.affine_ratio_result")
  expect_type(tab, "list")
  expect_equal(tab$status$target_id, 1:3)
  expect_true(is.na(tab$coefficients$estimate[1L]))
  expect_warning(ci <- confint(f), class = "rdistnn_partial_result_warning")
  expect_true(all(is.na(ci[1L, ])))
  expect_equal(tab$components$denominator_se, vapply(f$component_vcov, function(v) sqrt(v[2L, 2L]), numeric(1L)))
  expect_warning(expect_output(ret <- withVisible(print(f)), "zero_denominator"),
    class = "rdistnn_partial_result_warning")
  expect_false(ret$visible)
  .expect_identical_fast(ret$value, f)
  expect_output(ret <- withVisible(print(tab)), "zero_denominator")
  expect_false(ret$visible)
  expect_warning(degenerate <- .ratio_fixture(numerator = 4 * c(1, -1, 2, 1)),
    class = "rdistnn_partial_result_warning")
  expect_equal(degenerate$se, 0)
  expect_warning(ci <- confint(degenerate), class = "rdistnn_partial_result_warning")
  expect_true(all(is.na(ci)))
  no_cov <- .ratio_fixture(se = FALSE)
  point_table <- summary(no_cov)
  expect_equal(point_table$coefficients$estimate, no_cov$theta)
  expect_true(all(is.na(point_table$coefficients$std.error)))
  expect_error(vcov(no_cov), "covariance")
  expect_error(confint(no_cov), "covariance")
  })
})

test_that("prediction relocalizes frozen arrays and survives serialization", {
  f <- .ratio_fixture()
  expect_equal(predict(f, x0 = 0), f$theta)
  expect_equal(predict(f, x0 = 0, se.fit = TRUE), list(fit = f$theta, se.fit = f$se))
  at <- matrix(c(0, 5), ncol = 1L)
  direct <- .ratio_fixture(X0 = at, x0 = NULL)
  expect_equal(predict(f, X0 = at), direct$theta)
  p <- predict(f, X0 = at, se.fit = TRUE, interval = "confidence")
  expect_equal(p$fit[, "fit"], direct$theta)
  expect_equal(p$fit[, "lwr"], as.numeric(confint(direct)[, 1L]))
  expect_equal(p$se.fit, direct$se)
  path <- tempfile(fileext = ".rds")
  on.exit(unlink(path), add = TRUE)
  saveRDS(f, path)
  expect_equal(predict(readRDS(path), X0 = at, se.fit = TRUE), predict(f, X0 = at, se.fit = TRUE))
  expect_error(predict(f, newdata = at, X0 = at), "newdata")
  expect_error(predict(f), "x0|X0|newdata")
  expect_error(predict(f, x0 = 0, se.fit = 1), "se.fit")
  bad <- f
  bad$.model <- NULL
  expect_error(predict(bad, x0 = 0), "model")
})

test_that("data-frame CLATE prediction uses reporting columns only and never refits", {
  a <- .clate_fit_arguments()
  fit <- do.call(clate_dnn, .clate_frame_arguments(a, .clate_test_frame(a)))
  targets <- data.frame(unused = NA_real_, x2 = a$X0[, 2L], x1 = a$X0[, 1L])
  expect_equal(predict(fit, newdata = targets, se.fit = TRUE),
                list(fit = fit$theta, se.fit = fit$se))
  expect_error(predict(fit, newdata = targets["x1"]), "x2")
  point <- .ratio_fixture(se = FALSE)
  expect_equal(predict(point, x0 = 0, se.fit = TRUE)$se.fit, .ratio_fixture()$se)
  nearest <- .ratio_fixture(s = 4L, se = FALSE)
  expect_error(predict(nearest, x0 = 0, se.fit = TRUE), "s.*n|s < n")
})

test_that("numeric reporting selectors retain column identity for newdata", {
  d <- data.frame(junk = 10:13, x = 1:4, y = c(2, 1, 4, 3), d = c(1, -1, 2, 1))
  fit <- affine_ratio_dnn(data = d, x_cols = 2L, numerator_col = 3L,
                          denominator_col = 4L, x0 = 0, s = 2L)
  expect_equal(predict(fit, newdata = data.frame(x = 0, junk = 5)), 4)
  expect_equal(predict(fit, newdata = data.frame(x = 0)), 4)
  path <- tempfile(fileext = ".rds")
  on.exit(unlink(path), add = TRUE)
  saveRDS(fit, path)
  expect_equal(predict(readRDS(path), newdata = data.frame(x = 0), se.fit = TRUE),
               list(fit = 4, se.fit = sqrt(34)))
})

test_that("ratio reuse requires retained normalization controls", {
  fit <- .ratio_fixture(method = "asymp", normalize = FALSE)
  expect_equal(predict(fit, x0 = 0), fit$theta)
  expect_equal(diagnostics(fit)$n_eff, fit$neff)
  for (normalize in list(NULL, NA, 1, logical(), c(TRUE, FALSE))) {
    bad <- fit
    bad$.model$normalize <- normalize
    expect_error(predict(bad, x0 = 0), "normalize.*[Rr]efit")
    expect_error(diagnostics(bad), "normalize.*[Rr]efit")
  }
})

test_that("ratio reuse requires retained jackknife controls and supports point-only uncertainty", {
  fit <- .ratio_fixture()
  for (jackknife in list(NULL, NA_character_, "retired", character(),
                        c("fast", "loo"), 1)) {
    bad <- fit
    bad$meta$jackknife <- jackknife
    expect_error(predict(bad, x0 = 0), "jackknife.*[Rr]efit")
    expect_error(diagnostics(bad), "jackknife.*[Rr]efit")
  }
  point <- .ratio_fixture(se = FALSE)
  expect_identical(point$meta$jackknife, "none")
  expect_equal(predict(point, x0 = 0, se.fit = TRUE)$se.fit, .ratio_fixture()$se)
})
