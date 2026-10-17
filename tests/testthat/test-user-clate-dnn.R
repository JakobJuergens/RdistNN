test_that("CLATE wrapper equals explicit augmented signals and ratio engine", {
  a <- .clate_fit_arguments()
  sig <- do.call(iv_augmented_signals, a[names(.clate_test_predictions())])
  direct <- affine_ratio_dnn(X = a$X, X0 = a$X0,
                             numerator = sig$phi_y, denominator = sig$phi_d, s = a$s)
  fit <- do.call(clate_dnn, c(a, list(training = "full_sample")))
  expect_s3_class(fit, "clate_result")
  expect_s3_class(fit, "affine_ratio_result")
  expect_equal(fit$theta, direct$theta)
  expect_equal(fit$theta_vcov, direct$theta_vcov)
  expect_identical(fit$meta$component_roles, c(numerator = "reduced_form", denominator = "first_stage"))
  expect_identical(fit$meta$training, "full_sample")
})

test_that("perfect compliance reduces the estimate and uncertainty to CATE", {
  a <- .clate_fit_arguments()
  a$D <- a$A
  a$m_d0 <- rep(0, 8L)
  a$m_d1 <- rep(1, 8L)
  fit <- do.call(clate_dnn, a)
  signal <- a$m_y1 - a$m_y0 + a$A / a$e * (a$Y - a$m_y1) -
    (1 - a$A) / (1 - a$e) * (a$Y - a$m_y0)
  cate <- cate_dnn(X = a$X, X0 = a$X0, pseudo_outcome = signal, s = a$s)
  expect_equal(fit$theta, cate$tau, tolerance = 1e-12)
  expect_equal(fit$se, cate$se, tolerance = 1e-12)
  expect_equal(fit$components[, 2L], rep(1, 2L))
  a$A <- 1 - a$A
  a$e <- 1 - a$e
  a[c("m_y0", "m_y1")] <- a[c("m_y1", "m_y0")]
  a[c("m_d0", "m_d1")] <- a[c("m_d1", "m_d0")]
  flipped <- do.call(clate_dnn, a)
  expect_equal(flipped$theta, fit$theta)
  expect_equal(flipped$se, fit$se)
})

test_that("one-sided wrapper and optional nuisance payload remain aligned", {
  a <- .clate_fit_arguments()
  a$D <- a$D * a$A
  a$m_d0 <- NULL
  a$return <- c("nuisance", "signals")
  a$e[c(1, 2)] <- c(0, 1)
  fit <- do.call(clate_dnn, c(a, list(one_sided = TRUE, clip = .1, row_ids = letters[1:8])))
  expect_true(fit$meta$one_sided)
  expect_equal(fit$nuisance$propensity$used[c(1, 2)], c(.1, .9))
  expect_equal(fit$signals[, 2L], fit$nuisance$phi_d)
  expect_identical(fit$row_ids, letters[1:8])
  expect_identical(fit$.model$iv$A, as.integer(a$A))
})
