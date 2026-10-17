test_that("asymptotic DNN agrees with an analytical four-observation fixture", {
  fit <- .ratio_fixture(method = "asymp", return = c("W", "deletions"))
  expect_equal(as.numeric(fit$W), c(7, 5, 3, 1) / 16, tolerance = 1e-13)
  expect_equal(as.numeric(fit$components), c(17 / 8, 9 / 16), tolerance = 1e-13)
  expect_equal(fit$theta, 34 / 9, tolerance = 1e-12)
  want <- cbind(c(20, 25, 16, 17) / 9, c(2, 12, 3, 4) / 9)
  expect_equal(unname(fit$deletions[, , 1L]), want, tolerance = 1e-12)
  centered <- sweep(want, 2L, colMeans(want), "-")
  covariance <- 3 / 4 * crossprod(centered)
  expect_equal(unname(fit$component_vcov[[1L]]), covariance, tolerance = 1e-12)
  gradient <- c(1, -34 / 9) / (9 / 16)
  expect_equal(fit$se^2, as.numeric(crossprod(gradient, covariance %*% gradient)), tolerance = 1e-12)
  expect_equal(fit$se^2, 3075200 / 177147, tolerance = 1e-12)
})

test_that("effect estimators expose single-scale DNN only", {
  expect_error(.ratio_fixture(method = "tdnn"), "method")
  expect_error(do.call(clate_dnn, c(.clate_fit_arguments(), list(method = "tdnn"))), "method")
  expect_error(cate_dnn(X = matrix(1:8), x0 = 0, pseudo_outcome = 1:8,
    s = 3L, method = "tdnn"), "exact.*asymp")
})

test_that("general k uses the same exact and approximate weights as existing estimators", {
  expected <- list(exact = c(3, 3, 2, 0) / 8, asymp = c(47, 41, 29, 11) / 128)
  for (method in names(expected)) for (normalize in c(TRUE, FALSE)) {
    fit <- .ratio_fixture(method = method, k = 2L, s = 3L,
                          normalize = normalize, return = c("W", "neighbors"))
    expect_equal(as.numeric(fit$W), expected[[method]], tolerance = 1e-13)
    first <- np_dnn(X = matrix(1:4), x0 = 0, Y = c(2, 1, 4, 3), s = 3L,
                    k = 2L, method = method, normalize = normalize)
    second <- np_dnn(X = matrix(1:4), x0 = 0, Y = c(1, -1, 2, 1), s = 3L,
                     k = 2L, method = method, normalize = normalize)
    expect_equal(fit$components[1L, ], c(numerator = first$mu, denominator = second$mu))
    expect_equal(diag(fit$component_vcov[[1L]]), c(first$se^2, second$se^2), ignore_attr = TRUE)
    expect_identical(dim(fit$neighbors), c(2L, 1L))
    expect_identical(fit$meta$method, method)
    expect_identical(fit$meta$normalize, normalize)
    expect_equal(fit$meta$k, 2L)
  }
})

test_that("prediction and diagnostics preserve general-k asymptotic weighting", {
  fit <- .ratio_fixture(method = "asymp", k = 2L, s = 3L, normalize = FALSE,
                        return = c("W", "pseudo_values"))
  direct <- .ratio_fixture(method = "asymp", k = 2L, s = 3L, normalize = FALSE,
                           X0 = matrix(c(0, 5), ncol = 1L), x0 = NULL)
  expect_equal(predict(fit, X0 = direct$x0, se.fit = TRUE), list(fit = direct$theta, se.fit = direct$se))
  diag <- diagnostics(fit)
  expect_equal(diag$neff, 1 / sum(c(47, 41, 29, 11)^2 / 128^2))
  expect_equal(diag$component_covariance, fit$component_vcov[[1L]][1L, 2L])
  q <- fit$theta * fit$pseudo_values[, 2L, 1L] - fit$pseudo_values[, 1L, 1L]
  shares <- (q - mean(q))^2 / sum((q - mean(q))^2)
  expect_equal(diag$jk_top1_sq_share, max(shares), tolerance = 1e-12)
  expect_identical(diag$method, "asymp")
  expect_equal(diag$k, 2L)
})

test_that("general-k exact DNN matches literal subset averages and deleted covariance", {
  for (n in 3:7) for (s in 2:(n - 1L)) for (k in 2:s) {
    X <- matrix(seq_len(n) %% 3L)
    signal <- cbind(sin(seq_len(n)), 1 + seq_len(n) / n)
    fit <- suppressWarnings(affine_ratio_dnn(X = X, x0 = 0,
      numerator = signal[, 1L], denominator = signal[, 2L], k = k, s = s,
      return = "deletions"))
    full <- .clate_subset_average(X, 0, signal, s, k)
    loo <- t(vapply(seq_len(n), function(i) {
      .clate_subset_average(X[-i, , drop = FALSE], 0, signal[-i, , drop = FALSE], s, k)
    }, numeric(2L)))
    centered <- sweep(loo, 2L, colMeans(loo), "-")
    expect_equal(as.numeric(fit$components), full, tolerance = 1e-12)
    expect_equal(unname(fit$deletions[, , 1L]), loo, tolerance = 1e-12)
    expect_equal(unname(fit$component_vcov[[1L]]), (n - 1) / n * crossprod(centered), tolerance = 1e-12)
  }
})

test_that("both DNN schemes reduce CLATE estimates and inference to existing CATE", {
  a <- .clate_fit_arguments()
  a$D <- a$A
  a$m_d0 <- rep(0, 8L)
  a$m_d1 <- rep(1, 8L)
  a$return <- "signals"
  for (method in c("exact", "asymp")) for (k in 1:2) for (normalize in c(TRUE, FALSE)) {
    fit <- do.call(clate_dnn, c(a, list(method = method, k = k, normalize = normalize)))
    cate <- cate_dnn(X = a$X, X0 = a$X0, pseudo_outcome = fit$.model$signal[, 1L],
                      s = a$s, k = k, method = method, normalize = normalize)
    expect_equal(fit$theta, cate$tau, tolerance = 1e-12)
    expect_equal(fit$se, cate$se, tolerance = 1e-12)
    expect_equal(summary(fit)$coefficients, summary(cate)$coefficients, tolerance = 1e-12)
    expect_equal(confint(fit), confint(cate), tolerance = 1e-12)
    expect_equal(predict(fit, X0 = a$X0, se.fit = TRUE), predict(cate, X0 = a$X0, se.fit = TRUE))
    frame <- .clate_frame_arguments(a, .clate_test_frame(a))
    tabular <- do.call(clate_dnn, c(frame, list(method = method, k = k, normalize = normalize)))
    expect_equal(tabular$components, fit$components)
    expect_equal(tabular$theta_vcov, fit$theta_vcov)
  }
})

test_that("asymptotic and general-k paired engines agree across R compiled and frozen loops", {
  old <- options(RdistNN.use_compiled = FALSE)
  on.exit(options(old), add = TRUE)
  for (method in c("exact", "asymp")) for (k in 1:2) for (n in c(4L, 13L)) {
    args <- list(X = cbind(seq_len(n) %% 3L, seq_len(n) / n),
      X0 = rbind(c(0, 0), c(2, .5)), numerator = sin(seq_len(n)),
      denominator = 1 + seq_len(n) / n, method = method, k = k, s = n - 1L,
      return = c("deletions", "pseudo_values"))
    options(RdistNN.use_compiled = FALSE)
    fast <- do.call(affine_ratio_dnn, args)
    loop <- do.call(affine_ratio_dnn, c(args, list(jackknife = "loo")))
    options(RdistNN.use_compiled = TRUE)
    compiled <- do.call(affine_ratio_dnn, args)
    for (field in c("theta", "se", "deletions", "pseudo_values", "component_vcov", "theta_vcov", "status")) {
      expect_equal(fast[[field]], loop[[field]], tolerance = 1e-12)
      expect_equal(fast[[field]], compiled[[field]], tolerance = 1e-11)
    }
  }
})
