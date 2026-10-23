.retained_report_fits <- function(se = TRUE, method = "exact") {
  X <- matrix(seq_len(8L), ncol = 1L)
  X0 <- matrix(c(0, 9, 0), ncol = 1L)
  y <- seq_len(8L)^2
  tuning <- new_tuning("retained_report_state", eval = function(...) list(s = 3L),
    eval_grid = function(...) data.frame(s = c(2L, 3L)))
  common <- list(X = X, X0 = X0, se = se, method = method, tuning = tuning,
    tuning_mode = "grid", return = character())
  constructors <- list(
    np = function(...) np_dnn(Y = y, ...),
    cate = function(...) cate_dnn(pseudo_outcome = y, ...),
    ratio = function(...) affine_ratio_dnn(numerator = y, denominator = rep(1, 8L), ...),
    clate = function(...) clate_dnn(Y = y, D = rep(c(0, 1), 4L), A = rep(c(0, 1), 4L),
      e = rep(.5, 8L), m_y0 = rep(0, 8L), m_y1 = rep(1, 8L),
      m_d0 = rep(0, 8L), m_d1 = rep(1, 8L), ...))
  fits <- lapply(constructors, function(fn) .with_localization_policy_warnings(do.call(fn, common)))
  fits$plqr <- plqr_dnn(X = X, X0 = X0, Y = rep(1, 8L),
    D = matrix(seq_len(8L)), mD = matrix(0, 8L, 1L), g = rep(0, 8L),
    theta_eval = 0, density = rep(1, 8L), s = 3L, se = se, method = method, return = character())
  fits$tdnn <- np_tdnn(X = X, X0 = X0, Y = y, s1 = 2L, s2 = 3L,
    se = se, method = method, return = character())
  fits
}

.retained_report_routes <- function(fit) {
  routes <- list(coef = coef, summary = summary, confint = confint, vcov = vcov,
    estimates = plot_data_estimates, diagnostics = diagnostics,
    prediction = function(x) predict(x, x0 = .25))
  if (!is.null(fit$tuning_grid)) routes$candidate <- function(x) candidate_fit(x, 1L)
  routes
}

test_that("current canonical report arrays reject malformed state on public consumers", {
  fits <- .retained_report_fits()
  for (family in names(fits)) {
    estimate <- if (family == "plqr") "theta" else if (family %in% c("ratio", "clate")) "theta" else if (family == "cate") "tau" else "mu"
    uncertainty <- if (family == "plqr") "theta_se" else "se"
    for (case in c("missing_estimate", "short_estimate", "nonfinite_estimate",
                   "missing_se", "short_se", "negative_se", "infinite_se",
                   "missing_targets", "short_targets", "nonfinite_targets", "target_count")) {
      bad <- fits[[family]]
      if (case == "missing_estimate") bad[estimate] <- NULL
      if (case == "short_estimate") bad[[estimate]] <- if (family == "plqr") matrix(0, 1L, 1L) else 0
      if (case == "nonfinite_estimate") bad[[estimate]][1L] <- Inf
      if (case == "missing_se") bad[uncertainty] <- NULL
      if (case == "short_se") bad[[uncertainty]] <- if (family == "plqr") matrix(1, 1L, 1L) else 1
      if (case == "negative_se") bad[[uncertainty]][1L] <- -2
      if (case == "infinite_se") bad[[uncertainty]][1L] <- Inf
      if (case == "missing_targets") bad["x0"] <- NULL
      if (case == "short_targets") bad$x0 <- bad$x0[1L, , drop = FALSE]
      if (case == "nonfinite_targets") bad$x0[1L] <- NA_real_
      if (case == "target_count") bad$meta$m <- 1L
      before <- serialize(bad, NULL)
      for (route in .retained_report_routes(bad)) {
        error <- tryCatch({ .with_localization_policy_warnings(route(bad)); NULL }, error = identity)
        expect_true(inherits(error, "error") && grepl("[Rr]efit", conditionMessage(error)),
          info = paste(family, case))
      }
      .expect_identical_fast(serialize(bad, NULL), before, info = paste(family, case))
    }
  }
})

test_that("valid current arrays survive serialization, repeated targets and disabled uncertainty", {
  .with_partial_result_warning({
  for (se in c(TRUE, FALSE)) for (fit in .retained_report_fits(se)) {
    restored <- unserialize(serialize(fit, NULL))
    expect_equal(coef(restored), coef(fit))
    expect_equal(summary(restored), summary(fit))
    expect_equal(.with_localization_policy_warnings(plot_data_estimates(restored)), .with_localization_policy_warnings(plot_data_estimates(fit)))
    expect_equal(.with_localization_policy_warnings(predict(restored, x0 = .25)), .with_localization_policy_warnings(predict(fit, x0 = .25)))
    if (!is.null(fit$tuning_grid)) {
      candidate <- .with_localization_policy_warnings(candidate_fit(restored, 1L))
      reference <- .with_localization_policy_warnings(candidate_fit(fit, 1L))
      expect_equal(coef(candidate), coef(reference))
      expect_equal(summary(candidate), summary(reference))
      expect_identical(candidate$meta$candidate_id, reference$meta$candidate_id)
      numerical <- names(reference$.model)[vapply(reference$.model, is.numeric, logical(1L))]
      .expect_identical_fast(candidate$.model[numerical], reference$.model[numerical])
      .expect_identical_fast(candidate$.model$geometry, reference$.model$geometry)
      expect_s3_class(summary(candidate), class(summary(fit)))
      expect_equal(.with_localization_policy_warnings(predict(candidate, x0 = .25)),
        .with_localization_policy_warnings(predict(candidate_fit(fit, 1L), x0 = .25)))
    }
  }
  fit <- affine_ratio_dnn(X = matrix(1:4), X0 = matrix(c(0, 5)),
    numerator = rep(1, 4), denominator = c(0, 0, 1, 1), s = 4, se = FALSE)
  expect_true(is.na(coef(fit)[1L]))
  expect_s3_class(summary(fit), "summary.affine_ratio_result")
  expect_equal(.with_localization_policy_warnings(predict(fit, X0 = fit$x0)), unname(coef(fit)))
  failed <- affine_ratio_dnn(X = matrix(1:4), x0 = 0,
    numerator = rep(1, 4), denominator = rep(0, 4), s = 2)
  expect_true(all(is.na(summary(failed)$coefficients$std.error)))
  })
})

test_that("paired and PLQR covariance arrays retain target and component dimensions", {
  fits <- .retained_report_fits()[c("ratio", "clate", "plqr")]
  for (family in names(fits)) {
    fields <- if (family == "plqr") "theta_vcov" else c("theta_vcov", "component_vcov")
    for (field in fields) for (case in c("short_list", "wrong_matrix", "negative_diagonal")) {
      bad <- fits[[family]]
      if (case == "short_list") bad[[field]] <- bad[[field]][1L]
      if (case == "wrong_matrix") bad[[field]][[1L]] <- matrix(0, 3L, 3L)
      if (case == "negative_diagonal") bad[[field]][[1L]][1L, 1L] <- -1
      for (route in .retained_report_routes(bad)) {
        error <- tryCatch({ .with_localization_policy_warnings(route(bad)); NULL }, error = identity)
        expect_true(inherits(error, "error") && grepl("[Rr]efit", conditionMessage(error)),
          info = paste(family, field, case))
      }
    }
  }
})


test_that("moment-only PLQR candidate views keep explicit disabled uncertainty", {
  tuning <- new_tuning("retained_moment_candidate", eval = function(...) list(s = 3L),
    eval_grid = function(...) data.frame(s = c(2L, 3L)))
  fit <- plqr_dnn(X = matrix(1:8), x0 = 0, Y = rep(1, 8), D = matrix(1:8),
    mD = matrix(0, 8, 1), g = rep(0, 8), theta_eval = 0, se = FALSE,
    tuning = tuning, tuning_mode = "grid", return = character())
  candidate <- candidate_fit(fit, 1L)
  expect_true("moment_se" %in% names(candidate))
  expect_null(candidate$moment_se)
  expect_s3_class(summary(candidate), "summary.plqr_result")
})

test_that("scalar ratio candidate plot tables reject duplicate target rows", {
  fits <- .retained_report_fits()[c("ratio", "clate")]
  for (fit in fits) {
    bad <- fit
    bad$tuning_grid$estimates <- rbind(bad$tuning_grid$estimates,
      bad$tuning_grid$estimates[1L, , drop = FALSE])
    expect_error(.with_localization_policy_warnings(plot_data_estimates(bad, candidates = "all")), "duplicate target_id")
    expect_error(.with_localization_policy_warnings(plot_data_estimates(bad, candidates = bad$tuning_grid$candidates$candidate_id[[1L]])), "duplicate target_id")
    expect_equal(nrow(.with_localization_policy_warnings(plot_data_estimates(fit, candidates = "all"))),
      nrow(fit$tuning_grid$candidates) * nrow(fit$x0))
  }
})

test_that("paired uncertainty arrays agree with the retained covariance mode", {
  .with_partial_result_warning({
  fields <- c("se", "component_vcov", "theta_vcov", "theta_omega")
  for (fit in .retained_report_fits()[c("ratio", "clate")]) {
    expect_identical(fit$meta$covariance, "paired_frozen_score")
    expect_true(any(fit$status$wald_available))
    corruptions <- lapply(fields, function(field) {
      bad <- fit
      bad[field] <- list(NULL)
      bad
    })
    jointly_null <- fit
    jointly_null[fields] <- rep(list(NULL), length(fields))
    jointly_null$status$wald_available <- FALSE
    corruptions <- c(corruptions, list(jointly_null))
    for (bad in corruptions) for (route in .retained_report_routes(bad)) {
      expect_error(.with_localization_policy_warnings(route(bad)), "[Rr]efit")
    }
    for (mode in list(NULL, "", "unknown", c("none", "paired_frozen_score"))) {
      bad <- fit
      bad$meta["covariance"] <- list(mode)
      expect_error(summary(bad), "meta\\$covariance.*[Rr]efit")
    }
    disabled <- fit
    disabled$meta$covariance <- "none"
    expect_error(summary(disabled), "[Rr]efit")
  }
  for (fit in .retained_report_fits(FALSE)[c("ratio", "clate")]) {
    expect_identical(fit$meta$covariance, "none")
    expect_true(all(vapply(fit[fields], is.null, logical(1L))))
    expect_s3_class(summary(fit), "summary.affine_ratio_result")
    bad <- fit
    bad$status$wald_available[[1L]] <- TRUE
    expect_error(summary(bad), "[Rr]efit")
  }
  unavailable <- affine_ratio_dnn(X = matrix(1:4), x0 = 0,
    numerator = c(1, 2, 3, 4) * 1e-200, denominator = rep(1, 4), s = 2)
  expect_identical(unavailable$status$covariance_status, "numerical_failure")
  expect_true(is.na(unavailable$se))
  expect_identical(unavailable$meta$covariance, "paired_frozen_score")
  expect_s3_class(summary(unavailable), "summary.affine_ratio_result")
  expect_true(all(is.na(confint(unavailable))))
  })
})

test_that("ordinary scalar uncertainty agrees with retained jackknife mode", {
  for (fit in .retained_report_fits()[c("np", "cate", "tdnn")]) {
    expect_identical(fit$meta$jackknife, "fast")
    bad <- fit
    bad["se"] <- list(NULL)
    for (route in .retained_report_routes(bad)) {
      expect_error(.with_localization_policy_warnings(route(bad)), "[Rr]efit")
    }
    disabled <- fit
    disabled$meta$jackknife <- "none"
    expect_error(summary(disabled), "[Rr]efit")
    for (mode in list(NULL, "", "unknown", c("fast", "none"))) {
      bad <- fit
      bad$meta["jackknife"] <- list(mode)
      expect_error(summary(bad), "meta\\$jackknife.*[Rr]efit")
    }
  }
  for (fit in .retained_report_fits(FALSE)[c("np", "cate", "tdnn")]) {
    expect_identical(fit$meta$jackknife, "none")
    expect_null(fit$se)
    expect_s3_class(summary(fit), if (inherits(fit, "cate_result")) "summary.cate_result" else "summary.np_result")
  }
})

test_that("paired reports reject missing and inconsistent retained jackknife modes", {
  for (se in c(TRUE, FALSE)) for (fit in .retained_report_fits(se)[c("ratio", "clate")]) {
    modes <- c(list(NULL, "", "unknown", c("fast", "none")),
      list(if (se) "none" else "fast"))
    for (mode in modes) {
      bad <- fit
      bad$meta["jackknife"] <- list(mode)
      for (route in .retained_report_routes(bad)) {
        expect_error(.with_localization_policy_warnings(route(bad)),
          class = "rdistnn_state_error")
      }
    }
  }
})

test_that("paired retained weight methods reject missing, invalid and inconsistent controls", {
  fits <- .retained_report_fits(method = "asymp")[c("ratio", "clate")]
  exact <- .retained_report_fits(method = "exact")[c("ratio", "clate")]
  for (family in names(fits)) {
    fit <- fits[[family]]
    # This fixture distinguishes the former NULL-to-exact prediction fallback.
    expect_gt(max(abs(predict(fit, X0 = fit$x0) - predict(exact[[family]], X0 = fit$x0))), .01)
    cases <- c("missing", "null", "na", "vector", "invalid", "mismatch")
    for (location in c(".model", "meta")) for (case in cases) {
      bad <- fit
      if (case == "missing") bad[[location]]["method"] <- NULL
      if (case == "null") bad[[location]]["method"] <- list(NULL)
      if (case == "na") bad[[location]]$method <- NA_character_
      if (case == "vector") bad[[location]]$method <- c("exact", "asymp")
      if (case == "invalid") bad[[location]]$method <- "unrecognized"
      if (case == "mismatch") bad[[location]]$method <- "exact"
      before <- serialize(bad, NULL)
      routes <- c(.retained_report_routes(bad), list(
        prediction_se = function(x) predict(x, X0 = x$x0, se.fit = TRUE),
        prediction_interval = function(x) predict(x, X0 = x$x0, interval = "confidence")))
      for (route in routes) {
        error <- tryCatch({ .with_localization_policy_warnings(route(bad)); NULL }, error = identity)
        expect_true(inherits(error, "error") && grepl("[Rr]efit", conditionMessage(error)),
          info = paste(family, location, case))
      }
      .expect_identical_fast(serialize(bad, NULL), before, info = paste(family, location, case))
    }
  }
})

test_that("paired exact and asymptotic controls survive serialized point and SE prediction", {
  for (method in c("exact", "asymp")) for (se in c(TRUE, FALSE)) {
    fits <- .retained_report_fits(se = se, method = method)[c("ratio", "clate")]
    for (fit in fits) {
      expect_identical(fit$.model$method, method)
      expect_identical(fit$meta$method, method)
      restored <- unserialize(serialize(fit, NULL))
      point <- predict(fit, X0 = fit$x0)
      uncertainty <- predict(fit, X0 = fit$x0, se.fit = TRUE)
      expect_equal(point, unname(coef(fit)))
      expect_identical(predict(restored, X0 = fit$x0), point)
      expect_identical(predict(restored, X0 = fit$x0, se.fit = TRUE), uncertainty)
      if (se) expect_equal(uncertainty$se.fit, fit$se)
      candidate <- candidate_fit(restored, 1L)
      expect_identical(candidate$.model$method, method)
      expect_identical(candidate$meta$method, method)
      expect_equal(predict(candidate, X0 = fit$x0, se.fit = TRUE),
        predict(candidate_fit(fit, 1L), X0 = fit$x0, se.fit = TRUE))
    }
  }
})
