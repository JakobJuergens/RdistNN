.estimate_range_mean <- function(family, X, signal, target, method, se,
                                 k = 1L, s = k, s1 = 1L, s2 = 2L) {
  .with_localization_policy_warnings(switch(family,
    np = np_dnn(X = X, x0 = target, Y = signal, method = method, k = k, s = s, se = se),
    cate = cate_dnn(X = X, x0 = target, pseudo_outcome = signal, method = method, k = k, s = s, se = se),
    tdnn = np_tdnn(X = X, x0 = target, Y = signal, method = method, s1 = s1, s2 = s2, se = se)))
}

.estimate_range_outcome <- function(code, arg = "[Rr]escale") {
  out <- tryCatch(.with_localization_policy_warnings(code), error = identity)
  if (inherits(out, "error")) expect_match(conditionMessage(out), arg)
  out
}

test_that("actual overflowing weighted signals have mean errors and paired component status", {
  old <- options(RdistNN.use_compiled = FALSE)
  on.exit(options(old), add = TRUE)
  # The exact weighted value is 2 * xmax, genuinely outside double range.
  # A second finite component has exact weighted value 2 * 1 - 1 = 1.
  W <- matrix(c(2, -1))
  for (native in c(FALSE, TRUE)) {
    options(RdistNN.use_compiled = native)
    for (sign in c(-1, 1)) {
      signal <- c(sign * .Machine$double.xmax, 0)
      expect_error(.dnn_estimate_from_signal(signal, W, n = 2L, arg = "Y"), "Y.*[Rr]escale")
      for (overflow_component in 1:2) {
        paired <- matrix(1, 2L, 2L)
        paired[, overflow_component] <- signal
        before <- paired
        out <- .affine_ratio_evaluate(paired, W, ranks = NULL, rw = NULL,
          active = FALSE, jackknife = "none", tol = 0)
        expect_identical(out$status$ratio_status, "component_overflow")
        expect_false(out$status$wald_available)
        expect_identical(unname(out$components[1L, overflow_component]), NA_real_)
        expect_identical(unname(out$components[1L, 3L - overflow_component]), 1)
        expect_identical(out$theta, NA_real_)
        expect_null(out$deletions)
        expect_identical(paired, before)
      }
    }
    expect_identical(unname(.dnn_estimate_from_signal(
      c(.Machine$double.xmax / 4, 0), W, n = 2L)), matrix(.Machine$double.xmax / 2))
  }
})

test_that("signed TDNN rejects mathematically overflowing fits and new-target predictions", {
  old <- options(RdistNN.use_compiled = FALSE)
  on.exit(options(old), add = TRUE)
  X <- matrix(1:4)
  signal <- c(.Machine$double.xmax, .Machine$double.xmax, 0, 0)
  # At zero the exact weights are (13/16, 1/4, -1/32, -1/32), so the
  # weighted signal is (17/16) * xmax. At five it is -xmax/16, finite.
  for (native in c(FALSE, TRUE)) {
    options(RdistNN.use_compiled = native)
    for (se in c(FALSE, TRUE)) {
      expect_error(np_tdnn(X = X, x0 = 0, Y = signal,
        s1 = 1L, s2 = 3L, method = "exact", se = se), "Y.*[Rr]escale")
    }
    fit <- np_tdnn(X = X, x0 = 5, Y = signal, s1 = 1L, s2 = 3L,
      method = "exact", se = FALSE)
    expect_equal(unname(coef(fit)) / .Machine$double.xmax, -1 / 16, tolerance = 1e-14)
    for (se.fit in c(FALSE, TRUE)) {
      expect_error(predict(fit, x0 = 0, se.fit = se.fit), "Y.*[Rr]escale")
    }
    expect_equal(predict(unserialize(serialize(fit, NULL)), x0 = 5), unname(coef(fit)))
  }
})

test_that("extreme constant mean fits return finite state or actionable [Rr]escale errors", {
  old <- options(RdistNN.use_compiled = FALSE)
  on.exit(options(old), add = TRUE)
  X <- matrix(seq_len(49L))
  for (native in c(FALSE, TRUE)) for (method in c("exact", "asymp")) for (se in c(FALSE, TRUE)) {
    options(RdistNN.use_compiled = native)
    for (family in c("np", "cate", "tdnn")) {
      # n=49 exposed accumulator overflow although the exact represented-weight
      # average is finite. Permit improved accumulation, never an invalid fit.
      extreme <- .estimate_range_outcome(.estimate_range_mean(
        family, X, rep(.Machine$double.xmax, 49L), 0, method, se),
        if (family == "cate") "pseudo_outcome.*[Rr]escale" else "Y.*[Rr]escale")
      if (!inherits(extreme, "error")) {
        expect_true(all(is.finite(coef(extreme))))
        expect_equal(unname(coef(extreme)) / .Machine$double.xmax, 1, tolerance = 1e-14)
        expect_identical(summary(extreme)$coefficients$estimate, unname(coef(extreme)))
      }
      signal <- rep(.Machine$double.xmax / 4, 49L)
      fit <- .estimate_range_mean(family, X, signal, 0, method, se)
      expect_equal(unname(coef(fit)) / signal[1L], 1, tolerance = 1e-14)
      if (se) expect_equal(fit$se, 0)
      restored <- unserialize(serialize(fit, NULL))
      expect_equal(.with_localization_policy_warnings(predict(restored, x0 = 50)), unname(coef(fit)))
    }
  }
})

test_that("extreme paired fits preserve finite components and status-governed unavailable values", {
  .with_partial_result_warning({
  old <- options(RdistNN.use_compiled = FALSE)
  on.exit(options(old), add = TRUE)
  X <- matrix(seq_len(49L))
  for (native in c(FALSE, TRUE)) for (method in c("exact", "asymp")) for (se in c(FALSE, TRUE)) {
    options(RdistNN.use_compiled = native)
    control <- .with_localization_policy_warnings(affine_ratio_dnn(X = X, x0 = 0,
      numerator = rep(1, 49L), denominator = rep(1, 49L), s = 1L, method = method, se = se))
    for (overflow_component in 1:2) {
      signal <- matrix(1, 49L, 2L)
      signal[, overflow_component] <- .Machine$double.xmax
      fit <- .with_localization_policy_warnings(affine_ratio_dnn(X = X, x0 = 0,
        numerator = signal[, 1L], denominator = signal[, 2L], k = 1L, s = 1L,
        method = method, se = se, return = c("W", "signals", "deletions", "pseudo_values")))
      expect_true(fit$status$ratio_status %in% c("ok", "component_overflow"))
      unavailable <- identical(fit$status$ratio_status, "component_overflow")
      expect_true(all(is.finite(fit$components) | is.na(fit$components)))
      if (unavailable) {
        expect_false(fit$status$wald_available)
        expect_true(is.na(fit$components[1L, overflow_component]))
        expect_true(is.na(unname(coef(fit))))
      } else {
        expect_true(all(is.finite(fit$components)))
        expect_true(is.finite(unname(coef(fit))))
        expect_equal(fit$components[1L, overflow_component] / .Machine$double.xmax,
          1, tolerance = 1e-14)
      }
      expect_identical(fit$components[1L, 3L - overflow_component], control$components[1L, 3L - overflow_component])
      .expect_identical_fast(unname(fit$signals), signal)
      .expect_identical_fast(unname(fit$.model$signal), signal)
      expect_identical(summary(fit)$coefficients$estimate, unname(coef(fit)))
      expect_identical(.with_localization_policy_warnings(diagnostics(fit))$ratio_status, fit$status$ratio_status)
      expect_identical(.with_localization_policy_warnings(plot_data_estimates(fit))$estimate, unname(coef(fit)))
      restored <- unserialize(serialize(fit, NULL))
      expect_identical(.with_localization_policy_warnings(predict(restored, x0 = 50)), unname(coef(fit)))
      prediction <- .with_localization_policy_warnings(predict(restored, x0 = 50, se.fit = TRUE, interval = "confidence"))
      expect_identical(unname(prediction$fit[, "fit"]), unname(coef(fit)))
      expect_true(all(is.finite(prediction$se.fit) | is.na(prediction$se.fit)))
      if (se) {
        expect_true(all(is.finite(fit$deletions)))
        if (unavailable) expect_true(all(is.na(fit$pseudo_values[, overflow_component, ])))
        expect_true(all(is.finite(fit$pseudo_values[, 3L - overflow_component, ])))
        if (!fit$status$wald_available) expect_true(all(is.na(confint(fit))))
        if (unavailable) expect_true(all(is.na(vcov(fit))))
      } else expect_null(fit$se)
    }
    finite <- .with_localization_policy_warnings(affine_ratio_dnn(X = X, x0 = 0,
      numerator = rep(.Machine$double.xmax / 4, 49L), denominator = rep(1, 49L),
      s = 1L, method = method, se = se))
    expect_identical(finite$status$ratio_status, "ok")
    expect_equal(unname(coef(finite)) / (.Machine$double.xmax / 4), 1, tolerance = 1e-14)
  }
  })
})

test_that("extreme candidates remain finite or carry actionable failure and ratio status", {
  .with_partial_result_warning({
  old <- options(RdistNN.use_compiled = FALSE)
  on.exit(options(old), add = TRUE)
  X <- matrix(1:49)
  signal <- c(rep(.Machine$double.xmax, 47L), 0, 0)
  targets <- matrix(c(0, 50))
  tuning <- new_tuning("estimate_range_candidates", eval = function(...) list(s = 1L),
    eval_grid = function(...) data.frame(s = c(1L, 3L)))
  for (native in c(FALSE, TRUE)) {
    options(RdistNN.use_compiled = native)
    for (family in c("np", "cate", "ratio")) {
      common <- list(X = X, X0 = targets, k = 1L, tuning = tuning,
        tuning_mode = "grid", se = FALSE, return = "W")
      fit <- .with_localization_policy_warnings(switch(family,
        np = do.call(np_dnn, c(common, list(Y = signal))),
        cate = do.call(cate_dnn, c(common, list(pseudo_outcome = signal))),
        ratio = do.call(affine_ratio_dnn, c(common,
          list(numerator = signal, denominator = rep(1, 49L))))))
      expect_true(all(is.finite(coef(fit))))
      status <- fit$tuning_grid$candidates$status[2L]
      expect_true(status %in% c("ok", "failed"))
      if (identical(status, "failed")) {
        expect_false(identical(family, "ratio"))
        expect_match(fit$tuning_grid$candidates$reason[2L], "[Rr]escale")
        expect_error(candidate_fit(fit, 2L), "[Rr]escale")
      } else {
        candidate <- .with_localization_policy_warnings(candidate_fit(fit, 2L))
        expect_length(coef(candidate), 2L)
        if (identical(family, "ratio")) {
          expect_true(all(candidate$status$ratio_status %in% c("ok", "component_overflow")))
          unavailable <- candidate$status$ratio_status == "component_overflow"
          expect_identical(is.na(unname(coef(candidate))), unavailable)
          expect_identical(is.na(candidate$components[, 1L]), unavailable)
          expect_true(all(is.finite(candidate$components[, 2L])))
          expect_true(all(is.finite(candidate$components) | is.na(candidate$components)))
          expect_identical(as.numeric(candidate$.model$signal[, 1L]), signal)
          expect_identical(diagnostics(candidate)$ratio_status, candidate$status$ratio_status)
        } else {
          expect_true(all(is.finite(coef(candidate))))
        }
        restored <- unserialize(serialize(candidate, NULL))
        expect_identical(.with_localization_policy_warnings(predict(restored, X0 = targets)),
          unname(coef(candidate)))
      }
      expect_true(all(is.finite(summary(fit)$coefficients$estimate)))
    }
  }
  })
})

test_that("extreme new-target mean predictions return finite values or [Rr]escale errors", {
  old <- options(RdistNN.use_compiled = FALSE)
  on.exit(options(old), add = TRUE)
  X <- matrix(seq_len(49L))
  for (native in c(FALSE, TRUE)) for (method in c("exact", "asymp")) {
    options(RdistNN.use_compiled = native)
    for (family in c("np", "cate", "tdnn")) {
      if (family == "tdnn") {
        weights <- weights_tdnn(X = X, x0 = 0, s1 = 2L, s2 = 4L, method = method)$W
        k <- 1L; s <- 1L
      } else {
        k <- if (method == "exact") 1L else 3L
        s <- if (method == "exact") 3L else 5L
        weights <- weights_dnn(X = X, x0 = 0, k = k, s = s, method = method)$W
      }
      signal <- ifelse(weights[, 1L] > 0, .Machine$double.xmax, 0)
      fit <- .estimate_range_mean(family, X, signal, 50, method, se = FALSE,
        k = k, s = s, s1 = 2L, s2 = 4L)
      expect_true(all(is.finite(coef(fit))))
      for (se.fit in c(FALSE, TRUE)) {
        prediction <- .estimate_range_outcome(predict(fit, x0 = 0, se.fit = se.fit))
        if (!inherits(prediction, "error")) {
          values <- if (se.fit) prediction$fit else prediction
          expect_length(values, 1L)
          expect_true(all(is.finite(values)))
        }
      }
      expect_equal(predict(fit, x0 = 50), unname(coef(fit)))
    }
  }
})
