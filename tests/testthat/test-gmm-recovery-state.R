.gmm_recovery_fit <- function(model, se = TRUE, ...) {
  gmm_dnn(X = matrix(seq_len(8L)), X0 = matrix(c(1.25, 6.25)),
    model_data = data.frame(y = c(1, 3, 2, 5, 4, 7, 6, 8)),
    model = model, start = NULL, s = 3L, se = se, ...)
}

.gmm_recovery_model <- function() {
  new_gmm_model("recovery_state",
    setup = function(sample, n, ...) list(y = sample$y),
    moment_obs = function(sample, theta, ctx) ctx$y - theta[1L],
    jacobian_obs = function(sample, theta, ctx) rep(-1, length(ctx$y)),
    start = function(sample, ctx, w, x0) 0,
    coef_names = "location", moment_names = "mean")
}

.gmm_fault_model <- function(stage, signal) {
  model <- .gmm_recovery_model()
  if (stage == "start") model$start <- function(...) signal()
  if (stage %in% c("probe", "objective")) {
    model$moment_obs <- function(sample, theta, ctx) {
      if (stage == "probe" || theta[1L] != 0) return(signal())
      ctx$y - theta[1L]
    }
  }
  if (stage == "jacobian") model$jacobian_obs <- function(...) signal()
  model
}

test_that("GMM propagates opaque callback and resource failures at every boundary", {
  for (stage in c("start", "probe", "objective", "jacobian")) {
    for (class in c("opaque_callback_error", "simulated_resource_error")) {
      original <- structure(list(message = "callback sentinel", call = NULL),
        class = c(class, "error", "condition"))
      model <- .gmm_fault_model(stage, function() stop(original))
      observed <- tryCatch(.gmm_recovery_fit(model), error = identity)
      expect_identical(observed, original, info = paste(stage, class))
    }
  }
})

test_that("GMM preserves callback warnings promoted by the caller", {
  old <- options(warn = 2)
  on.exit(options(old), add = TRUE)
  for (stage in c("start", "probe", "objective", "jacobian")) {
    model <- .gmm_fault_model(stage, function() warning("callback warning sentinel"))
    observed <- tryCatch(.gmm_recovery_fit(model), error = identity)
    expect_s3_class(observed, "error")
    if (inherits(observed, "error")) {
      expect_match(conditionMessage(observed), "callback warning sentinel", info = stage)
    }
  }
})

test_that("GMM callback warnings keep caller handlers and warning options", {
  old <- options("warn")
  on.exit(options(old), add = TRUE)
  original <- structure(list(message = "callback notice", call = NULL),
    class = c("callback_notice", "warning", "condition"))
  model <- .gmm_recovery_model()
  model$jacobian_obs <- function(sample, theta, ctx) {
    warning(original)
    rep(-1, length(ctx$y))
  }
  for (level in c(-1L, 0L, 1L, 2L)) {
    options(warn = level)
    seen <- list()
    fit <- withCallingHandlers(.gmm_recovery_fit(model), callback_notice = function(w) {
      seen[[length(seen) + 1L]] <<- w
      invokeRestart("muffleWarning")
    })
    expect_identical(seen, list(original, original))
    expect_true(all(is.finite(coef(fit))))
    expect_identical(getOption("warn"), level)
  }
})

test_that("GMM propagates interrupts without replacing them by missing targets", {
  original <- structure(list(message = "interrupt sentinel", call = NULL),
    class = c("interrupt", "condition"))
  for (stage in c("start", "probe", "objective", "jacobian")) {
    model <- .gmm_fault_model(stage, function() stop(original))
    observed <- tryCatch(.gmm_recovery_fit(model), interrupt = identity)
    expect_identical(observed, original, info = stage)
  }
})

test_that("malformed trial moment returns stop GMM with a callback condition", {
  model <- .gmm_fault_model("objective", function() "invalid moment")
  observed <- tryCatch(.gmm_recovery_fit(model), error = identity)
  expect_s3_class(observed, "rdistnn_callback_error")
  expect_identical(observed$code, "callback_invalid")
  expect_identical(observed$argument, "model$moment_obs()")
})

test_that("GMM penalties apply only to arithmetic overflow after finite moments", {
  model <- .gmm_recovery_model()
  model$moment_obs <- function(sample, theta, ctx) rep(1e200, 8L)
  value <- RdistNN:::.gmm_objective(0, model, NULL, NULL, rep(1 / 8, 8L), 8L)
  expect_identical(value, .Machine$double.xmax)
  model$moment_obs <- function(...) stop("opaque objective failure")
  expect_error(RdistNN:::.gmm_objective(0, model, NULL, NULL, rep(1 / 8, 8L), 8L),
    "opaque objective failure")
})

test_that("owned GMM numerical target failures retain rows and aggregate once", {
  model <- .gmm_recovery_model()
  fit_targets <- function(targets) gmm_dnn(X = matrix(seq_len(8L)),
    X0 = matrix(targets), model_data = data.frame(y = c(rep(0, 6L), 1, 2)),
    model = model, s = 3L, omega_ridge = 0,
    return = c("g_obs", "omega", "jacobian"))
  capture_partial <- function(expr) {
    seen <- list()
    value <- withCallingHandlers(expr, rdistnn_partial_result_warning = function(w) {
      seen[[length(seen) + 1L]] <<- w
      invokeRestart("muffleWarning")
    })
    list(value = value, warnings = seen)
  }
  for (targets in list(c(1.25, 8.25), c(1.25, 1.125))) {
    observed <- capture_partial(fit_targets(targets))
    fit <- observed$value
    failed <- if (targets[2L] == 8.25) c(TRUE, FALSE) else c(TRUE, TRUE)
    expect_identical(unname(is.na(coef(fit)[1L, ])), failed)
    expect_length(observed$warnings, 1L)
    expect_identical(observed$warnings[[1L]]$code, "partial_result")
    expect_equal(observed$warnings[[1L]]$affected, sum(failed))
    expect_identical(observed$warnings[[1L]]$examples, as.character(which(failed)))
    expect_match(conditionMessage(observed$warnings[[1L]]), "[Rr]escale.*convergence")
    expect_true(all(grepl("Could not invert", fit$convergence$message[failed])))
    expect_equal(nrow(summary(fit)$coefficients), 2L)
    restored <- unserialize(serialize(fit, NULL))
    expect_identical(coef(restored), fit$theta)
    expect_identical(vcov(restored), fit$vcov)
    expect_identical(capture_partial(fit_targets(targets))$warnings[[1L]]$affected,
      observed$warnings[[1L]]$affected)
  }
})

test_that("GMM singular-system preflight matches the existing solve cutoff", {
  safe <- RdistNN:::.gmm_safe_solve
  eps <- .Machine$double.eps
  for (a in c(0, eps / 2, eps, 2 * eps, .5)) {
    A <- diag(c(1, a))
    expect_equal(safe(A), tryCatch(solve(A), error = function(e) NULL))
    expect_equal(safe(A, c(1, 1)), tryCatch(solve(A, c(1, 1)), error = function(e) NULL))
  }
})

test_that("unexpected GMM linear algebra failures propagate", {
  local_mocked_bindings(solve = function(...) stop("opaque solver failure"), .package = "base")
  expect_error(RdistNN:::.gmm_safe_solve(diag(2L)), "opaque solver failure")
})

test_that("unexpected GMM eigenvalue failures propagate", {
  local_mocked_bindings(eigen = function(...) stop("opaque eigenvalue failure"), .package = "base")
  expect_error(RdistNN:::.gmm_condition_number(diag(2L)), "opaque eigenvalue failure")
})

test_that("GMM methods reject corrupted current arrays with refit guidance", {
  fit <- .gmm_recovery_fit(.gmm_recovery_model())
  corruptions <- list(
    theta = "this is corrupted",
    theta_first = matrix(0, 2L, 2L),
    se = matrix(-1, 1L, 2L),
    vcov = list(matrix("bad"), NULL),
    J = "bad", df = c(1.5, 1.5), pval_J = c(2, 2), neff = c(-1, -1),
    convergence = data.frame(target = "wrong"))
  consumers <- list(coef = coef, vcov = vcov, summary = summary,
    confint = confint, print = function(x) capture.output(print(x)),
    plot_data_estimates = plot_data_estimates,
    plot_weights = plot_weights, plot_estimates = plot_estimates,
    plot = function(x) plot(x))
  for (field in names(corruptions)) {
    bad <- fit
    bad[[field]] <- corruptions[[field]]
    for (consumer in names(consumers)) {
      observed <- tryCatch(consumers[[consumer]](bad), error = identity)
      expect_s3_class(observed, "rdistnn_state_error")
      if (inherits(observed, "rdistnn_state_error")) {
        expect_identical(observed$code, "state_invalid")
        expect_match(conditionMessage(observed), "[Rr]efit")
      }
    }
  }
})

test_that("GMM current-state validation rejects absent and misaligned arrays", {
  fit <- .gmm_recovery_fit(.gmm_recovery_model(), return = c("omega", "jacobian", "g_obs"))
  for (field in c("theta", "theta_first", "se", "vcov", "J", "df", "neff", "convergence")) {
    bad <- fit
    bad[[field]] <- NULL
    expect_error(coef(bad), class = "rdistnn_state_error")
  }
  for (field in c("theta", "theta_first", "se")) {
    bad <- fit
    colnames(bad[[field]]) <- rev(colnames(bad[[field]]))
    expect_error(coef(bad), class = "rdistnn_state_error")
  }
  for (field in c("omega", "jacobian", "g_obs")) {
    bad <- fit
    bad[[field]][1L] <- list(matrix("corrupt"))
    expect_error(coef(bad), class = "rdistnn_state_error")
  }
  for (field in c("q", "p_theta", "se", "coef_names", "moment_names")) {
    bad <- fit
    bad$meta[[field]] <- NULL
    expect_error(summary(bad), class = "rdistnn_state_error")
  }
})

test_that("current GMM serialization preserves public arrays and missing uncertainty", {
  for (se in c(TRUE, FALSE)) {
    fit <- .gmm_recovery_fit(.gmm_recovery_model(), se = se)
    restored <- unserialize(serialize(fit, NULL))
    expect_identical(coef(restored), fit$theta)
    expect_identical(summary(restored)$coefficients, summary(fit)$coefficients)
    expect_identical(plot_data_estimates(restored), plot_data_estimates(fit))
    if (se) expect_identical(vcov(restored), vcov(fit)) else {
      expect_error(vcov(restored), class = "rdistnn_extra_unavailable")
      expect_error(confint(restored), class = "rdistnn_extra_unavailable")
    }
  }
})

.gmm_fd_arithmetic_model <- function() {
  new_gmm_model("finite_moment_derivative_range",
    setup = function(sample, n, ...) NULL,
    moment_obs = function(sample, theta, ctx) rep(1e308 * tanh(theta[1L] / 1e-7), 8L),
    start = function(...) 0, coef_names = "location", moment_names = "mean")
}

test_that("optional GMM Jacobians preserve arithmetic-unavailable current fits", {
  old <- options(RdistNN.use_compiled = FALSE)
  on.exit(options(old), add = TRUE)
  model <- .gmm_fd_arithmetic_model()
  expect_true(all(is.finite(model$moment_obs(NULL, 1e-5, NULL))))
  expect_true(all(is.finite(model$moment_obs(NULL, -1e-5, NULL))))
  for (se in c(TRUE, FALSE)) for (keep in c(TRUE, FALSE)) {
    notices <- list()
    fit <- withCallingHandlers(gmm_dnn(X = matrix(seq_len(8L)), x0 = 0,
      model = model, s = 3L, se = se, jacobian = "fd",
      return = if (keep) "jacobian" else character()),
      rdistnn_partial_result_warning = function(w) {
        notices[[length(notices) + 1L]] <<- w
        invokeRestart("muffleWarning")
      })
    expect_identical(as.numeric(coef(fit)), 0)
    expect_length(notices, if (se) 1L else 0L)
    if (se) {
      expect_true(all(is.na(fit$se)))
      expect_null(vcov(fit))
      expect_match(conditionMessage(notices[[1L]]), "information.*numeric range.*[Rr]escale")
    } else {
      expect_null(fit$se)
      expect_error(vcov(fit), class = "rdistnn_extra_unavailable")
    }
    if (keep) expect_identical(as.numeric(fit$jacobian[[1L]]), Inf) else {
      expect_false("jacobian" %in% names(fit))
    }
    restored <- unserialize(serialize(fit, NULL))
    expect_identical(coef(restored), fit$theta)
    expect_identical(summary(restored)$coefficients, summary(fit)$coefficients)
    expect_identical(plot_data_estimates(restored), plot_data_estimates(fit))
  }
})

test_that("unavailable GMM diagnostics do not license nonfinite callback returns", {
  for (se in c(TRUE, FALSE)) {
    model <- .gmm_fd_arithmetic_model()
    model$jacobian_obs <- function(...) rep(Inf, 8L)
    error <- expect_error(gmm_dnn(X = matrix(seq_len(8L)), x0 = 0,
      model = model, s = 3L, se = se, jacobian = "analytic", return = "jacobian"),
      class = "rdistnn_callback_error")
    expect_identical(error$code, "callback_invalid")
    expect_identical(error$argument, "model$jacobian_obs()")
  }
  fit <- .gmm_recovery_fit(.gmm_recovery_model(), return = "jacobian")
  expect_true(all(is.finite(fit$se)))
  fit$jacobian[[1L]][1L] <- Inf
  expect_error(coef(fit), class = "rdistnn_state_error")
})
