# Recovery requires an explicitly declared optional failure; ordinary host
# conditions must retain their normal control flow.
.reporting_grid <- function(callback) {
  candidates <- data.frame(candidate_id = c("default", "optional"),
    computed_key = c("s_2", "s_3"), s = c(2L, 3L),
    is_default = c(TRUE, FALSE), feasible = TRUE, reason = NA_character_, status = "ok")
  RdistNN:::.dnn_grid_eval(candidates, matrix(1:5, 5, 1),
    weight_spec_fn = function(row) rank_weight_dnn_exact(k = 1L, s = row$s),
    estimate_fn = function(W, row) {
      if (!row$is_default) callback()
      data.frame(target_id = 1L, component = NA_character_, estimate = sum(W))
    })
}

.reporting_method_failure <- function(message) {
  stop(structure(list(message = message, call = NULL),
    class = c("rdistnn_simulation_method_failure", "rdistnn_error", "error", "condition")))
}

.reporting_warnings <- function(expr) {
  seen <- list()
  value <- withCallingHandlers(expr, warning = function(w) {
    seen[[length(seen) + 1L]] <<- w
    invokeRestart("muffleWarning")
  })
  list(value = value, warnings = seen)
}

test_that("optional grid callbacks propagate unexpected errors and promoted warnings", {
  expect_error(.reporting_grid(function() stop("unexpected callback")), "unexpected callback")
  withr::local_options(warn = 2)
  expect_error(.reporting_grid(function() warning("callback warning")), "callback warning")
  # Calling handlers get first refusal even when the caller promotes warnings.
  handled <- .reporting_warnings(.reporting_grid(function() warning("handled warning")))
  expect_identical(handled$value$candidates$status, c("ok", "ok"))
  expect_length(handled$warnings, 1L)
})

test_that("only declared optional candidate failures retain failed rows", {
  failure <- structure(list(message = "declared candidate failure", call = NULL),
    class = c("rdistnn_candidate_failure", "rdistnn_error", "error", "condition"))
  out <- .reporting_warnings(.reporting_grid(function() stop(failure)))
  expect_identical(out$value$candidates$status, c("ok", "failed"))
  expect_equal(nrow(out$value$estimates), 2L)
  expect_true(is.na(out$value$estimates$estimate[[2L]]))
  expect_length(out$warnings, 1L)
  expect_s3_class(out$warnings[[1L]], "rdistnn_partial_result_warning")
})

test_that("dashboard propagates explicit invalid options and malformed retained weights", {
  skip_if_not_installed("ggplot2")
  fit <- np_dnn(X = matrix(1:10), x0 = 0, Y = 1:10, s = 3L)
  expect_error(plot_dashboard(fit, style = "not_a_style"), "one of")
  fit$W[1, 1] <- NA_real_
  expect_error(plot_dashboard(fit), class = "rdistnn_state_error")
})

test_that("dashboard omits only missing optional extras and retains list fallback", {
  skip_if_not_installed("ggplot2")
  fit <- np_dnn(X = matrix(1:10), x0 = 0, Y = 1:10, s = 3L, return = character())
  expect_error(RdistNN:::.dnn_diag_require_W(fit), class = "rdistnn_extra_unavailable")
  original <- RdistNN:::.plot_pkg_installed
  local_mocked_bindings(.plot_pkg_installed = function(pkg) {
    if (pkg == "patchwork") FALSE else original(pkg)
  }, .package = "RdistNN")
  expect_message(out <- plot_dashboard(fit), "named list")
  expect_named(out, "estimates")
})

test_that("ratio partial output keeps every target and one bounded status summary", {
  fit_call <- function(se) affine_ratio_dnn(X = matrix(1:8), X0 = matrix(seq(-2, -1, length.out = 100)),
    numerator = 1:8, denominator = rep(0, 8), s = 2L, se = se)
  out <- .reporting_warnings(fit_call(FALSE))
  expect_length(out$warnings, 1L)
  expect_s3_class(out$warnings[[1L]], "rdistnn_partial_result_warning")
  expect_lt(nchar(conditionMessage(out$warnings[[1L]])), 800L)
  expect_lte(length(out$warnings[[1L]]$examples), 5L)
  expect_equal(out$warnings[[1L]]$affected, 100L)
  expect_match(conditionMessage(out$warnings[[1L]]), "status")
  expect_equal(nrow(out$value$status), 100L)
  expect_true(all(out$value$status$ratio_status == "zero_denominator"))
  expect_true(all(is.finite(out$value$components[, 1L])))
  expect_true(all(is.na(out$value$theta)))
  again <- .reporting_warnings(fit_call(FALSE))
  expect_length(again$warnings, 1L)
})

test_that("defined negative ratios and deliberately disabled covariance do not warn", {
  for (warn in c(-1, 0, 1, 2)) {
    withr::local_options(warn = warn)
    out <- .reporting_warnings(affine_ratio_dnn(X = matrix(1:8), x0 = 0,
      numerator = 1:8, denominator = rep(-1, 8), s = 2L, se = FALSE))
    expect_length(out$warnings, 0L)
    expect_true(out$value$status$negative_denominator)
    expect_true(is.finite(out$value$theta))
  }
})

test_that("grid does not recover generic numerical, resource, interrupt or internal conditions", {
  conditions <- list(
    structure(list(message = "numerical", call = NULL),
      class = c("rdistnn_numerical_error", "rdistnn_error", "error", "condition")),
    structure(list(message = "allocation", call = NULL),
      class = c("resourceError", "error", "condition")),
    structure(list(message = "interrupt", call = NULL),
      class = c("interrupt", "condition")),
    structure(list(message = "invariant", call = NULL),
      class = c("rdistnn_internal_error", "rdistnn_error", "error", "condition")))
  for (cnd in conditions) {
    caught <- tryCatch(.reporting_grid(function() stop(cnd)), condition = identity)
    expect_identical(caught, cnd)
  }
})

test_that("grid callbacks do not change options or RNG while preserving warning handling", {
  set.seed(551)
  rng <- .Random.seed
  opts <- options()
  for (warn in c(-1, 0, 1, 2)) {
    withr::local_options(warn = warn)
    out <- .reporting_warnings(.reporting_grid(function() warning("visible warning")))
    expect_length(out$warnings, 1L)
    expect_identical(out$value$candidates$status, c("ok", "ok"))
    expect_equal(getOption("warn"), warn)
  }
  expect_identical(.Random.seed, rng)
  expect_identical(options()[setdiff(names(opts), "warn")], opts[setdiff(names(opts), "warn")])
})

test_that("dashboard unexpected render failures and promoted warnings propagate", {
  skip_if_not_installed("ggplot2")
  fit <- np_dnn(X = matrix(1:10), x0 = 0, Y = 1:10, s = 3L)
  local_mocked_bindings(plot_estimates = function(...) stop("renderer failure"), .package = "RdistNN")
  expect_error(plot_dashboard(fit), "renderer failure")
})

test_that("dashboard calling handlers can muffle warnings under warn two", {
  skip_if_not_installed("ggplot2")
  fit <- np_dnn(X = matrix(1:10), x0 = 0, Y = 1:10, s = 3L)
  original <- RdistNN:::plot_estimates
  local_mocked_bindings(plot_estimates = function(...) {
    warning("panel qualification")
    original(...)
  }, .package = "RdistNN")
  withr::local_options(warn = 2)
  expect_error(plot_dashboard(fit), "panel qualification")
  handled <- .reporting_warnings(plot_dashboard(fit))
  expect_length(handled$warnings, 1L)
})

test_that("malformed optional grid callback output is a callback contract error", {
  candidates <- data.frame(candidate_id = c("default", "optional"),
    computed_key = c("s_2", "s_3"), s = c(2L, 3L),
    is_default = c(TRUE, FALSE), feasible = TRUE, reason = NA_character_, status = "ok")
  expect_error(RdistNN:::.dnn_grid_eval(candidates, matrix(1:5, 5, 1),
    weight_spec_fn = function(row) rank_weight_dnn_exact(k = 1L, s = row$s),
    estimate_fn = function(W, row) {
      if (row$is_default) data.frame(target_id = 1L, component = NA_character_, estimate = 1) else NULL
    }), class = "rdistnn_callback_error")
})

test_that("simulation recovers declared failures while internal failures propagate", {
  original <- RdistNN:::.get_sim_method
  local_mocked_bindings(.get_sim_method = function(name) {
    if (name == "causal_forest") return(function(...) .reporting_method_failure("declared fit failure"))
    if (name == "dnn_oracle") return(function(...) RdistNN:::.dnn_abort("internal_invariant", "package invariant"))
    original(name)
  }, .package = "RdistNN")
  out <- .reporting_warnings(sim_cate_one_rep(n = 30L, kx = 2L, k = 1L, s = 4L,
    design = "dgp1_linear", methods = "causal_forest", seed = 27L))
  expect_identical(out$value$status, "failed")
  expect_length(out$warnings, 1L)
  expect_s3_class(out$warnings[[1L]], "rdistnn_method_failed_warning")
  expect_error(sim_cate_one_rep(n = 30L, kx = 2L, k = 1L, s = 4L,
    design = "dgp1_linear", methods = "dnn_oracle", seed = 27L), "package invariant")
})

test_that("simulation promoted warnings propagate and outer calling handlers retain control", {
  original <- RdistNN:::.get_sim_method
  local_mocked_bindings(.get_sim_method = function(name) {
    if (name == "causal_forest") return(function(...) {
      warning("external qualification")
      RdistNN:::.sim_cate_failure_rows("causal_forest", 30L, 2L, 1L, 4L, NULL, "default", .05)
    })
    original(name)
  }, .package = "RdistNN")
  withr::local_options(warn = 2)
  call <- function() sim_cate_one_rep(n = 30L, kx = 2L, k = 1L, s = 4L,
    design = "dgp1_linear", methods = "causal_forest", seed = 28L)
  expect_error(call(), "external qualification")
  out <- .reporting_warnings(call())
  expect_identical(conditionMessage(out$warnings[[1L]]), "external qualification")
  expect_identical(out$value$method, "causal_forest")
})

test_that("simulation batch propagates DGP errors and preserves caller RNG", {
  set.seed(556)
  before <- .Random.seed
  expect_error(sim_cate_run(designs = "bad_design", n_vals = 30L, kx_vals = 2L,
    n_reps = 1L, master_seed = 57L, methods = "dnn_oracle", verbose = FALSE), "design")
  expect_identical(.Random.seed, before)
})

test_that("IV nuisance adapters propagate promoted warnings and internal errors", {
  original <- RdistNN:::.fit_sim_iv_nuisances_crossfit
  local_mocked_bindings(.fit_sim_iv_nuisances_crossfit = function(...) {
    warning("IV nuisance qualification")
    original(...)
  }, .package = "RdistNN")
  withr::local_options(warn = 2)
  call <- function() sim_iv_one_rep(n = 40L, kx = 2L, q = 1L, k = 1L, s = 4L,
    design = "dgp1_iv_justid_linear", methods = "dnn_iv_crossfit", seed = 29L, n_folds = 2L)
  expect_error(call(), "IV nuisance qualification")
  handled <- .reporting_warnings(call())
  expect_true(any(vapply(handled$warnings, conditionMessage, character(1)) == "IV nuisance qualification"))
  expect_equal(nrow(handled$value), 1L)
})

test_that("tuning callback errors in simulation metadata propagate", {
  tuning <- new_tuning(name = "broken callback", eval = function(...) stop("tuning invariant"))
  expect_error(RdistNN:::.sim_cate_default_lambda(tuning, 30L, 2L, 1L), "tuning invariant")
})

test_that("simulation batches aggregate declared failures and retain reasons", {
  original <- RdistNN:::.get_sim_method
  local_mocked_bindings(.get_sim_method = function(name) {
    if (name == "causal_forest") return(function(...) .reporting_method_failure("declared method failure"))
    original(name)
  }, .package = "RdistNN")
  set.seed(557)
  rng <- .Random.seed
  out <- .reporting_warnings(sim_cate_run(designs = "dgp1_linear", n_vals = 30L,
    kx_vals = 2L, n_reps = 3L, master_seed = 58L,
    methods = c("dnn_oracle", "causal_forest"), verbose = FALSE))
  expect_length(out$warnings, 1L)
  expect_s3_class(out$warnings[[1L]], "rdistnn_method_failed_warning")
  expect_equal(out$warnings[[1L]]$affected, 3L)
  expect_equal(sum(out$value$status == "failed"), 3L)
  expect_equal(nrow(attr(out$value, "method_failures")), 3L)
  expect_true(all(attr(out$value, "method_failures")$reason == "declared method failure"))
  expect_identical(.Random.seed, rng)
})

test_that("ratio inference consumers summarize only requested unavailable targets", {
  fit <- .reporting_warnings(affine_ratio_dnn(X = matrix(1:4), X0 = matrix(c(0, 5)),
    numerator = c(2, 1, 4, 3), denominator = c(0, 0, -1, 2), s = 3L))$value
  bad <- .reporting_warnings(confint(fit, parm = 1L))
  expect_length(bad$warnings, 1L)
  expect_s3_class(bad$warnings[[1L]], "rdistnn_partial_result_warning")
  expect_true(all(is.na(bad$value)))
  good <- .reporting_warnings(confint(fit, parm = 2L))
  expect_length(good$warnings, 0L)
  expect_true(all(is.finite(good$value)))
  shown <- .reporting_warnings(summary(fit))
  expect_length(shown$warnings, 1L)
  expect_equal(nrow(shown$value$coefficients), 2L)
})

test_that("nested external IV regression recovery does not erase warning promotion", {
  original <- stats::lm
  local_mocked_bindings(lm = function(...) {
    warning("external regression qualification")
    original(...)
  }, .package = "stats")
  withr::local_options(warn = 2)
  call <- function() sim_iv_one_rep(n = 40L, kx = 2L, q = 1L, k = 1L, s = 4L,
    design = "dgp1_iv_justid_linear", methods = "2sls_global", seed = 30L)
  expect_error(call(), "external regression qualification")
  out <- .reporting_warnings(call())
  expect_true(any(vapply(out$warnings, conditionMessage, character(1)) == "external regression qualification"))
  expect_true(is.finite(out$value$theta_hat))
})

test_that("built-in mean adapters declare only their known candidate range failures", {
  original <- RdistNN:::.dnn_estimate_from_signal
  parent <- NULL
  local_mocked_bindings(.dnn_estimate_from_signal = function(signal, W, n, arg = "signal") {
    if (W[1L, 1L] > .4) {
      parent <<- tryCatch(RdistNN:::.dnn_abort("numerical_range", "known candidate mean range"),
        rdistnn_numerical_error = identity)
      stop(parent)
    }
    original(signal, W, n, arg)
  }, .package = "RdistNN")
  tun <- new_tuning("range candidate", eval = function(...) list(s = 2L),
    eval_grid = function(...) data.frame(s = c(2L, 3L)))
  for (fit in list(
    function() np_dnn(X = matrix(1:6), x0 = 0, Y = 1:6, tuning = tun, tuning_mode = "grid", se = FALSE),
    function() cate_dnn(X = matrix(1:6), x0 = 0, pseudo_outcome = 1:6, tuning = tun, tuning_mode = "grid", se = FALSE))) {
    out <- .reporting_warnings(fit())
    expect_identical(out$value$tuning_grid$candidates$status, c("ok", "failed"))
    expect_length(out$warnings, 1L)
    expect_s3_class(out$warnings[[1L]], "rdistnn_partial_result_warning")
  }
})

test_that("grid rank-weight provider violations retain callback identity", {
  candidates <- data.frame(candidate_id = "default", computed_key = "s_2", s = 2L,
    is_default = TRUE, feasible = TRUE, reason = NA_character_, status = "ok")
  expect_error(RdistNN:::.dnn_grid_eval(candidates, matrix(1:5, 5, 1),
    weight_spec_fn = function(row) NULL,
    estimate_fn = function(W, row) data.frame(target_id = 1L, component = NA_character_, estimate = 1)),
    class = "rdistnn_callback_error")
})

test_that("absent optional GRF methods produce one complete batch failure summary", {
  local_mocked_bindings(.sim_cate_grf_available = function() FALSE, .package = "RdistNN")
  out <- .reporting_warnings(sim_cate_run(designs = "dgp1_linear", n_vals = 30L,
    kx_vals = 2L, n_reps = 3L, master_seed = 59L,
    methods = c("causal_forest", "llf"), verbose = FALSE))
  expect_length(out$warnings, 1L)
  expect_s3_class(out$warnings[[1L]], "rdistnn_method_failed_warning")
  expect_equal(out$warnings[[1L]]$affected, 6L)
  expect_equal(sum(out$value$status == "failed"), 6L)
  expect_equal(nrow(attr(out$value, "method_failures")), 6L)
})


test_that("opaque external method errors propagate without failure replacement", {
  original <- RdistNN:::.get_sim_method
  local_mocked_bindings(.get_sim_method = function(name) {
    if (name == "causal_forest") return(function(...) stop("opaque external failure"))
    original(name)
  }, .package = "RdistNN")
  expect_error(sim_cate_one_rep(n = 30L, kx = 2L, k = 1L, s = 4L,
    design = "dgp1_linear", methods = "causal_forest", seed = 31L), "opaque external failure")
})

test_that("opaque nested IV regression errors propagate without NA substitution", {
  local_mocked_bindings(lm = function(...) stop("opaque regression failure"), .package = "stats")
  expect_error(sim_iv_one_rep(n = 40L, kx = 2L, q = 1L, k = 1L, s = 4L,
    design = "dgp1_iv_justid_linear", methods = "2sls_global", seed = 32L),
    "opaque regression failure")
})

test_that("simulation linear-system checks preserve solve's existing numerical cutoff", {
  safe <- RdistNN:::.sim_iv_safe_solve
  eps <- .Machine$double.eps
  for (a in c(0, eps / 2, eps, 2 * eps, .5)) {
    A <- diag(c(1, a))
    reference <- tryCatch(solve(A), error = function(e) NULL)
    expect_equal(safe(A), reference)
    reference_b <- tryCatch(solve(A, c(1, 1)), error = function(e) NULL)
    expect_equal(safe(A, c(1, 1)), reference_b)
  }
  for (bad in c(NA_real_, NaN, Inf, -Inf)) {
    A <- diag(c(1, bad))
    reference <- tryCatch(solve(A), error = identity)
    actual <- tryCatch(safe(A), error = identity)
    if (inherits(reference, "error")) {
      expect_s3_class(actual, "error")
      expect_identical(conditionMessage(actual), conditionMessage(reference))
    } else expect_equal(actual, reference)
  }
})

test_that("unexpected simulation solver failures propagate", {
  local_mocked_bindings(solve = function(...) stop("opaque solver failure"), .package = "base")
  expect_error(RdistNN:::.sim_iv_safe_solve(diag(2)), "opaque solver failure")
})


test_that("unexpected simulation condition-number failures propagate", {
  local_mocked_bindings(rcond = function(...) stop("opaque condition-number failure"), .package = "base")
  expect_error(RdistNN:::.sim_iv_safe_solve(diag(2)), "opaque condition-number failure")
})

test_that("recognized oracle IV numerical failures retain reasons and aggregate", {
  call <- function(s = 40L) sim_iv_one_rep(n = 40L, kx = 2L, q = 2L, k = 1L,
    s = s, design = "dgp2_iv_homog_overid", methods = "iv_oracle", seed = 10L)
  failed <- .reporting_warnings(call())
  expect_true(is.na(failed$value$theta_hat))
  expect_true(is.na(failed$value$se))
  expect_length(failed$warnings, 1L)
  expect_s3_class(failed$warnings[[1L]], "rdistnn_method_failed_warning")
  detail <- attr(failed$value, "method_failures")
  expect_equal(nrow(detail), 1L)
  expect_identical(detail$method, "iv_oracle")
  expect_true(nzchar(detail$reason))
  expect_lte(nchar(detail$reason), 240L)
  finite <- .reporting_warnings(call(s = 5L))
  expect_true(is.finite(finite$value$theta_hat))
  expect_length(finite$warnings, 0L)
  expect_null(attr(finite$value, "method_failures"))
  batch <- .reporting_warnings(sim_iv_run(designs = "dgp2_iv_homog_overid",
    n_vals = 40L, kx_vals = 2L, q_vals = 2L, k = 1L, n_reps = 3L,
    s_fn = function(n, kx) n, methods = "iv_oracle", master_seed = 10L, verbose = FALSE))
  expect_length(batch$warnings, 1L)
  expect_equal(batch$warnings[[1L]]$affected, 3L)
  expect_true(all(is.na(batch$value$theta_hat)))
  expect_equal(nrow(attr(batch$value, "method_failures")), 3L)
  expect_false(anyNA(attr(batch$value, "method_failures")$seed))
})
