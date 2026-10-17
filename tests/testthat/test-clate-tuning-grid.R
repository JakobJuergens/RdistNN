test_that("CLATE grid mode aligns candidate schemas and retains the default result", {
  a <- .clate_fit_arguments()
  a$s <- NULL
  a$tuning <- .clate_alignment_tuning()
  for (method in c("exact", "asymp")) for (engine in c("fast", "loo")) {
    scalar <- do.call(clate_dnn, c(a, list(method = method, jackknife = engine)))
    expect_warning(grid <- do.call(clate_dnn, c(a, list(method = method, jackknife = engine,
                                                        tuning_mode = "grid"))), "infeasible")
    for (field in c("theta", "se", "components", "component_vcov", "theta_vcov", "status", "x0")) {
      expect_equal(grid[[field]], scalar[[field]], tolerance = 1e-12)
    }
    payload <- grid$tuning_grid
    expect_named(payload, c("candidates", "estimates", "diagnostics", "meta"))
    expect_true(all(c("candidate_id", "computed_key", "is_default", "feasible", "reason", "status", "s") %in% names(payload$candidates)))
    expect_true(all(c("candidate_id", "target_id", "component", "estimate", "se", "status", "reason") %in% names(payload$estimates)))
    expect_equal(payload$meta$n_unique, 3L)
    expect_equal(payload$meta$n_skipped, 1L)
    for (i in which(payload$candidates$feasible)) {
      row <- payload$candidates[i, ]
      args <- a
      args$tuning <- NULL
      args$s <- row$s
      fit <- do.call(clate_dnn, c(args, list(method = method, jackknife = engine)))
      tab <- payload$estimates[payload$estimates$candidate_id == row$candidate_id, ]
      expect_equal(tab$estimate, fit$theta, tolerance = 1e-12)
      expect_equal(tab$se, fit$se, tolerance = 1e-12)
    }
    skipped <- payload$estimates$status == "infeasible"
    expect_true(all(is.na(payload$estimates$estimate[skipped])))
    expect_true(all(is.na(payload$estimates$se[skipped])))
    expect_s3_class(diagnostics(grid, candidates = "all"), "dnn_diagnostics")
    expect_s3_class(plot_data_estimates(grid, candidates = "all"), "data.frame")
  }
})

test_that("grid evaluation freezes geometry and keeps valid point estimates without covariance", {
  calls <- 0L
  distance <- new_distance("counted_clate", eval = function(X, x0, ctx = NULL) {
    calls <<- calls + 1L
    abs(X[, 1L] - x0[1L]) + seq_len(nrow(X)) / nrow(X)^2
  })
  a <- .clate_fit_arguments()
  a$s <- NULL
  a$tuning <- .clate_alignment_tuning()
  expect_warning(fit <- do.call(clate_dnn, c(a, list(distance = distance,
    tuning_mode = "grid", method = "asymp"))), "infeasible")
  expect_equal(calls, nrow(a$X0))
  no_cov <- do.call(clate_dnn, c(a, list(tuning_mode = "grid", se = FALSE)))
  expect_null(no_cov$component_vcov)
  expect_true(all(is.na(no_cov$tuning_grid$estimates$se)))
  expect_true(all(is.finite(no_cov$tuning_grid$estimates$estimate)))
  expect_error(.ratio_fixture(tuning_mode = "grid"), "tuning_mode")
})
