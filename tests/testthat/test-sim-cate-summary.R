# Simulation aggregation, nuisance routes, and failed-replicate counts.
test_that("summarise_sim_results: correct aggregation with method column", {
  raw <- sim_cate_run(
    designs     = "dgp1_linear",
    n_vals      = 100L,
    kx_vals     = 2L,
    n_reps      = 4L,
    master_seed = 7L,
    verbose     = FALSE
  )

  tab <- summarise_sim_results(raw)
  expect_s3_class(tab, "data.frame")
  # 1 method x 1 design x 1 n x 1 kx = 1 row
  expect_equal(nrow(tab), 1L)
  # tau_hat remains raw; the summary only records its empirical dispersion.
  expect_true(all(c("method", "design", "n", "kx", "n_reps", "n_failed",
                    "bias", "rmse", "se_mean", "cover", "ci_len",
                    "emp_sd", "se_to_sd", "neff")
                  %in% names(tab)))
  expect_false("tau_hat" %in% names(tab))
  expect_false("tau0"    %in% names(tab))
  # method column should come first
  expect_equal(names(tab)[1L], "method")
  expect_equal(tab$n_reps, 4L)
  expect_true(tab$rmse >= 0)
  expect_true(tab$cover >= 0 && tab$cover <= 1)
})

test_that("summarise_sim_results: multi-method returns one row per method", {
  raw <- sim_cate_run(
    designs     = "dgp1_linear",
    n_vals      = 100L,
    kx_vals     = 2L,
    n_reps      = 3L,
    master_seed = 7L,
    methods     = c("dnn_oracle", "ols_interact"),
    verbose     = FALSE
  )

  tab <- summarise_sim_results(raw)
  expect_equal(nrow(tab), 2L)
  expect_true(all(c("dnn_oracle", "ols_interact") %in% tab$method))
})

test_that("summarise_sim_results separates and retains GRF nuisance routes", {
  raw <- data.frame(
    method = "causal_forest",
    design = "dgp1_linear",
    n = 100L,
    kx = 2L,
    candidate_id = NA_character_,
    lambda = NA_real_,
    s = NA_integer_,
    bias = c(0.1, 0.2, -0.1, -0.2),
    sq_err = c(0.01, 0.04, 0.01, 0.04),
    se = c(0.2, 0.2, 0.3, 0.3),
    cover = c(1, 1, 0, 1),
    ci_len = c(0.8, 0.8, 1.2, 1.2),
    status = "ok",
    method_variant = rep(
      c("grf_oracle_nuisance", "grf_internal_nuisance"), each = 2L
    ),
    forest_nuisance = rep(c("oracle", "internal"), each = 2L),
    y_hat_supplied = rep(c(TRUE, FALSE), each = 2L),
    w_hat_supplied = rep(c(TRUE, FALSE), each = 2L),
    stringsAsFactors = FALSE
  )

  tab <- summarise_sim_results(raw)
  oracle <- tab[tab$forest_nuisance == "oracle", , drop = FALSE]
  internal <- tab[tab$forest_nuisance == "internal", , drop = FALSE]

  expect_equal(nrow(tab), 2L)
  expect_equal(oracle$n_reps, 2L)
  expect_identical(oracle$method_variant, "grf_oracle_nuisance")
  expect_identical(oracle$forest_nuisance, "oracle")
  expect_true(oracle$y_hat_supplied)
  expect_true(oracle$w_hat_supplied)
  expect_equal(internal$n_reps, 2L)
  expect_identical(internal$method_variant, "grf_internal_nuisance")
  expect_identical(internal$forest_nuisance, "internal")
  expect_false(internal$y_hat_supplied)
  expect_false(internal$w_hat_supplied)
})

test_that("summarise_sim_results preserves non-GRF and legacy raw behavior", {
  raw <- data.frame(
    method = rep("dnn_oracle", 2L),
    design = "dgp1_linear",
    n = 100L,
    kx = 2L,
    candidate_id = NA_character_,
    lambda = NA_real_,
    s = 12L,
    bias = c(0.1, -0.1),
    sq_err = c(0.01, 0.01),
    se = c(0.2, 0.2),
    cover = c(1, 1),
    ci_len = c(0.8, 0.8),
    status = "ok",
    method_variant = NA_character_,
    forest_nuisance = NA_character_,
    y_hat_supplied = NA,
    w_hat_supplied = NA,
    stringsAsFactors = FALSE
  )

  tab <- summarise_sim_results(raw)

  expect_equal(nrow(tab), 1L)
  expect_identical(tab$method, "dnn_oracle")
  expect_equal(tab$n_reps, 2L)
  expect_true(is.na(tab$method_variant))
  expect_true(is.na(tab$forest_nuisance))
  expect_true(is.na(tab$y_hat_supplied))
  expect_true(is.na(tab$w_hat_supplied))
  expect_type(tab$method_variant, "character")
  expect_type(tab$forest_nuisance, "character")
  expect_type(tab$y_hat_supplied, "logical")
  expect_type(tab$w_hat_supplied, "logical")

  route_cols <- c(
    "method_variant", "forest_nuisance", "y_hat_supplied", "w_hat_supplied"
  )
  legacy <- raw[, setdiff(names(raw), route_cols), drop = FALSE]
  legacy$method <- "causal_forest"
  legacy_tab <- NULL
  expect_no_warning(legacy_tab <- summarise_sim_results(legacy))

  expect_equal(nrow(legacy_tab), 1L)
  expect_identical(legacy_tab$method, "causal_forest")
  expect_true(all(is.na(legacy_tab[, route_cols, drop = FALSE])))
  expect_type(legacy_tab$method_variant, "character")
  expect_type(legacy_tab$forest_nuisance, "character")
  expect_type(legacy_tab$y_hat_supplied, "logical")
  expect_type(legacy_tab$w_hat_supplied, "logical")
})

test_that("summarise_sim_results: error on missing columns", {
  bad <- data.frame(design = "x", n = 100, kx = 2, bias = 0)
  expect_error(summarise_sim_results(bad), "missing required columns")
})

test_that("summarise_sim_results reports branch, radius, and studentized shape", {
  raw <- sim_cate_run(
    designs = "dgp2_nonlinear",
    n_vals = 160L,
    kx_vals = 2L,
    n_reps = 12L,
    master_seed = 707L,
    k = 1L,
    s_fn = function(n, kx) 24L,
    methods = c("dnn_fullsample", "dnn_oracle"),
    verbose = FALSE
  )
  tab <- summarise_sim_results(raw)
  required <- c(
    "bias_to_se", "localized_truth_mean", "branch_drift",
    "estimator_remainder", "standardized_branch_drift",
    "physical_radius_90_mean", "physical_radius_90_p25",
    "physical_radius_90_p50", "physical_radius_90_p75",
    "studentized_mean", "studentized_sd", "studentized_skewness",
    "studentized_excess_kurtosis", "studentized_q025", "studentized_q975"
  )
  expect_true(all(required %in% names(tab)))
  expect_true(all(is.finite(as.matrix(tab[, required]))))
  expect_equal(tab$bias, tab$branch_drift + tab$estimator_remainder, tolerance = 1e-12)
})

# ── summarise_sim_results: n_failed counting ──────────────────────────────────

test_that("summarise_sim_results: n_failed counts NA-metric rows correctly", {
  # Build a raw data.frame with 2 good reps and 1 failed (NA) rep
  raw_good <- sim_cate_run(
    designs     = "dgp1_linear",
    n_vals      = 100L,
    kx_vals     = 2L,
    n_reps      = 2L,
    master_seed = 1L,
    verbose     = FALSE
  )
  # Inject one NA row (simulates a rep-level failure)
  na_row        <- raw_good[1L, ]
  na_row$bias   <- NA_real_
  na_row$sq_err <- NA_real_
  na_row$se     <- NA_real_
  na_row$cover  <- NA_real_
  na_row$ci_len <- NA_real_
  raw_with_fail <- rbind(raw_good, na_row)

  tab <- summarise_sim_results(raw_with_fail)
  expect_equal(nrow(tab), 1L)
  expect_equal(tab$n_reps,   3L)
  expect_equal(tab$n_failed, 1L)
  # Metrics should be computed from the 2 good rows only
  expect_true(is.finite(tab$bias))
  expect_true(is.finite(tab$rmse))
  expect_true(is.finite(tab$emp_sd))
  expect_true(is.finite(tab$se_to_sd))
})

test_that("summarise_sim_results: n_failed = 0 when no NAs present", {
  raw <- sim_cate_run(
    designs     = "dgp2_nonlinear",
    n_vals      = 100L,
    kx_vals     = 2L,
    n_reps      = 3L,
    master_seed = 5L,
    verbose     = FALSE
  )
  tab <- summarise_sim_results(raw)
  expect_equal(tab$n_failed, 0L)
})
