# Method registry, repeated runs, reproducibility, and failure metadata.
test_that("sim_cate_methods: returns character vector of known methods", {
  m <- sim_cate_methods()
  expect_type(m, "character")
  expect_true(length(m) >= 4L)
  expect_true(all(c("dnn_oracle", "dnn_crossfit", "local_constant",
                    "ols_interact") %in% m))
})

test_that("sim_cate_run: smoke test with 3 reps over a small grid", {
  raw <- sim_cate_run(
    designs     = c("dgp1_linear", "dgp2_nonlinear"),
    n_vals      = 100L,
    kx_vals     = 2L,
    n_reps      = 3L,
    master_seed = 1L,
    verbose     = FALSE
  )

  expect_s3_class(raw, "data.frame")
  # 2 designs x 1 n x 1 kx x 3 reps x 1 method (dnn_oracle) = 6
  expect_equal(nrow(raw), 6L)
  expect_true(all(c("method", "design", "n", "kx", "s", "seed",
                    "tau_hat", "tau0", "bias", "sq_err", "se",
                    "cover", "ci_len", "neff") %in% names(raw)))
  expect_true(all(is.finite(raw$se)))
  expect_true(all(raw$se > 0))
})

test_that("sim_cate_run: replicate seeds are unique and reproducible", {
  args <- list(
    designs = "dgp1_linear",
    n_vals = 100L,
    kx_vals = 2L,
    n_reps = 6L,
    master_seed = 123L,
    verbose = FALSE
  )
  raw1 <- do.call(sim_cate_run, args)
  raw2 <- do.call(sim_cate_run, args)

  expect_equal(raw1$seed, raw2$seed)
  expect_equal(length(unique(raw1$seed)), 6L)
})

test_that("sim_cate_run: no future.apply RNG misuse warning", {
  skip_if_not_installed("future.apply")

  warns <- character()
  invisible(withCallingHandlers(
    sim_cate_run(
      designs = "dgp1_linear",
      n_vals = 80L,
      kx_vals = 2L,
      n_reps = 2L,
      master_seed = 321L,
      verbose = FALSE
    ),
    warning = function(w) {
      warns <<- c(warns, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  ))

  expect_false(any(grepl("UNRELIABLE VALUE|future\\.seed", warns, ignore.case = TRUE)))
})

test_that("grf-missing benchmark warnings are emitted without call metadata", {
  testthat::local_mocked_bindings(
    .sim_cate_grf_available = function() FALSE,
    .package = "RdistNN"
  )

  m_cf <- getFromNamespace(".method_causal_forest", "RdistNN")
  m_llf <- getFromNamespace(".method_llf", "RdistNN")
  metric_names <- getFromNamespace(".sim_cate_metric_names", "RdistNN")()
  dgp <- sim_cate_data("dgp1_linear", n = 40L, kx = 2L, seed = 11L)
  x_cols <- grep("^X", names(dgp$data), value = TRUE)

  w_cf <- NULL
  cf_result <- withCallingHandlers(
    m_cf(dgp, x_cols = x_cols, alpha = 0.05),
    warning = function(w) {
      w_cf <<- w
      invokeRestart("muffleWarning")
    }
  )
  expect_s3_class(w_cf, "warning")
  expect_match(conditionMessage(w_cf), "grf not available; causal_forest returns NA\\.")
  expect_null(conditionCall(w_cf))
  expect_identical(cf_result$status, "failed")
  expect_true(all(is.na(cf_result[, metric_names, drop = FALSE])))
  expect_identical(cf_result$method_variant, "grf_oracle_nuisance")
  expect_identical(cf_result$forest_nuisance, "oracle")
  expect_true(cf_result$y_hat_supplied)
  expect_true(cf_result$w_hat_supplied)

  w_llf <- NULL
  llf_result <- withCallingHandlers(
    m_llf(dgp, x_cols = x_cols, alpha = 0.05),
    warning = function(w) {
      w_llf <<- w
      invokeRestart("muffleWarning")
    }
  )
  expect_s3_class(w_llf, "warning")
  expect_match(conditionMessage(w_llf), "grf not available; llf returns NA\\.")
  expect_null(conditionCall(w_llf))
  expect_identical(llf_result$status, "failed")
  expect_true(all(is.na(llf_result[, metric_names, drop = FALSE])))
  expect_identical(llf_result$method_variant, "grf_local_linear_oracle_nuisance")
  expect_identical(llf_result$forest_nuisance, "oracle")
  expect_true(llf_result$y_hat_supplied)
  expect_true(llf_result$w_hat_supplied)
})

test_that("method runtime failures retain registered GRF identity", {
  original_get <- getFromNamespace(".get_sim_method", "RdistNN")
  testthat::local_mocked_bindings(
    .get_sim_method = function(name) {
      if (identical(name, "causal_forest")) {
        return(function(...) stop("injected GRF runtime failure", call. = FALSE))
      }
      original_get(name)
    },
    .package = "RdistNN"
  )

  failed <- NULL
  expect_warning(
    failed <- sim_cate_one_rep(
      n = 40L, kx = 2L, k = 1L, s = 6L,
      design = "dgp1_linear", methods = "causal_forest", seed = 122L
    ),
    "injected GRF runtime failure"
  )
  metric_names <- getFromNamespace(".sim_cate_metric_names", "RdistNN")()

  expect_identical(failed$status, "failed")
  expect_true(all(is.na(failed[, metric_names, drop = FALSE])))
  expect_identical(failed$method_variant, "grf_oracle_nuisance")
  expect_identical(failed$forest_nuisance, "oracle")
  expect_true(failed$y_hat_supplied)
  expect_true(failed$w_hat_supplied)
})

test_that("outer replicate failures retain registered GRF identity", {
  failed <- NULL
  expect_warning(
    failed <- sim_cate_run(
      designs = "bad_design", n_vals = 40L, kx_vals = 2L,
      n_reps = 1L, master_seed = 123L, methods = "causal_forest",
      verbose = FALSE
    ),
    "rep seed"
  )
  metric_names <- getFromNamespace(".sim_cate_metric_names", "RdistNN")()

  expect_identical(failed$status, "failed")
  expect_true(all(is.na(failed[, metric_names, drop = FALSE])))
  expect_identical(failed$method_variant, "grf_oracle_nuisance")
  expect_identical(failed$forest_nuisance, "oracle")
  expect_true(failed$y_hat_supplied)
  expect_true(failed$w_hat_supplied)
})

test_that("unknown simulation methods are rejected before failure shaping", {
  resolver <- getFromNamespace(".sim_cate_method_metadata", "RdistNN")
  expect_error(resolver("not_registered"), "Unknown simulation method")
  expect_error(
    sim_cate_one_rep(
      n = 40L, kx = 2L, k = 1L, s = 6L,
      design = "dgp1_linear", methods = "not_registered", seed = 124L
    ),
    "Unknown simulation method"
  )
  expect_error(
    sim_cate_run(
      designs = "bad_design", n_vals = 40L, kx_vals = 2L,
      n_reps = 1L, methods = "not_registered", verbose = FALSE
    ),
    "Unknown simulation method"
  )
})

test_that("GRF comparator metadata identifies supplied oracle nuisances", {
  skip_if_not_installed("grf")
  method <- getFromNamespace(".method_causal_forest", "RdistNN")
  dgp <- sim_cate_data("dgp1_linear", n = 80L, kx = 2L, seed = 119L)
  x_cols <- grep("^X", names(dgp$data), value = TRUE)
  result <- method(dgp, x_cols = x_cols, alpha = 0.05)

  expect_identical(result$method_variant, "grf_oracle_nuisance")
  expect_identical(result$forest_nuisance, "oracle")
  expect_true(result$y_hat_supplied)
  expect_true(result$w_hat_supplied)
})

test_that("sim_cate_run: multi-method produces one row per method per rep", {
  raw <- sim_cate_run(
    designs     = "dgp1_linear",
    n_vals      = 100L,
    kx_vals     = 2L,
    n_reps      = 2L,
    master_seed = 1L,
    methods     = c("dnn_oracle", "ols_interact"),
    verbose     = FALSE
  )

  expect_s3_class(raw, "data.frame")
  # 1 design x 1 n x 1 kx x 2 reps x 2 methods = 4
  expect_equal(nrow(raw), 4L)
  expect_true(all(c("dnn_oracle", "ols_interact") %in% raw$method))
})

test_that("sim_cate_run raw schema includes typed GRF nuisance-route fields", {
  schema <- getFromNamespace(".sim_cate_raw_schema", "RdistNN")()
  route_fields <- c(
    "method_variant", "forest_nuisance", "y_hat_supplied", "w_hat_supplied"
  )

  expect_true(all(route_fields %in% schema))
  expect_lt(max(match(route_fields, schema)), match("status", schema))
})

test_that("sim_cate_run preserves GRF oracle nuisance metadata", {
  skip_if_not_installed("grf")
  raw <- sim_cate_run(
    designs = "dgp1_linear",
    n_vals = 80L,
    kx_vals = 2L,
    n_reps = 1L,
    master_seed = 120L,
    methods = "causal_forest",
    verbose = FALSE
  )
  successful <- raw[raw$status == "ok", , drop = FALSE]

  expect_gt(nrow(successful), 0L)
  expect_true(all(successful$method_variant == "grf_oracle_nuisance"))
  expect_true(all(successful$forest_nuisance == "oracle"))
  expect_true(all(successful$y_hat_supplied))
  expect_true(all(successful$w_hat_supplied))
  expect_type(raw$method_variant, "character")
  expect_type(raw$forest_nuisance, "character")
  expect_type(raw$y_hat_supplied, "logical")
  expect_type(raw$w_hat_supplied, "logical")
})

test_that("sim_cate_run keeps GRF metadata and typed missing DNN metadata", {
  skip_if_not_installed("grf")
  raw <- sim_cate_run(
    designs = "dgp1_linear",
    n_vals = 80L,
    kx_vals = 2L,
    n_reps = 1L,
    master_seed = 121L,
    methods = c("dnn_oracle", "causal_forest"),
    verbose = FALSE
  )
  grf_row <- raw[raw$method == "causal_forest" & raw$status == "ok", , drop = FALSE]
  dnn_row <- raw[raw$method == "dnn_oracle" & raw$status == "ok", , drop = FALSE]

  expect_gt(nrow(grf_row), 0L)
  expect_true(all(grf_row$method_variant == "grf_oracle_nuisance"))
  expect_true(all(grf_row$forest_nuisance == "oracle"))
  expect_true(all(grf_row$y_hat_supplied))
  expect_true(all(grf_row$w_hat_supplied))
  expect_gt(nrow(dnn_row), 0L)
  expect_true(all(is.na(dnn_row$method_variant)))
  expect_true(all(is.na(dnn_row$forest_nuisance)))
  expect_true(all(is.na(dnn_row$y_hat_supplied)))
  expect_true(all(is.na(dnn_row$w_hat_supplied)))
  expect_type(dnn_row$method_variant, "character")
  expect_type(dnn_row$forest_nuisance, "character")
  expect_type(dnn_row$y_hat_supplied, "logical")
  expect_type(dnn_row$w_hat_supplied, "logical")
})

# ── sim_cate_methods registry ─────────────────────────────────────────────────

test_that("sim_cate_methods: includes all six core paper methods", {
  m <- sim_cate_methods()
  expect_true(all(c("dnn_oracle", "dnn_oracle_mahal", "dnn_fullsample", "dnn_crossfit",
                    "local_constant", "ols_interact") %in% m))
})

# ── sim_cate_run with all 5 designs ──────────────────────────────────────────

test_that("sim_cate_run: smoke test covers all 6 designs", {
  raw <- sim_cate_run(
    designs     = c("dgp1_linear", "dgp2_nonlinear", "dgp3_sparse_overlap",
                    "dgp4_clustered", "dgp5_correlated", "dgp6_nuisance_efficiency"),
    n_vals      = 100L,
    kx_vals     = 2L,
    n_reps      = 2L,
    master_seed = 1L,
    verbose     = FALSE
  )
  # 6 designs x 1 n x 1 kx x 2 reps x 1 method (default dnn_oracle) = 12
  expect_equal(nrow(raw), 12L)
  expect_true(all(c("dgp1_linear", "dgp2_nonlinear", "dgp3_sparse_overlap",
                    "dgp4_clustered", "dgp5_correlated",
                    "dgp6_nuisance_efficiency") %in% raw$design))
  expect_true(all(is.finite(raw$se)))
})

test_that("sim_cate_run: dgp4 and dgp5 produce finite results with dnn_oracle_mahal", {
  raw <- sim_cate_run(
    designs     = c("dgp4_clustered", "dgp5_correlated"),
    n_vals      = 100L,
    kx_vals     = 2L,
    n_reps      = 2L,
    master_seed = 2L,
    methods     = c("dnn_oracle", "dnn_oracle_mahal"),
    verbose     = FALSE
  )
  # 2 designs x 2 reps x 2 methods = 8
  expect_equal(nrow(raw), 8L)
  expect_true(all(is.finite(raw$tau_hat)))
  expect_true(all(raw$se > 0))
})
