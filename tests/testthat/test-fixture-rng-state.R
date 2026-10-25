test_that("loading fixture helpers preserves the caller RNG state", {
  .test_rng_scope({
    for (file in c("helper-distances.R", "helper-jackknife.R", "helper-user-cate-dnn.R", "helper-user-plqr-dnn.R", "helper-user-np-dnn.R")) {
      set.seed(271L)
      before <- .Random.seed
      sys.source(test_path(file), envir = new.env(parent = environment()))
      expect_identical(.Random.seed, before, info = file)
      rm(".Random.seed", envir = .GlobalEnv)
      sys.source(test_path(file), envir = new.env(parent = environment()))
      expect_false(exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE), info = file)
    }
  })
})

test_that("shared random fixtures are deterministic without consuming caller draws", {
  .test_rng_scope({
    fixtures <- list(
      query = function() .make_X0(3L, 2L),
      distance_training = function() gen_X(),
      distance_targets = function() gen_X0(),
      distance_target = function() gen_x0(),
      np_interface = function() .make_np_iface_data(),
      tdnn_interface = function() .make_tdnn_iface_data(),
      np = function() .make_np_data(),
      tdnn = function() .make_tdnn_data(),
      cate = function() .make_data(),
      cate_frame = function() .make_df_data(),
      plqr = function() .make_plqr_data(),
      plqr_inference = function() .make_plqr_inference_data()
    )
    for (fixture in fixtures) {
      set.seed(271L)
      before <- .Random.seed
      first <- fixture()
      expect_identical(.Random.seed, before)
      set.seed(991L)
      expect_identical(fixture(), first)
      rm(".Random.seed", envir = .GlobalEnv)
      expect_identical(fixture(), first)
      expect_false(exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
    }
  })
})

test_that("test RNG scopes restore kinds and seed presence after errors", {
  .test_rng_scope({
    for (seeded in c(TRUE, FALSE)) {
      set.seed(271L)
      before <- .Random.seed
      kind <- RNGkind()
      if (!seeded) rm(".Random.seed", envir = .GlobalEnv)
      expect_error(.test_rng_scope({
        RNGkind("L'Ecuyer-CMRG", normal.kind = "Inversion")
        set.seed(99L)
        stop("fixture failure")
      }), "fixture failure")
      expect_identical(RNGkind(), kind)
      if (seeded) expect_identical(.Random.seed, before)
      else expect_false(exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
    }
  })
})

test_that("test-file execution restores caller RNG after inline seed setup", {
  .test_rng_scope({
    cases <- list(
      list(file = "test-api-distances.R",
        desc = "as_distance() wraps a function and is_distance() recognizes it"),
      list(file = "test-reporting-recovery.R",
        desc = "grid callbacks do not change options or RNG while preserving warning handling")
    )
    for (case in cases) for (seeded in c(TRUE, FALSE)) {
      set.seed(271L)
      if (!seeded) rm(".Random.seed", envir = .GlobalEnv)
      before <- get0(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
      kind <- RNGkind()
      # Run the real testthat lifecycle, including its deferred block cleanup.
      result <- testthat::test_file(test_path(case$file), desc = case$desc,
        reporter = "silent", load_helpers = FALSE)
      after <- get0(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
      .expect_identical_fast(after, before, info = case$file)
      expect_identical(RNGkind(), kind, info = case$file)
      status <- as.data.frame(result)
      expect_equal(sum(status$failed), 0L, info = case$file)
      expect_false(any(status$error), info = case$file)
    }
  })
})

test_that("suite-local fixtures load and repeat without changing caller RNG state", {
  .test_rng_scope({
    cases <- list(
      list("test-api-ranks.R", ".rand_ranks_mat", list(n = 8L, m = 3L)),
      list("test-candidate-fit.R", ".candidate_view_np_case"),
      list("test-plot-builders-contracts.R", ".pbc_np_fit"),
      list("test-plot-builders-contracts.R", ".pbc_cate_targets", list(m = 3L)),
      list("test-plot-builders-contracts.R", ".pbc_cate_fit"),
      list("test-plot-builders.R", ".pb_np_profile_fit"),
      list("test-plot-builders.R", ".pb_cate_fit"),
      list("test-plot-data-contracts.R", ".pdc_np_fit"),
      list("test-plot-data.R", ".pd_np_fit"),
      list("test-plot-diagnostics-contracts.R", ".plot_diag_contract_cate_df"),
      list("test-plot-diagnostics.R", ".plot_diag_cate_df"),
      list("test-plot-diagnostics.R", ".plot_diag_matrix_cate_fit"),
      list("test-plot-specialized.R", ".plot_special_cate_fit"),
      list("test-rank-weight-table-contract.R", ".rank_table_rank_matrix", list(n = 8L, m = 3L)),
      list("test-rank-weight-tdnn-asymp.R", ".make_rank_matrix", list(n = 8L, m = 3L)),
      list("test-score-cate.R", ".get_cate_fx"),
      list("test-score-npr.R", "reg_demo"),
      list("test-score-plqr.R", ".get_plqr_fx"),
      list("test-tdnn-jackknife.R", ".make_X_jk", list(n = 8L)),
      list("test-tuning-grid-diagnostics-plots.R", ".grid_plot_np_data"),
      list("test-user-cate-methods.R", ".make_cate_fit"),
      list("test-user-np-dnn-tuning-grid.R", ".np_grid_np_data")
    )
    # Load the real suite definitions while leaving its tests for their own run.
    # This checks behavior without matching source text or depending on file order.
    load_suite <- function(file) {
      env <- new.env(parent = environment())
      env$test_that <- function(...) invisible(NULL)
      sys.source(test_path(file), envir = env)
      env
    }
    snapshot <- function(value) {
      fit <- if (is.list(value) && !is.null(value$fit)) value$fit else value
      if (inherits(fit, c("rdistnn_np_result", "rdistnn_cate_result"))) {
        # Compare scientific inputs and outputs, not callback environments.
        return(list(X = fit$.model$x_mat, x0 = fit$x0,
          signal = if (is.null(fit$.model$y_vec)) fit$.model$psi else fit$.model$y_vec,
          estimate = if (is.null(fit$mu)) fit$tau else fit$mu, se = fit$se))
      }
      value
    }
    for (file in unique(c(vapply(cases, `[[`, "", 1L),
                          "test-user-cate-dnn-core.R", "test-user-plqr-dnn.R"))) {
      set.seed(271L)
      before <- .Random.seed
      load_suite(file)
      .expect_identical_fast(.Random.seed, before, info = file)
      rm(".Random.seed", envir = .GlobalEnv)
      load_suite(file)
      expect_false(exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE), info = file)
    }
    for (case in cases) {
      env <- load_suite(case[[1L]])
      fixture <- get(case[[2L]], envir = env, inherits = FALSE)
      args <- if (length(case) >= 3L) case[[3L]] else list()
      run <- function() snapshot(do.call(fixture, args))
      set.seed(271L)
      before <- .Random.seed
      first <- run()
      .expect_identical_fast(.Random.seed, before, info = case[[2L]])
      set.seed(991L)
      .expect_identical_fast(run(), first, info = case[[2L]])
      rm(".Random.seed", envir = .GlobalEnv)
      .expect_identical_fast(run(), first, info = case[[2L]])
      expect_false(exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE), info = case[[2L]])
    }
  })
})
