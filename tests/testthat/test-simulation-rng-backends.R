.simulation_execution_scope <- function(code, with_future = FALSE) {
  run <- function() {
    opts <- options()
    old_plan <- if (with_future) future::plan("list") else NULL
    on.exit({
      # Restoring the caller's plan cleans up only the backend created here.
      tryCatch({
        if (with_future) future::plan(old_plan, .cleanup = TRUE)
      }, finally = {
        added <- setdiff(names(options()), names(opts))
        if (length(added)) options(stats::setNames(rep(list(NULL), length(added)), added))
        options(opts)
      })
    })
    if (with_future) future::plan(future::sequential, .cleanup = FALSE)
    force(code)
  }
  # Restore RNG state after future and options cleanup, including on errors.
  .test_rng_scope(run())
}

.simulation_base_runner <- function(fun) {
  # Force the existing optional-dependency branch without replacing its map or
  # numerical code. The same tests run the real branch in minimal installations.
  environment(fun) <- list2env(list(requireNamespace = function(package, ...) {
    if (identical(package, "future.apply")) FALSE
    else base::requireNamespace(package, ...)
  }), parent = environment(fun))
  fun
}

.simulation_rng_runs <- function(fallback = FALSE, stochastic = FALSE) {
  cate <- if (fallback) .simulation_base_runner(sim_cate_run) else sim_cate_run
  iv <- if (fallback) .simulation_base_runner(sim_iv_run) else sim_iv_run
  args <- list(n_vals = c(32L, 36L), kx_vals = 2L, n_reps = 2L,
    master_seed = 19L, verbose = FALSE)
  if (stochastic) args$tuning <- function(...) sample(3:8, 1L)
  else args$s_fn <- function(n, kx) 5L
  list(
    cate = do.call(cate, c(args, list(designs = "dgp2_nonlinear", methods = "dnn_oracle"))),
    iv = do.call(iv, c(args, list(designs = "dgp2_iv_homog_overid", q_vals = 2L,
      methods = "dnn_iv_oracle")))
  )
}

test_that("simulation fallback repeats for a fixed RNG configuration", {
  has_future_apply <- requireNamespace("future.apply", quietly = TRUE)
  .simulation_execution_scope({
    previous <- NULL
    for (kind in c("Mersenne-Twister", "L'Ecuyer-CMRG")) {
      RNGkind(kind, normal.kind = "Inversion", sample.kind = "Rejection")
      configured_kind <- RNGkind()
      set.seed(271L)
      before <- .Random.seed
      first <- .simulation_rng_runs(fallback = TRUE, stochastic = TRUE)
      .expect_identical_fast(.Random.seed, before)
      expect_identical(RNGkind(), configured_kind)
      .expect_identical_fast(.simulation_rng_runs(fallback = TRUE, stochastic = TRUE), first)
      rm(".Random.seed", envir = .GlobalEnv)
      .expect_identical_fast(.simulation_rng_runs(fallback = TRUE, stochastic = TRUE), first)
      expect_false(exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
      expect_true(all(is.finite(first$cate$tau_hat)))
      expect_true(all(is.finite(first$iv$theta_hat)))
      if (!is.null(previous)) expect_false(identical(first, previous))
      previous <- first
      if (!has_future_apply) {
        .expect_identical_fast(.simulation_rng_runs(stochastic = TRUE), first)
      }
    }
  })
})

test_that("future mapping and fallback have distinct documented RNG contracts", {
  skip_if_not_installed("future")
  skip_if_not_installed("future.apply")
  .simulation_execution_scope({
    for (kind in c("Mersenne-Twister", "L'Ecuyer-CMRG")) {
      RNGkind(kind, normal.kind = "Inversion", sample.kind = "Rejection")
      configured_kind <- RNGkind()
      set.seed(271L)
      before <- .Random.seed
      fixed <- .simulation_rng_runs()
      .expect_identical_fast(.simulation_rng_runs(), fixed)
      .expect_identical_fast(.Random.seed, before)
      expect_identical(RNGkind(), configured_kind)
      fallback <- .simulation_rng_runs(fallback = TRUE)
      for (study in names(fixed)) {
        expect_identical(fixed[[study]]$seed, fallback[[study]]$seed)
        if (identical(kind, "L'Ecuyer-CMRG")) {
          .expect_identical_fast(fixed[[study]], fallback[[study]])
        } else {
          estimate <- if (identical(study, "cate")) "tau_hat" else "theta_hat"
          expect_false(identical(fixed[[study]][[estimate]], fallback[[study]][[estimate]]))
        }
      }
      stochastic <- .simulation_rng_runs(stochastic = TRUE)
      .expect_identical_fast(.simulation_rng_runs(stochastic = TRUE), stochastic)
      fallback_stochastic <- .simulation_rng_runs(fallback = TRUE, stochastic = TRUE)
      for (study in names(stochastic)) {
        expect_identical(stochastic[[study]]$s[1:2], fallback_stochastic[[study]]$s[1:2])
        expect_false(identical(stochastic[[study]]$s, fallback_stochastic[[study]]$s))
      }
    }
  }, with_future = TRUE)
})

test_that("oracle simulations agree across future sequential and multisession plans", {
  skip_if_not_installed("future")
  skip_if_not_installed("future.apply")
  installed <- file.exists(file.path(.libPaths(), "RdistNN", "Meta", "package.rds"))
  if (!any(installed)) {
    skip("Multisession workers need an installed RdistNN; source-only loading is unavailable to workers.")
  }
  .simulation_execution_scope({
    sequential <- list()
    for (kind in c("Mersenne-Twister", "L'Ecuyer-CMRG")) {
      RNGkind(kind, normal.kind = "Inversion", sample.kind = "Rejection")
      sequential[[kind]] <- lapply(c(FALSE, TRUE), function(stochastic) {
        .simulation_rng_runs(stochastic = stochastic)
      })
    }
    launch_error <- tryCatch({
      future::plan(future::multisession, workers = 2L)
      future::value(future::future(TRUE))
      NULL
    }, error = identity)
    if (!is.null(launch_error)) {
      skip(paste("Two-worker multisession is unavailable:", conditionMessage(launch_error)))
    }
    # Avoid silently comparing a source checkout with an unrelated installation.
    worker_signature <- future::value(future::future({
      lapply(c("sim_cate_one_rep", "sim_iv_one_rep"), function(name) {
        deparse(body(getExportedValue("RdistNN", name)))
      })
    }))
    local_signature <- lapply(c("sim_cate_one_rep", "sim_iv_one_rep"), function(name) {
      deparse(body(getExportedValue("RdistNN", name)))
    })
    expect_identical(worker_signature, local_signature)
    for (kind in names(sequential)) {
      RNGkind(kind, normal.kind = "Inversion", sample.kind = "Rejection")
      configured_kind <- RNGkind()
      set.seed(271L)
      before <- .Random.seed
      for (i in seq_along(sequential[[kind]])) {
        stochastic <- c(FALSE, TRUE)[[i]]
        actual <- .simulation_rng_runs(stochastic = stochastic)
        .expect_identical_fast(actual, sequential[[kind]][[i]])
        .expect_identical_fast(.simulation_rng_runs(stochastic = stochastic), actual)
        .expect_identical_fast(.Random.seed, before)
        expect_identical(RNGkind(), configured_kind)
      }
    }
  }, with_future = TRUE)
})

test_that("simulation test scopes restore execution state after errors", {
  .simulation_execution_scope({
    options(RdistNN.test_scope_existing = "before")
    for (with_future in unique(c(FALSE, requireNamespace("future", quietly = TRUE)))) {
      for (seeded in c(TRUE, FALSE)) {
        set.seed(271L)
        before <- .Random.seed
        kind <- RNGkind()
        old_plan <- if (with_future) future::plan("list") else NULL
        if (!seeded) rm(".Random.seed", envir = .GlobalEnv)
        expect_error(.simulation_execution_scope({
          options(RdistNN.test_scope_existing = "changed", RdistNN.test_scope_added = TRUE)
          RNGkind("L'Ecuyer-CMRG", normal.kind = "Inversion")
          set.seed(99L)
          stop("backend fixture failure")
        }, with_future = with_future), "backend fixture failure")
        expect_identical(getOption("RdistNN.test_scope_existing"), "before")
        expect_null(getOption("RdistNN.test_scope_added"))
        expect_identical(RNGkind(), kind)
        if (seeded) .expect_identical_fast(.Random.seed, before)
        else expect_false(exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
        if (with_future) .expect_identical_fast(future::plan("list"), old_plan)
      }
    }
  })
})
