.rng_seeded_calls <- function() {
  calls <- list(
    np = function(seed) np_demo_data(16, seed = seed)$data,
    cate = function(seed) cate_demo_data(16, seed = seed)$data,
    plqr = function(seed) plqr_demo_data(16, seed = seed)$data,
    clate = function(seed) clate_demo_data(16, seed = seed)$data,
    iv = function(seed) iv_demo_data(16, seed = seed)$data,
    np_subsets = function(seed) np_demo_subsets(16, K = 4, s = 3, seed = seed),
    cate_subsets = function(seed) cate_demo_subsets(16, K = 4, s = 3, seed = seed),
    plqr_subsets = function(seed) plqr_demo_subsets(16, K = 4, s = 3, seed = seed)
  )
  for (design in c("dgp1_linear", "dgp2_nonlinear", "dgp3_sparse_overlap",
    "dgp4_clustered", "dgp5_correlated", "dgp6_nuisance_efficiency")) {
    calls[[design]] <- local({
      selected <- design
      function(seed) sim_cate_data(selected, 16, seed = seed)$data
    })
  }
  for (design in c("dgp1_iv_justid_linear", "dgp2_iv_homog_overid", "dgp3_iv_heterog_overid",
    "dgp4_iv_weak_local", "dgp5_iv_correlated_instruments", "dgp6_iv_nonlinear_firststage")) {
    calls[[design]] <- local({
      selected <- design
      function(seed) sim_iv_data(selected, 16, seed = seed)$data
    })
  }
  calls
}

test_that("seeded generators preserve caller state and retain their draws", {
  .test_rng_scope({
    for (generate in .rng_seeded_calls()) {
      set.seed(19)
      expected <- generate(NULL)
      set.seed(271)
      before <- .Random.seed
      expect_identical(generate(19), expected)
      expect_identical(.Random.seed, before)
      rm(".Random.seed", envir = .GlobalEnv)
      expect_identical(generate(19), expected)
      expect_false(exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
      set.seed(19)
      before_null <- .Random.seed
      generate(NULL)
      expect_false(identical(.Random.seed, before_null))
    }
  })
})

test_that("seeded replicate and grid entry points restore caller state", {
  .test_rng_scope({
    calls <- list(
      function(seed) sim_cate_one_rep(40, 2, 1, 5, "dgp2_nonlinear", seed = seed),
      function(seed) sim_iv_one_rep(40, 2, 2, 1, 5, "dgp2_iv_homog_overid", seed = seed),
      function(seed) sim_cate_run("dgp2_nonlinear", n_vals = 40, n_reps = 2,
        master_seed = seed, verbose = FALSE),
      function(seed) sim_iv_run("dgp2_iv_homog_overid", n_vals = 40, q_vals = 2,
        n_reps = 2, master_seed = seed, verbose = FALSE)
    )
    for (run in calls) {
      set.seed(271)
      before <- .Random.seed
      first <- run(19)
      expect_identical(.Random.seed, before)
      expect_identical(run(19), first)
      rm(".Random.seed", envir = .GlobalEnv)
      expect_identical(run(19), first)
      expect_false(exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
      set.seed(271)
      run(NULL)
      expect_false(identical(.Random.seed, before))
    }
  })
})

test_that("seed scopes restore state after errors and RNG-kind changes", {
  .test_rng_scope({
    errors <- list(
      function() np_demo_data(12, seed = 19, inject = "bad"),
      function() cate_demo_data(12, seed = 19, inject = "bad"),
      function() plqr_demo_data(12, seed = 19, inject = "bad"),
      function() sim_cate_data("dgp4_clustered", 12, seed = 19, inner_sd = -1),
      function() sim_cate_one_rep(40, 2, 1, design = "dgp2_nonlinear", seed = 19,
        tuning = function(...) {
          runif(1)
          stop("callback failure")
        }),
      function() sim_cate_run("dgp2_nonlinear", n_vals = 40, n_reps = 1,
        s_fn = function(n, kx) stop("test failure"), verbose = FALSE),
      function() sim_iv_run("dgp2_iv_homog_overid", n_vals = 40, n_reps = 1,
        s_fn = function(n, kx) stop("test failure"), verbose = FALSE)
    )
    for (run in errors) {
      set.seed(271)
      before <- .Random.seed
      expect_error(run())
      expect_identical(.Random.seed, before)
      rm(".Random.seed", envir = .GlobalEnv)
      expect_error(run())
      expect_false(exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
    }
    set.seed(271, kind = "L'Ecuyer-CMRG")
    before <- .Random.seed
    kind <- RNGkind()
    expect_error(sim_cate_run("dgp2_nonlinear", n_vals = 40, n_reps = 1,
      s_fn = function(n, kx) {
        RNGkind("Mersenne-Twister")
        stop("test failure")
      }, verbose = FALSE), "test failure")
    expect_identical(RNGkind(), kind)
    expect_identical(.Random.seed, before)
    np_demo_data(12, seed = 19)
    expect_identical(RNGkind(), kind)
    expect_identical(.Random.seed, before)
  })
})

test_that("seed scopes include successful user tuning callbacks", {
  .test_rng_scope({
    for (had_seed in c(TRUE, FALSE)) {
      set.seed(271)
      before <- .Random.seed
      if (!had_seed) rm(".Random.seed", envir = .GlobalEnv)
      result <- sim_cate_one_rep(40, 2, 1, design = "dgp2_nonlinear", seed = 19,
        tuning = function(...) {
          runif(1)
          5L
        })
      expect_true(is.finite(result$tau_hat))
      if (had_seed) expect_identical(.Random.seed, before)
      else expect_false(exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
    }
  })
})

test_that("replicate seeds determine stochastic tuning independently of caller state", {
  .test_rng_scope({
    set.seed(19)
    selected_s <- sample(3:20, 1L)
    expected <- sim_cate_one_rep(40, 2, 1, s = selected_s,
      design = "dgp2_nonlinear", seed = 19)

    for (ambient_seed in c(1L, 2L, NA_integer_)) {
      set.seed(if (is.na(ambient_seed)) 3L else ambient_seed)
      before <- .Random.seed
      if (is.na(ambient_seed)) rm(".Random.seed", envir = .GlobalEnv)
      actual <- sim_cate_one_rep(40, 2, 1, design = "dgp2_nonlinear", seed = 19,
        tuning = function(...) sample(3:20, 1L))
      expect_identical(actual, expected)
      if (is.na(ambient_seed)) {
        expect_false(exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
      } else {
        expect_identical(.Random.seed, before)
      }
    }
  })
})

test_that("deterministic public and native calls preserve present and absent seeds", {
  old <- options(RdistNN.use_compiled = TRUE, RdistNN.num_threads = 2L)
  on.exit(options(old), add = TRUE)
  calls <- list(
    public_distance = function() distance_eval(distance_l2(), matrix(1:3), 0),
    native_distance = function() rustat_dist_l2(matrix(1:3), 0),
    subsets = function() rustat_subset_row_sums_g(matrix(c(3L, 1L, 3L), 1L), c(.5, -1, 2)),
    exact_table = function() rustat_exact_rank_table(3L, 1L, 2L),
    accumulation = function() rustat_compensated_sum(c(2^100, 1, -2^100)),
    plqr = function() {
      ctx <- rustat_plqr_ctx_create(c(.5, 1.5, 2.5), matrix(c(1, 2, 3)),
        rep(0, 3L), matrix(0, 3L, 1L), .5, .5, 1L)
      rustat_plqr_g_all(ctx)
    }
  )
  expected <- list(c(1, 2, 3), c(1, 2, 3), 4.5, c(2/3, 1/3, 0), 1, c(-.5, 1, 1.5))
  .test_rng_scope({
    for (i in seq_along(calls)) for (had_seed in c(TRUE, FALSE)) {
      set.seed(271)
      before <- .Random.seed
      if (!had_seed) rm(".Random.seed", envir = .GlobalEnv)
      value <- calls[[i]]()
      expect_equal(as.numeric(value), expected[[i]], tolerance = 1e-15, info = names(calls)[i])
      if (had_seed) expect_identical(.Random.seed, before)
      else expect_false(exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE), info = names(calls)[i])
    }
  })
})
