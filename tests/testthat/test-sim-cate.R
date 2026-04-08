test_that("sim_cate_data: output schema is correct for each design", {
  designs <- c("dgp1_linear", "dgp2_nonlinear", "dgp3_sparse_overlap",
               "dgp4_clustered", "dgp5_correlated")
  expected_cols <- c(
    "X1", "X2", "D", "Y",
    "e_oracle", "m_oracle", "m0_oracle", "m1_oracle",
    "tau_oracle", "rho_oracle"
  )

  for (des in designs) {
    dgp <- sim_cate_data(des, n = 50L, kx = 2L, seed = 1L)

    expect_type(dgp, "list")
    expect_named(dgp, c("data", "x0", "tau_x0", "truth"))
    expect_s3_class(dgp$data, "data.frame")
    expect_equal(nrow(dgp$data), 50L)
    expect_true(all(expected_cols %in% names(dgp$data)),
                info = paste("missing cols in", des))

    expect_length(dgp$x0, 2L)
    expect_type(dgp$tau_x0, "double")
    expect_named(dgp$truth, c("tau_fn", "e_fn", "m0_fn", "m_fn"))
  }
})

test_that("sim_cate_data: deterministic with fixed seed", {
  dgp1 <- sim_cate_data("dgp2_nonlinear", n = 80L, kx = 2L, seed = 7L)
  dgp2 <- sim_cate_data("dgp2_nonlinear", n = 80L, kx = 2L, seed = 7L)
  expect_equal(dgp1$data, dgp2$data)
  expect_equal(dgp1$tau_x0, dgp2$tau_x0)
})

test_that("sim_cate_data: different seeds produce different data", {
  dgp1 <- sim_cate_data("dgp1_linear", n = 80L, kx = 2L, seed = 1L)
  dgp2 <- sim_cate_data("dgp1_linear", n = 80L, kx = 2L, seed = 2L)
  expect_false(identical(dgp1$data$Y, dgp2$data$Y))
})

test_that("sim_cate_data: propensity is clipped to (clip, 1-clip)", {
  dgp <- sim_cate_data("dgp3_sparse_overlap", n = 200L, kx = 2L, seed = 1L,
                       clip = 0.05)
  expect_true(all(dgp$data$e_oracle >= 0.05))
  expect_true(all(dgp$data$e_oracle <= 0.95))
})

test_that("sim_cate_data: m1_oracle = m0_oracle + tau_oracle", {
  dgp <- sim_cate_data("dgp2_nonlinear", n = 100L, kx = 2L, seed = 3L)
  expect_equal(dgp$data$m1_oracle,
               dgp$data$m0_oracle + dgp$data$tau_oracle,
               tolerance = 1e-12)
})

test_that("sim_cate_data: unknown design gives informative error", {
  expect_error(
    sim_cate_data("bad_design", n = 10L, kx = 2L),
    "Unknown design"
  )
})

test_that("sim_cate_data: kx < 2 is rejected", {
  expect_error(sim_cate_data("dgp1_linear", n = 10L, kx = 1L))
})

test_that("sim_cate_data: sigma validation is explicit", {
  expect_error(
    sim_cate_data("dgp1_linear", n = 10L, kx = 2L, sigma = -0.1),
    "sigma"
  )
  expect_error(
    sim_cate_data("dgp1_linear", n = 10L, kx = 2L, sigma = c(0.1, 0.2)),
    "sigma"
  )
  expect_error(
    sim_cate_data("dgp1_linear", n = 10L, kx = 2L, sigma = NA_real_),
    "sigma"
  )
})

test_that("sim_cate_data: clip validation is explicit", {
  expect_error(
    sim_cate_data("dgp1_linear", n = 10L, kx = 2L, clip = 0),
    "clip"
  )
  expect_error(
    sim_cate_data("dgp1_linear", n = 10L, kx = 2L, clip = 0.5),
    "clip"
  )
  expect_error(
    sim_cate_data("dgp1_linear", n = 10L, kx = 2L, clip = Inf),
    "clip"
  )
})

test_that("sim_cate_data: x0 validation checks length and finiteness", {
  expect_error(
    sim_cate_data("dgp1_linear", n = 10L, kx = 2L, x0 = c(0, 0, 0)),
    "x0 length must equal ncol\\(X\\)"
  )
  expect_error(
    sim_cate_data("dgp1_linear", n = 10L, kx = 2L, x0 = c(0, NA_real_)),
    "x0 must be finite"
  )
})

test_that("sim_cate_one_rep: returns a data.frame with correct structure (single method)", {
  r <- sim_cate_one_rep(
    n = 100L, kx = 2L, k = 3L,
    s = floor(100^(3/4)),
    design = "dgp1_linear", seed = 42L
  )
  expect_s3_class(r, "data.frame")
  expect_equal(nrow(r), 1L)
  expect_true(all(c("method", "tau_hat", "tau0", "bias", "sq_err",
                    "se", "cover", "ci_len", "neff") %in% names(r)))
  expect_equal(r$method, "dnn_oracle")
  expect_true(all(is.finite(unlist(r[, c("tau_hat", "tau0", "bias",
                                          "sq_err", "se", "cover",
                                          "ci_len", "neff")]))))
  expect_true(r$se > 0)
  expect_true(r$cover %in% c(0, 1))
  expect_true(r$neff > 0)
})

test_that("sim_cate_one_rep: returns one row per method when multiple methods given", {
  r <- sim_cate_one_rep(
    n = 100L, kx = 2L, k = 3L,
    s = floor(100^(3/4)),
    design = "dgp1_linear", seed = 42L,
    methods = c("dnn_oracle", "ols_interact")
  )
  expect_s3_class(r, "data.frame")
  expect_equal(nrow(r), 2L)
  expect_equal(r$method, c("dnn_oracle", "ols_interact"))
})

test_that("sim_cate_one_rep: deterministic with fixed seed", {
  args <- list(n = 80L, kx = 2L, k = 3L, s = 12L,
               design = "dgp2_nonlinear", seed = 99L)
  r1 <- do.call(sim_cate_one_rep, args)
  r2 <- do.call(sim_cate_one_rep, args)
  expect_equal(r1, r2)
})

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
  # tau_hat and tau0 are per-rep; they are NOT aggregated into the summary
  expect_true(all(c("method", "design", "n", "kx", "n_reps", "n_failed",
                    "bias", "rmse", "se_mean", "cover", "ci_len", "neff")
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

test_that("summarise_sim_results: error on missing columns", {
  bad <- data.frame(design = "x", n = 100, kx = 2, bias = 0)
  expect_error(summarise_sim_results(bad), "missing required columns")
})

# ── DGP 4: clustered covariate design ────────────────────────────────────────

test_that("sim_cate_data dgp4_clustered: propensity is clipped to (clip, 1-clip)", {
  dgp <- sim_cate_data("dgp4_clustered", n = 200L, kx = 2L, seed = 1L, clip = 0.05)
  expect_true(all(dgp$data$e_oracle >= 0.05))
  expect_true(all(dgp$data$e_oracle <= 0.95))
})

test_that("sim_cate_data dgp4_clustered: m1_oracle = m0_oracle + tau_oracle", {
  dgp <- sim_cate_data("dgp4_clustered", n = 100L, kx = 2L, seed = 3L)
  expect_equal(dgp$data$m1_oracle,
               dgp$data$m0_oracle + dgp$data$tau_oracle,
               tolerance = 1e-12)
})

test_that("sim_cate_data dgp4_clustered: mixture creates dense cluster near x0", {
  # With mix_prob = 0.7 and inner_sd = 0.4, the majority of obs should be near
  # x0 = (0, 0), giving a bimodal X1 distribution.  We verify the dense inner
  # cluster exists: at least 50% of obs should have |X1| < 1 (inner cluster has
  # sd=0.4 so nearly all of them land there, outer is at X1 ≈ 3 by default).
  dgp <- sim_cate_data("dgp4_clustered", n = 500L, kx = 2L, seed = 7L)
  prop_near <- mean(abs(dgp$data$X1) < 1)
  expect_gt(prop_near, 0.5)
})

test_that("sim_cate_data dgp4_clustered: outer cluster is well separated from x0", {
  # Outer cluster is centred at (outer_x1=3, 0).  With n large, some fraction of
  # obs should have X1 > 2.
  dgp <- sim_cate_data("dgp4_clustered", n = 500L, kx = 2L, seed = 8L)
  expect_gt(mean(dgp$data$X1 > 2), 0.05)
})

test_that("sim_cate_data dgp4_clustered: tau_x0 matches tau_fn at x0", {
  dgp <- sim_cate_data("dgp4_clustered", n = 50L, kx = 2L, seed = 2L)
  tau_at_x0 <- dgp$truth$tau_fn(matrix(dgp$x0, nrow = 1L))
  expect_equal(dgp$tau_x0, as.numeric(tau_at_x0), tolerance = 1e-12)
})

# ── DGP 5: correlated covariate design ───────────────────────────────────────

test_that("sim_cate_data dgp5_correlated: propensity is clipped to (clip, 1-clip)", {
  dgp <- sim_cate_data("dgp5_correlated", n = 200L, kx = 2L, seed = 1L, clip = 0.05)
  expect_true(all(dgp$data$e_oracle >= 0.05))
  expect_true(all(dgp$data$e_oracle <= 0.95))
})

test_that("sim_cate_data dgp5_correlated: m1_oracle = m0_oracle + tau_oracle", {
  dgp <- sim_cate_data("dgp5_correlated", n = 100L, kx = 2L, seed = 3L)
  expect_equal(dgp$data$m1_oracle,
               dgp$data$m0_oracle + dgp$data$tau_oracle,
               tolerance = 1e-12)
})

test_that("sim_cate_data dgp5_correlated: default x0 is (0.75, 0, ...)", {
  dgp <- sim_cate_data("dgp5_correlated", n = 50L, kx = 2L, seed = 1L)
  expect_equal(dgp$x0[1L], 0.75, tolerance = 1e-12)
  expect_equal(dgp$x0[2L], 0.0,  tolerance = 1e-12)
})

test_that("sim_cate_data dgp5_correlated: X1 and X2 are strongly correlated", {
  # Design uses rho = 0.85; empirical correlation should be close
  dgp <- sim_cate_data("dgp5_correlated", n = 2000L, kx = 2L, seed = 10L)
  r   <- cor(dgp$data$X1, dgp$data$X2)
  expect_gt(r, 0.75)
  expect_lt(r, 0.95)
})

test_that("sim_cate_data dgp5_correlated: tau_x0 matches tau_fn at x0", {
  dgp <- sim_cate_data("dgp5_correlated", n = 50L, kx = 2L, seed = 2L)
  tau_at_x0 <- dgp$truth$tau_fn(matrix(dgp$x0, nrow = 1L))
  expect_equal(dgp$tau_x0, as.numeric(tau_at_x0), tolerance = 1e-12)
})

# ── Benchmark methods via sim_cate_one_rep ────────────────────────────────────

test_that("sim_cate_one_rep: dnn_crossfit returns finite metrics", {
  r <- sim_cate_one_rep(
    n = 120L, kx = 2L, k = 3L,
    s = floor(120L^(3/4)),
    design = "dgp2_nonlinear",
    methods = "dnn_crossfit",
    seed = 42L
  )
  expect_equal(nrow(r), 1L)
  expect_equal(r$method, "dnn_crossfit")
  expect_true(is.finite(r$tau_hat))
  expect_true(r$se > 0)
  expect_true(r$cover %in% c(0, 1))
})

test_that("sim_cate_one_rep: local_constant returns finite metrics", {
  r <- sim_cate_one_rep(
    n = 120L, kx = 2L, k = 3L,
    s = floor(120L^(3/4)),
    design = "dgp2_nonlinear",
    methods = "local_constant",
    seed = 42L
  )
  expect_equal(nrow(r), 1L)
  expect_equal(r$method, "local_constant")
  expect_true(is.finite(r$tau_hat))
  expect_true(r$se > 0)
  expect_true(r$cover %in% c(0, 1))
})

test_that("sim_cate_one_rep: dnn_oracle_mahal returns finite metrics on dgp5_correlated", {
  r <- sim_cate_one_rep(
    n = 120L, kx = 2L, k = 3L,
    s = floor(120L^(3/4)),
    design = "dgp5_correlated",
    methods = "dnn_oracle_mahal",
    seed = 42L
  )
  expect_equal(nrow(r), 1L)
  expect_equal(r$method, "dnn_oracle_mahal")
  expect_true(is.finite(r$tau_hat))
  expect_true(r$se > 0)
  expect_true(r$cover %in% c(0, 1))
})

test_that("sim_cate_one_rep: all five methods run together without error", {
  r <- sim_cate_one_rep(
    n = 100L, kx = 2L, k = 3L,
    s = floor(100L^(3/4)),
    design = "dgp2_nonlinear",
    methods = c("dnn_oracle", "dnn_oracle_mahal", "dnn_crossfit",
                "local_constant", "ols_interact"),
    seed = 7L
  )
  expect_equal(nrow(r), 5L)
  expect_equal(r$method,
               c("dnn_oracle", "dnn_oracle_mahal", "dnn_crossfit",
                 "local_constant", "ols_interact"))
  expect_true(all(is.finite(r$tau_hat)))
  expect_true(all(r$se > 0))
})

# ── sim_cate_methods registry ─────────────────────────────────────────────────

test_that("sim_cate_methods: includes all five paper methods", {
  m <- sim_cate_methods()
  expect_true(all(c("dnn_oracle", "dnn_oracle_mahal", "dnn_crossfit",
                    "local_constant", "ols_interact") %in% m))
})

# ── sim_cate_run with all 5 designs ──────────────────────────────────────────

test_that("sim_cate_run: smoke test covers all 5 designs", {
  raw <- sim_cate_run(
    designs     = c("dgp1_linear", "dgp2_nonlinear", "dgp3_sparse_overlap",
                    "dgp4_clustered", "dgp5_correlated"),
    n_vals      = 100L,
    kx_vals     = 2L,
    n_reps      = 2L,
    master_seed = 1L,
    verbose     = FALSE
  )
  # 5 designs x 1 n x 1 kx x 2 reps x 1 method (default dnn_oracle) = 10
  expect_equal(nrow(raw), 10L)
  expect_true(all(c("dgp1_linear", "dgp2_nonlinear", "dgp3_sparse_overlap",
                    "dgp4_clustered", "dgp5_correlated") %in% raw$design))
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
