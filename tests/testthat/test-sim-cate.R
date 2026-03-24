test_that("sim_cate_data: output schema is correct for each design", {
  designs <- c("dgp1_linear", "dgp2_nonlinear", "dgp3_sparse_overlap")
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
