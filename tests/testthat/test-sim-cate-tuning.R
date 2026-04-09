test_that("sim_cate_run tuning matches legacy s_fn for paper default", {
  base_args <- list(
    designs = "dgp1_linear",
    n_vals = 120L,
    kx_vals = 2L,
    n_reps = 2L,
    master_seed = 11L,
    verbose = FALSE
  )

  raw_fn <- suppressWarnings(do.call(sim_cate_run, c(base_args, list(
    s_fn = function(n, kx) floor(n^((kx + 1L) / (kx + 2L)))
  ))))
  raw_tun <- suppressWarnings(do.call(sim_cate_run, c(base_args, list(
    tuning = tuning_dnn_paper()
  ))))

  expect_equal(raw_tun$s, raw_fn$s)
  expect_equal(raw_tun$tau_hat, raw_fn$tau_hat, tolerance = 1e-12)
  expect_equal(raw_tun$se, raw_fn$se, tolerance = 1e-12)
})

test_that("sim_cate_run rejects simultaneous tuning and s_fn", {
  expect_error(
    sim_cate_run(
      designs = "dgp1_linear",
      n_vals = 80L,
      kx_vals = 2L,
      n_reps = 1L,
      master_seed = 1L,
      tuning = tuning_dnn_paper(),
      s_fn = function(n, kx) floor(n^((kx + 1L) / (kx + 2L))),
      verbose = FALSE
    ),
    "cannot provide both `tuning` and `s_fn`"
  )
})

test_that("sim_cate_run remains deterministic with tuning and fixed master_seed", {
  args <- list(
    designs = "dgp2_nonlinear",
    n_vals = 100L,
    kx_vals = 2L,
    n_reps = 3L,
    master_seed = 123L,
    tuning = tuning_dnn_paper(),
    verbose = FALSE
  )
  raw1 <- suppressWarnings(do.call(sim_cate_run, args))
  raw2 <- suppressWarnings(do.call(sim_cate_run, args))
  expect_equal(raw1, raw2)
})

test_that("sim_cate_run accepts tuning as registered name and bare function", {
  base_args <- list(
    designs = "dgp2_nonlinear",
    n_vals = 120L,
    kx_vals = 2L,
    n_reps = 2L,
    master_seed = 33L,
    k = 2L,
    verbose = FALSE
  )

  nm <- paste0("sim_fixed_", sample.int(1e6, 1L))
  register_tuning(nm, function(...) 4L)

  raw_name <- suppressWarnings(do.call(sim_cate_run, c(base_args, list(tuning = nm))))
  raw_fun <- suppressWarnings(do.call(sim_cate_run, c(base_args, list(tuning = function(...) 4L))))
  raw_fn <- suppressWarnings(do.call(sim_cate_run, c(base_args, list(s_fn = function(n, kx) 4L))))

  expect_equal(raw_name$s, raw_fn$s)
  expect_equal(raw_fun$s, raw_fn$s)
  expect_equal(raw_name$tau_hat, raw_fn$tau_hat, tolerance = 1e-12)
  expect_equal(raw_fun$tau_hat, raw_fn$tau_hat, tolerance = 1e-12)
})

test_that("sim_cate_run resolves tuning per (n, kx) cell", {
  tun_cell <- as_tuning(
    function(n, kx, ...) max(2L, floor((n + 2L * kx) / 40L)),
    name = "cellwise_rule"
  )

  raw <- suppressWarnings(sim_cate_run(
    designs = "dgp1_linear",
    n_vals = c(80L, 120L),
    kx_vals = c(2L, 3L),
    n_reps = 1L,
    master_seed = 44L,
    k = 2L,
    tuning = tun_cell,
    verbose = FALSE
  ))

  cells <- unique(raw[, c("n", "kx", "s")])
  expected_s <- mapply(
    function(n, kx) max(2L, floor((n + 2L * kx) / 40L)),
    cells$n, cells$kx
  )
  expect_equal(cells$s, as.integer(expected_s))
})

test_that("sim_cate_run default s_fn path is unchanged by tuning feature", {
  base <- list(
    designs = "dgp1_linear",
    n_vals = 100L,
    kx_vals = 2L,
    n_reps = 2L,
    master_seed = 55L,
    k = 2L,
    verbose = FALSE
  )

  raw_default <- suppressWarnings(do.call(sim_cate_run, base))
  raw_explicit <- suppressWarnings(do.call(sim_cate_run, c(base, list(
    s_fn = function(n, kx) floor(n^((kx + 1L) / (kx + 2L)))
  ))))

  expect_equal(raw_default$s, raw_explicit$s)
  expect_equal(raw_default$tau_hat, raw_explicit$tau_hat, tolerance = 1e-12)
  expect_equal(raw_default$se, raw_explicit$se, tolerance = 1e-12)
})

test_that("infeasible tuning values do not crash sim_cate_run and produce NA rows", {
  raw <- expect_warning(
    sim_cate_run(
      designs = "dgp1_linear",
      n_vals = 80L,
      kx_vals = 2L,
      n_reps = 1L,
      master_seed = 66L,
      k = 3L,
      tuning = as_tuning(function(...) 1L, name = "infeasible_small_s"),
      verbose = FALSE
    ),
    "method 'dnn_oracle' failed"
  )

  expect_true(all(raw$s == 1L))
  expect_true(all(is.na(raw$tau_hat)))
  expect_true(all(is.na(raw$se)))
  expect_true(all(is.na(raw$bias)))
  expect_true(all(is.na(raw$cover)))
})

test_that("tuning_dnn_paper kx_override controls cell-wise s variation", {
  base <- list(
    designs = "dgp2_nonlinear",
    n_vals = 100L,
    kx_vals = c(2L, 4L),
    n_reps = 1L,
    master_seed = 77L,
    k = 2L,
    verbose = FALSE
  )

  raw_no_override <- suppressWarnings(do.call(sim_cate_run, c(base, list(
    tuning = tuning_dnn_paper()
  ))))
  raw_override <- suppressWarnings(do.call(sim_cate_run, c(base, list(
    tuning = tuning_dnn_paper(kx_override = 2L)
  ))))

  expect_gt(length(unique(raw_no_override$s)), 1L)
  expect_equal(length(unique(raw_override$s)), 1L)
})
