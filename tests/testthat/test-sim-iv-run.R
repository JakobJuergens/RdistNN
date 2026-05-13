test_that("sim_iv_methods lists DNN and classical IV methods", {
  methods <- sim_iv_methods()
  expect_true(all(c(
    "dnn_iv_oracle", "dnn_iv_fullsample", "dnn_iv_crossfit",
    "iv_oracle", "2sls_global", "2sls_interact", "local_2sls"
  ) %in% methods))
})

test_that("sim_iv_one_rep returns full IV metric schema", {
  res <- sim_iv_one_rep(
    n = 100L,
    kx = 2L,
    q = 2L,
    k = 2L,
    s = 9L,
    design = "dgp2_iv_homog_overid",
    methods = c("dnn_iv_oracle", "iv_oracle"),
    seed = 10L
  )

  expect_equal(nrow(res), 2L)
  expect_true(all(c(
    "method", "theta_hat", "theta0", "bias", "sq_err", "se", "cover",
    "ci_len", "neff", "J", "df", "pval_J", "reject_J"
  ) %in% names(res)))
  expect_true(all(is.finite(res$theta_hat)))
  expect_true(all(res$df == 1))
})

test_that("sim_iv_run is deterministic and skips duplicate just-identified q cells", {
  run_once <- function() {
    sim_iv_run(
      designs = "dgp1_iv_justid_linear",
      n_vals = 90L,
      kx_vals = 2L,
      q_vals = c(1L, 2L, 3L),
      n_reps = 2L,
      master_seed = 11L,
      k = 2L,
      s_fn = function(n, kx) 8L,
      methods = "iv_oracle",
      verbose = FALSE
    )
  }

  raw1 <- run_once()
  raw2 <- run_once()

  expect_equal(raw1, raw2)
  expect_equal(unique(raw1$q), 1L)
  expect_equal(nrow(raw1), 2L)
})

test_that("sim_iv_run expands IV cells, methods, and aligned seeds", {
  raw <- sim_iv_run(
    designs = c("dgp1_iv_justid_linear", "dgp2_iv_homog_overid"),
    n_vals = c(90L, 100L),
    kx_vals = 2L,
    q_vals = c(2L, 3L),
    n_reps = 2L,
    master_seed = 13L,
    k = 2L,
    s_fn = function(n, kx) 8L,
    methods = c("iv_oracle", "2sls_global"),
    verbose = FALSE
  )

  expect_equal(nrow(raw), 24L)
  expect_equal(sort(unique(raw$q[raw$design == "dgp1_iv_justid_linear"])), 1L)
  expect_equal(sort(unique(raw$q[raw$design == "dgp2_iv_homog_overid"])), c(2L, 3L))

  key <- paste(raw$design, raw$n, raw$kx, raw$q, raw$seed)
  method_count <- tapply(raw$method, key, function(x) length(unique(x)))
  expect_true(all(method_count == 2L))
})

test_that("sim_iv_one_rep reports NA rows when a method fails internally", {
  res <- sim_iv_one_rep(
    n = 8L,
    kx = 2L,
    q = 2L,
    k = 1L,
    s = 3L,
    design = "dgp2_iv_homog_overid",
    methods = "dnn_iv_crossfit",
    seed = 14L,
    n_folds = 20L
  )

  expect_equal(nrow(res), 1L)
  expect_equal(res$method, "dnn_iv_crossfit")
  expect_true(is.na(res$theta_hat))
  expect_true(is.na(res$se))
  expect_true(is.na(res$cover))
})

test_that("sim_iv_run validates grid controls", {
  expect_error(sim_iv_run(n_reps = 0L, verbose = FALSE), "n_reps")
  expect_error(sim_iv_run(q_vals = 4L, n_reps = 1L, verbose = FALSE), "q_vals")
  expect_error(
    sim_iv_run(
      designs = "dgp2_iv_homog_overid",
      n_vals = 80L,
      q_vals = 2L,
      n_reps = 1L,
      tuning = tuning_dnn_paper(),
      s_fn = function(n, kx) 8L,
      verbose = FALSE
    ),
    "tuning/s_fn"
  )
  expect_error(
    sim_iv_one_rep(
      n = 80L, kx = 2L, q = 2L, k = 1L, s = 8L,
      design = "dgp2_iv_homog_overid", methods = "bad_method"
    ),
    "Unknown IV simulation method"
  )
})

test_that("summarise_iv_results handles just-identified J columns", {
  raw <- sim_iv_run(
    designs = "dgp1_iv_justid_linear",
    n_vals = 90L,
    kx_vals = 2L,
    q_vals = c(1L, 2L),
    n_reps = 2L,
    master_seed = 12L,
    k = 2L,
    s_fn = function(n, kx) 8L,
    methods = "iv_oracle",
    verbose = FALSE
  )

  out <- summarise_iv_results(raw)
  expect_equal(nrow(out), 1L)
  expect_true(is.na(out$J_median))
  expect_true(is.na(out$pval_J_median))
  expect_true(is.na(out$reject_J))
})

test_that("summarise_iv_results counts failures and validates schema", {
  raw <- sim_iv_run(
    designs = "dgp2_iv_homog_overid",
    n_vals = 90L,
    kx_vals = 2L,
    q_vals = 2L,
    n_reps = 3L,
    master_seed = 15L,
    k = 2L,
    s_fn = function(n, kx) 8L,
    methods = "iv_oracle",
    verbose = FALSE
  )
  raw$bias[1L] <- NA_real_
  raw$sq_err[1L] <- NA_real_
  raw$se[1L] <- NA_real_
  raw$cover[1L] <- NA_real_
  raw$ci_len[1L] <- NA_real_

  out <- summarise_iv_results(raw)
  expect_equal(out$n_reps, 3L)
  expect_equal(out$n_failed, 1L)
  expect_true(is.finite(out$bias))
  expect_true(is.finite(out$rmse))

  expect_error(summarise_iv_results(raw[, setdiff(names(raw), "bias")]),
               "missing required columns")
})
