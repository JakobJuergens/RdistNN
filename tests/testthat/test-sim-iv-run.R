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
