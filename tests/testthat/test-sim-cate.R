.test_aipw_score <- function(dat, e_hat, m0_hat, m1_hat) {
  pmax_e <- pmax(e_hat, 1e-8)
  pmax_1e <- pmax(1 - e_hat, 1e-8)
  m1_hat - m0_hat +
    dat$D * (dat$Y - m1_hat) / pmax_e -
    (1 - dat$D) * (dat$Y - m0_hat) / pmax_1e
}

.test_fit_aipw_nuisances <- function(dat, x_cols, train = rep(TRUE, nrow(dat)),
                                     clip = 1e-3) {
  e_fmla <- stats::as.formula(paste("D ~", paste(x_cols, collapse = " + ")))
  m_fmla <- stats::as.formula(paste("Y ~", paste(x_cols, collapse = " + ")))
  train_dat <- dat[train, , drop = FALSE]

  e_fit <- stats::glm(e_fmla, data = train_dat, family = stats::binomial())
  e_hat <- stats::predict(e_fit, newdata = dat, type = "response")
  e_hat <- pmin(pmax(e_hat, clip), 1 - clip)

  m0_fit <- stats::lm(m_fmla, data = train_dat[train_dat$D == 0L, , drop = FALSE])
  m1_fit <- stats::lm(m_fmla, data = train_dat[train_dat$D == 1L, , drop = FALSE])
  list(
    e_hat = e_hat,
    m0_hat = stats::predict(m0_fit, newdata = dat),
    m1_hat = stats::predict(m1_fit, newdata = dat)
  )
}

test_that("sim_cate_data: output schema is correct for each design", {
  designs <- c("dgp1_linear", "dgp2_nonlinear", "dgp3_sparse_overlap",
               "dgp4_clustered", "dgp5_correlated", "dgp6_nuisance_efficiency")
  expected_cols <- c(
    "X1", "X2", "D", "Y",
    "e_oracle", "m_oracle", "m0_oracle", "m1_oracle",
    "tau_oracle", "rho_oracle", "aipw_oracle"
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
    expect_named(
      dgp$truth,
      c("tau_fn", "e_fn", "m0_fn", "sigma", "score_sd_x0", "m_fn")
    )
    expect_equal(dgp$truth$sigma, 0.5)
    expect_true(is.finite(dgp$truth$score_sd_x0))
    expect_gt(dgp$truth$score_sd_x0, 0)
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

test_that("sim_cate_data: aipw_oracle matches the dormant-source AIPW score", {
  dgp <- sim_cate_data("dgp2_nonlinear", n = 100L, kx = 2L, seed = 3L)
  dat <- dgp$data
  aipw_manual <- .test_aipw_score(
    dat = dat,
    e_hat = dat$e_oracle,
    m0_hat = dat$m0_oracle,
    m1_hat = dat$m1_oracle
  )

  expect_equal(dat$aipw_oracle, aipw_manual, tolerance = 1e-12)
  expect_false(isTRUE(all.equal(dat$aipw_oracle, dat$rho_oracle)))
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
    "x0.*NA"
  )
  expect_error(
    sim_cate_data("dgp1_linear", n = 10L, kx = 2L, x0 = c(0, Inf)),
    "x0.*non-finite"
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

test_that("grf-missing benchmark warnings are emitted without call metadata", {
  skip_if(requireNamespace("grf", quietly = TRUE), "grf installed; missing-package warning path not active")

  m_cf <- getFromNamespace(".method_causal_forest", "RdistNN")
  m_llf <- getFromNamespace(".method_llf", "RdistNN")
  dgp <- sim_cate_data("dgp1_linear", n = 40L, kx = 2L, seed = 11L)
  x_cols <- grep("^X", names(dgp$data), value = TRUE)

  w_cf <- NULL
  withCallingHandlers(
    m_cf(dgp, x_cols = x_cols, alpha = 0.05),
    warning = function(w) {
      w_cf <<- w
      invokeRestart("muffleWarning")
    }
  )
  expect_s3_class(w_cf, "warning")
  expect_match(conditionMessage(w_cf), "grf not available; causal_forest returns NA\\.")
  expect_null(conditionCall(w_cf))

  w_llf <- NULL
  withCallingHandlers(
    m_llf(dgp, x_cols = x_cols, alpha = 0.05),
    warning = function(w) {
      w_llf <<- w
      invokeRestart("muffleWarning")
    }
  )
  expect_s3_class(w_llf, "warning")
  expect_match(conditionMessage(w_llf), "grf not available; llf returns NA\\.")
  expect_null(conditionCall(w_llf))
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

test_that("sim_cate_data dgp2/dgp3/dgp4/dgp6: default x0 is (1, 0, ...)", {
  for (design in c("dgp2_nonlinear", "dgp3_sparse_overlap",
                   "dgp4_clustered", "dgp6_nuisance_efficiency")) {
    dgp <- sim_cate_data(design, n = 50L, kx = 2L, seed = 1L)
    expect_equal(dgp$x0[1L], 1.0, tolerance = 1e-12)
    expect_equal(dgp$x0[2L], 0.0, tolerance = 1e-12)
  }
})

test_that("sim_cate_data dgp4_clustered: mixture creates a dense spike near x0", {
  dgp <- sim_cate_data("dgp4_clustered", n = 1200L, kx = 2L, seed = 7L)
  d2 <- (dgp$data$X1 - dgp$x0[1L])^2 + (dgp$data$X2 - dgp$x0[2L])^2
  prop_near <- mean(d2 < 0.6^2)
  expect_gt(prop_near, 0.7)
})

test_that("sim_cate_data dgp4_clustered: both remote clusters are populated", {
  dgp <- sim_cate_data("dgp4_clustered", n = 1500L, kx = 2L, seed = 8L)
  pos_cluster <- mean(dgp$data$X1 > 2.5 & dgp$data$X2 > 0.5)
  neg_cluster <- mean(dgp$data$X1 < -0.5 & dgp$data$X2 < -0.5)
  expect_gt(pos_cluster, 0.04)
  expect_gt(neg_cluster, 0.04)
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

# ── DGP 6: nuisance-efficiency design ────────────────────────────────────────

test_that("sim_cate_data dgp6_nuisance_efficiency: propensity is clipped to (clip, 1-clip)", {
  dgp <- sim_cate_data("dgp6_nuisance_efficiency", n = 200L, kx = 2L, seed = 1L, clip = 0.05)
  expect_true(all(dgp$data$e_oracle >= 0.05))
  expect_true(all(dgp$data$e_oracle <= 0.95))
})

test_that("sim_cate_data dgp6_nuisance_efficiency: oracle marginal mean is linear", {
  dgp <- sim_cate_data("dgp6_nuisance_efficiency", n = 200L, kx = 2L, seed = 4L)
  m_linear <- 1 + 1.5 * (dgp$data$X1 - 1) - dgp$data$X2
  expect_equal(dgp$data$m_oracle, m_linear, tolerance = 1e-12)
})

test_that("sim_cate_data dgp6_nuisance_efficiency: tau_x0 matches tau_fn at x0", {
  dgp <- sim_cate_data("dgp6_nuisance_efficiency", n = 50L, kx = 2L, seed = 2L)
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

test_that("sim_cate_one_rep: dnn_fullsample returns finite metrics", {
  r <- sim_cate_one_rep(
    n = 120L, kx = 2L, k = 3L,
    s = floor(120L^(3/4)),
    design = "dgp2_nonlinear",
    methods = "dnn_fullsample",
    seed = 42L
  )
  expect_equal(nrow(r), 1L)
  expect_equal(r$method, "dnn_fullsample")
  expect_true(is.finite(r$tau_hat))
  expect_true(r$se > 0)
  expect_true(r$cover %in% c(0, 1))
})

test_that("sim_cate_one_rep: dnn_oracle matches direct oracle AIPW fit", {
  n <- 120L
  kx <- 2L
  k <- 3L
  s <- floor(n^(3/4))
  seed <- 42L
  dgp <- sim_cate_data("dgp2_nonlinear", n = n, kx = kx, seed = seed)

  fit <- cate_dnn(
    data = dgp$data,
    x_cols = paste0("X", seq_len(kx)),
    y_col = "Y",
    t_col = "D",
    e_col = "e_oracle",
    m0_col = "m0_oracle",
    m1_col = "m1_oracle",
    x0 = dgp$x0,
    k = k,
    s = s,
    se = TRUE,
    jackknife = "fast",
    return = "ess"
  )
  r <- sim_cate_one_rep(
    n = n, kx = kx, k = k,
    s = s,
    design = "dgp2_nonlinear",
    methods = "dnn_oracle",
    seed = seed
  )

  expect_equal(nrow(r), 1L)
  expect_equal(r$method, "dnn_oracle")
  expect_equal(r$tau_hat, fit$tau, tolerance = 1e-12)
  expect_equal(r$se, fit$se, tolerance = 1e-12)
  expect_equal(r$neff, fit$ess, tolerance = 1e-12)
})

test_that("sim_cate_one_rep: dnn_fullsample matches direct full-sample AIPW fit", {
  n <- 120L
  kx <- 2L
  k <- 3L
  s <- floor(n^(3/4))
  seed <- 42L
  clip <- 1e-3
  dgp <- sim_cate_data("dgp2_nonlinear", n = n, kx = kx, seed = seed, clip = clip)
  x_cols <- paste0("X", seq_len(kx))
  dat <- dgp$data

  nuis <- .test_fit_aipw_nuisances(dat, x_cols = x_cols, clip = clip)
  aipw_hat <- .test_aipw_score(
    dat = dat,
    e_hat = nuis$e_hat,
    m0_hat = nuis$m0_hat,
    m1_hat = nuis$m1_hat
  )

  fit <- cate_dnn(
    X = as.matrix(dat[, x_cols]),
    x0 = dgp$x0,
    pseudo_outcome = aipw_hat,
    k = k,
    s = s,
    se = TRUE,
    jackknife = "fast",
    return = "ess"
  )
  r <- sim_cate_one_rep(
    n = n, kx = kx, k = k,
    s = s,
    design = "dgp2_nonlinear",
    methods = "dnn_fullsample",
    seed = seed,
    clip = clip
  )

  expect_equal(nrow(r), 1L)
  expect_equal(r$method, "dnn_fullsample")
  expect_equal(r$tau_hat, fit$tau, tolerance = 1e-12)
  expect_equal(r$se, fit$se, tolerance = 1e-12)
  expect_equal(r$neff, fit$ess, tolerance = 1e-12)
})

test_that("sim_cate_one_rep: dnn_crossfit matches direct cross-fitted AIPW fit", {
  n <- 120L
  kx <- 2L
  k <- 3L
  s <- floor(n^(3/4))
  seed <- 42L
  clip <- 1e-3
  n_folds <- 5L

  dgp <- sim_cate_data("dgp2_nonlinear", n = n, kx = kx, seed = seed, clip = clip)
  x_cols <- paste0("X", seq_len(kx))
  dat <- dgp$data

  fold_ids <- sample(rep_len(seq_len(n_folds), n))
  e_cf <- m0_cf <- m1_cf <- numeric(n)
  for (fold_k in seq_len(n_folds)) {
    te <- fold_ids == fold_k
    nuis <- .test_fit_aipw_nuisances(
      dat = dat,
      x_cols = x_cols,
      train = !te,
      clip = clip
    )
    e_cf[te] <- nuis$e_hat[te]
    m0_cf[te] <- nuis$m0_hat[te]
    m1_cf[te] <- nuis$m1_hat[te]
  }
  aipw_cf <- .test_aipw_score(
    dat = dat,
    e_hat = e_cf,
    m0_hat = m0_cf,
    m1_hat = m1_cf
  )

  fit <- cate_dnn(
    X = as.matrix(dat[, x_cols]),
    x0 = dgp$x0,
    pseudo_outcome = aipw_cf,
    k = k,
    s = s,
    se = TRUE,
    jackknife = "fast",
    return = "ess"
  )
  r <- sim_cate_one_rep(
    n = n, kx = kx, k = k,
    s = s,
    design = "dgp2_nonlinear",
    methods = "dnn_crossfit",
    seed = seed,
    clip = clip,
    n_folds = n_folds
  )

  expect_equal(nrow(r), 1L)
  expect_equal(r$method, "dnn_crossfit")
  expect_equal(r$tau_hat, fit$tau, tolerance = 1e-12)
  expect_equal(r$se, fit$se, tolerance = 1e-12)
  expect_equal(r$neff, fit$ess, tolerance = 1e-12)
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

test_that("sim_cate_one_rep: all six core methods run together without error", {
  r <- sim_cate_one_rep(
    n = 100L, kx = 2L, k = 3L,
    s = floor(100L^(3/4)),
    design = "dgp2_nonlinear",
    methods = c("dnn_oracle", "dnn_oracle_mahal", "dnn_fullsample", "dnn_crossfit",
                "local_constant", "ols_interact"),
    seed = 7L
  )
  expect_equal(nrow(r), 6L)
  expect_equal(r$method,
               c("dnn_oracle", "dnn_oracle_mahal", "dnn_fullsample", "dnn_crossfit",
                 "local_constant", "ols_interact"))
  expect_true(all(is.finite(r$tau_hat)))
  expect_true(all(r$se > 0))
})

test_that("sim_cate_one_rep exposes the theorem-facing DNN decomposition", {
  r <- sim_cate_one_rep(
    n = 160L, kx = 2L, k = 1L, s = 24L,
    design = "dgp2_nonlinear",
    methods = c("dnn_oracle", "dnn_fullsample", "dnn_crossfit", "local_constant"),
    seed = 706L
  )
  diagnostics <- c(
    "localized_truth", "branch_drift", "estimator_remainder",
    "studentized_error", "score_sd0", "standardized_branch_drift",
    "physical_radius_90", "treated_ess", "control_ess", "top5_score_share"
  )
  expect_true(all(diagnostics %in% names(r)))
  expect_true(all(is.finite(as.matrix(r[, diagnostics]))))
  expect_equal(r$bias, r$branch_drift + r$estimator_remainder, tolerance = 1e-12)
  expect_equal(r$studentized_error, r$bias / r$se, tolerance = 1e-12)

  dnn <- r[grepl("^dnn_", r$method), , drop = FALSE]
  expect_equal(length(unique(dnn$localized_truth)), 1L)
  expect_equal(length(unique(dnn$branch_drift)), 1L)
  expect_equal(length(unique(dnn$physical_radius_90)), 1L)
  expect_equal(length(unique(dnn$treated_ess)), 1L)
  expect_equal(length(unique(dnn$control_ess)), 1L)
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
