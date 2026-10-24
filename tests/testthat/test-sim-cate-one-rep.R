# Single-replicate estimators and direct-fit equivalence.
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
  withr::local_preserve_seed()
  n <- 120L
  kx <- 2L
  k <- 3L
  s <- floor(n^(3/4))
  seed <- 42L
  clip <- 1e-3
  n_folds <- 5L

  # Match the replicate's continuous data-then-fold stream explicitly.
  set.seed(seed)
  dgp <- sim_cate_data("dgp2_nonlinear", n = n, kx = kx, seed = NULL, clip = clip)
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
    "localized_truth", "localization_drift", "estimator_remainder",
    "studentized_error", "score_sd0", "standardized_localization_drift",
    "physical_radius_90", "treated_ess", "control_ess", "top5_score_share"
  )
  expect_true(all(diagnostics %in% names(r)))
  expect_true(all(is.finite(as.matrix(r[, diagnostics]))))
  expect_equal(r$bias, r$localization_drift + r$estimator_remainder, tolerance = 1e-12)
  expect_equal(r$studentized_error, r$bias / r$se, tolerance = 1e-12)

  dnn <- r[grepl("^dnn_", r$method), , drop = FALSE]
  expect_equal(length(unique(dnn$localized_truth)), 1L)
  expect_equal(length(unique(dnn$localization_drift)), 1L)
  expect_equal(length(unique(dnn$physical_radius_90)), 1L)
  expect_equal(length(unique(dnn$treated_ess)), 1L)
  expect_equal(length(unique(dnn$control_ess)), 1L)
})
test_that("simulation method rows require the current explicit status and schema", {
  metric <- getFromNamespace(".sim_metric", "RdistNN")
  validate <- getFromNamespace(".sim_cate_normalize_method_rows", "RdistNN")
  row <- metric(tau_hat = 1, se_hat = .2, tau0 = 0, alpha = .05)
  expect_identical(validate(row), row)
  for (field in c("status", "tau_hat", "method_variant", "w_hat_supplied")) {
    bad <- row[, setdiff(names(row), field), drop = FALSE]
    expect_error(validate(bad), paste0("missing required column\\(s\\): ", field))
  }
  expect_error(validate(as.list(row)), "must be a data.frame")
  row$status <- "failed"
  expect_identical(validate(row)$status, "failed")
})

test_that("estimated CATE simulation signals honor a requested tiny clip", {
  .test_rng_scope({
    RNGkind("Mersenne-Twister", "Inversion", "Rejection")
    clip <- 1e-10
    cases <- list(
      list(method = "dnn_fullsample", n = 8L, seed = 1L, s = 3L),
      list(method = "dnn_crossfit", n = 24L, seed = 3L, s = 4L)
    )
    for (case in cases) {
      # Restart the actual data-then-fold stream. These small samples deliberately
      # produce separated propensity fits, so suppress their expected GLM warnings.
      set.seed(case$seed)
      dgp <- sim_cate_data("dgp3_sparse_overlap", n = case$n, kx = 2L,
        seed = NULL, clip = clip)
      dat <- dgp$data
      x_cols <- c("X1", "X2")
      if (case$method == "dnn_fullsample") {
        nuis <- suppressWarnings(.test_fit_aipw_nuisances(dat, x_cols, clip = clip))
      } else {
        fold_ids <- sample(rep_len(1:4, case$n))
        nuis <- list(e_hat = numeric(case$n), m0_hat = numeric(case$n),
          m1_hat = numeric(case$n))
        for (fold in 1:4) {
          test <- fold_ids == fold
          trained <- suppressWarnings(.test_fit_aipw_nuisances(
            dat, x_cols, train = !test, clip = clip))
          for (name in names(nuis)) nuis[[name]][test] <- trained[[name]][test]
        }
      }
      expect_true(any(nuis$e_hat < 1e-8 | nuis$e_hat > 1 - 1e-8))
      # Independent score formula uses the requested-clipped predictions directly;
      # the old generic test score helper imposes its own 1e-8 floor.
      signal <- nuis$m1_hat - nuis$m0_hat +
        dat$D * (dat$Y - nuis$m1_hat) / nuis$e_hat -
        (1 - dat$D) * (dat$Y - nuis$m0_hat) / (1 - nuis$e_hat)
      expected <- cate_dnn(X = as.matrix(dat[, x_cols]), x0 = dgp$x0,
        pseudo_outcome = signal, k = 1L, s = case$s, se = TRUE)
      actual <- suppressWarnings(sim_cate_one_rep(
        n = case$n, kx = 2L, k = 1L, s = case$s,
        design = "dgp3_sparse_overlap", methods = case$method,
        seed = case$seed, clip = clip, n_folds = 4L))
      expect_identical(actual$status, "ok")
      expect_equal(actual$tau_hat, expected$tau, tolerance = 1e-12)
      expect_equal(actual$se, expected$se, tolerance = 1e-12)
    }
  })
})

test_that("invalid CATE cross-fit fold counts propagate argument errors", {
  for (folds in c(0, .5, 2.5, 41)) {
    expect_error(sim_cate_one_rep(
      n = 40L, kx = 2L, k = 1L, s = 4L,
      design = "dgp2_nonlinear", methods = c("dnn_oracle", "dnn_crossfit"),
      seed = 42L, n_folds = folds), "n_folds", class = "rdistnn_input_error")
  }
})

test_that("CATE cross-fitting accepts two folds through leave-one-out folds", {
  for (folds in c(2L, 24L)) {
    actual <- sim_cate_one_rep(
      n = 24L, kx = 2L, k = 1L, s = 4L,
      design = "dgp1_linear", methods = "dnn_crossfit",
      seed = 42L, n_folds = folds)
    expect_identical(actual$status, "ok")
    expect_true(is.finite(actual$tau_hat))
    expect_true(is.finite(actual$se) && actual$se > 0)
  }
})
