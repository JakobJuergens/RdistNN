skip_if_not(
  all(vapply(
    c("sim_cate_one_rep", "sim_cate_run", "summarise_sim_results", "new_tuning"),
    exists, logical(1),
    mode = "function"
  )),
  "Layer 8 simulation-grid dependencies not found; skipping file"
)

.layer8_tuning <- function(default_s = 6L, s_grid = c(4L, 6L, 8L)) {
  new_tuning(
    name = "layer8_grid",
    eval = function(...) list(s = default_s, meta = list(kappa = 1)),
    eval_grid = function(...) {
      data.frame(
        kappa = c(0.75, 1, 1.25),
        s = s_grid,
        stringsAsFactors = FALSE
      )
    }
  )
}

.layer8_raw_schema <- c(
  "method", "design", "n", "kx", "k", "tuning_source",
  "candidate_id", "kappa", "s", "seed",
  "tau_hat", "tau0", "bias", "sq_err", "se", "cover", "ci_len", "neff",
  "status"
)

test_that("sim_cate_one_rep grid mode emits DNN candidate rows matching scalar runs", {
  tun <- .layer8_tuning()
  n <- 64L
  kx <- 2L
  k <- 1L
  seed <- 8801L

  grid <- sim_cate_one_rep(
    n = n, kx = kx, k = k,
    design = "dgp2_nonlinear",
    methods = "dnn_oracle",
    seed = seed,
    tuning = tun,
    tuning_mode = "grid"
  )

  expect_s3_class(grid, "data.frame")
  expect_equal(nrow(grid), 3L)
  expect_true(all(c("candidate_id", "kappa", "s", "status") %in% names(grid)))
  expect_equal(grid$method, rep("dnn_oracle", 3L))
  expect_equal(grid$status, rep("ok", 3L))
  expect_false(anyNA(grid$candidate_id))

  for (i in seq_len(nrow(grid))) {
    scalar <- sim_cate_one_rep(
      n = n, kx = kx, k = k, s = grid$s[[i]],
      design = "dgp2_nonlinear",
      methods = "dnn_oracle",
      seed = seed
    )
    expect_true(is.na(scalar$candidate_id))
    expect_equal(grid$tau_hat[[i]], scalar$tau_hat[[1L]], tolerance = 0)
    expect_equal(grid$se[[i]], scalar$se[[1L]], tolerance = 0)
    expect_equal(grid$neff[[i]], scalar$neff[[1L]], tolerance = 0)
  }
})

test_that("sim_cate_one_rep grid mode reuses one cross-fit score across candidates", {
  tun <- .layer8_tuning()
  n <- 72L
  kx <- 2L
  k <- 1L
  seed <- 8802L

  grid <- sim_cate_one_rep(
    n = n, kx = kx, k = k,
    design = "dgp2_nonlinear",
    methods = "dnn_crossfit",
    seed = seed,
    n_folds = 4L,
    tuning = tun,
    tuning_mode = "grid"
  )

  expect_equal(nrow(grid), 3L)
  expect_equal(grid$status, rep("ok", 3L))

  for (i in seq_len(nrow(grid))) {
    scalar <- sim_cate_one_rep(
      n = n, kx = kx, k = k, s = grid$s[[i]],
      design = "dgp2_nonlinear",
      methods = "dnn_crossfit",
      seed = seed,
      n_folds = 4L
    )
    expect_equal(grid$tau_hat[[i]], scalar$tau_hat[[1L]], tolerance = 0)
    expect_equal(grid$se[[i]], scalar$se[[1L]], tolerance = 0)
    expect_equal(grid$neff[[i]], scalar$neff[[1L]], tolerance = 0)
  }
})

test_that("sim_cate_run grid mode uses raw schema and does not duplicate benchmarks", {
  raw <- sim_cate_run(
    designs = "dgp1_linear",
    n_vals = 60L,
    kx_vals = 2L,
    n_reps = 2L,
    master_seed = 8803L,
    k = 1L,
    tuning = .layer8_tuning(),
    tuning_mode = "grid",
    methods = c("dnn_oracle", "ols_interact"),
    verbose = FALSE
  )

  expect_s3_class(raw, "data.frame")
  expect_equal(names(raw), .layer8_raw_schema)
  expect_equal(nrow(raw), 8L)
  expect_equal(unique(raw$tuning_source), "tuning")

  dnn <- raw[raw$method == "dnn_oracle", , drop = FALSE]
  bench <- raw[raw$method == "ols_interact", , drop = FALSE]
  expect_equal(nrow(dnn), 6L)
  expect_equal(nrow(bench), 2L)
  expect_equal(as.integer(table(dnn$seed)), rep(3L, 2L))
  expect_true(all(is.na(bench$candidate_id)))
  expect_true(all(is.na(bench$kappa)))
  expect_true(all(is.na(bench$s)))
  expect_equal(bench$status, rep("ok", 2L))
})

test_that("sim_cate_run grid mode records infeasible candidates without duplicating benchmarks", {
  raw <- NULL
  expect_warning(
    raw <- sim_cate_run(
      designs = "dgp1_linear",
      n_vals = 50L,
      kx_vals = 2L,
      n_reps = 1L,
      master_seed = 8804L,
      k = 2L,
      tuning = .layer8_tuning(default_s = 4L, s_grid = c(1L, 4L, 6L)),
      tuning_mode = "grid",
      methods = c("dnn_oracle", "ols_interact"),
      verbose = FALSE
    ),
    "Skipped infeasible"
  )

  dnn <- raw[raw$method == "dnn_oracle", , drop = FALSE]
  bench <- raw[raw$method == "ols_interact", , drop = FALSE]
  expect_equal(nrow(dnn), 3L)
  expect_equal(nrow(bench), 1L)
  expect_equal(sum(dnn$status == "infeasible"), 1L)
  expect_true(all(is.na(dnn$tau_hat[dnn$status == "infeasible"])))
  expect_true(all(is.na(dnn$se[dnn$status == "infeasible"])))
  expect_true(is.na(bench$s))
})

test_that("sim_cate_run grid mode shapes total-failure rows by method and candidate", {
  raw <- NULL
  expect_warning(
    raw <- sim_cate_run(
      designs = "bad_design",
      n_vals = 50L,
      kx_vals = 2L,
      n_reps = 1L,
      master_seed = 8805L,
      k = 1L,
      tuning = .layer8_tuning(),
      tuning_mode = "grid",
      methods = c("dnn_oracle", "ols_interact"),
      verbose = FALSE
    ),
    "rep seed"
  )

  expect_equal(names(raw), .layer8_raw_schema)
  expect_equal(nrow(raw), 4L)
  expect_equal(sum(raw$method == "dnn_oracle"), 3L)
  expect_equal(sum(raw$method == "ols_interact"), 1L)
  expect_equal(raw$status, rep("failed", 4L))
  expect_true(all(is.na(raw$tau_hat)))
  expect_true(all(is.na(raw$s[raw$method == "ols_interact"])))
})

test_that("summarise_sim_results groups by candidate fields with NA-safe benchmark groups", {
  raw <- data.frame(
    method = c("dnn_oracle", "dnn_oracle", "dnn_oracle", "dnn_oracle", "ols_interact"),
    design = "dgp1_linear",
    n = 50L,
    kx = 2L,
    k = 1L,
    tuning_source = "tuning",
    candidate_id = c("kappa_0p75", "kappa_0p75", "kappa_0p75", "kappa_1", NA_character_),
    kappa = c(0.75, 0.75, 0.75, 1, NA_real_),
    s = c(4L, 4L, 4L, 6L, NA_integer_),
    seed = c(1L, 2L, 3L, 1L, 1L),
    tau_hat = c(1, NA, 3, 2, 3),
    tau0 = 1,
    bias = c(0, NA, 2, 1, 2),
    sq_err = c(0, NA, 4, 1, 4),
    se = c(0.2, NA, 0.4, 0.3, 0.4),
    cover = c(1, NA, 1, 0, 0),
    ci_len = c(0.4, NA, 0.8, 0.6, 0.8),
    neff = c(10, NA, 14, 12, NA),
    status = c("ok", "failed", "ok", "ok", "ok"),
    stringsAsFactors = FALSE
  )

  tab <- summarise_sim_results(raw)
  expect_equal(nrow(tab), 3L)
  expect_true(all(c("candidate_id", "kappa", "s") %in% names(tab)))

  cand_075 <- tab[!is.na(tab$candidate_id) & tab$candidate_id == "kappa_0p75", , drop = FALSE]
  cand_1 <- tab[!is.na(tab$candidate_id) & tab$candidate_id == "kappa_1", , drop = FALSE]
  bench <- tab[is.na(tab$candidate_id), , drop = FALSE]

  expect_equal(cand_075$n_reps, 3L)
  expect_equal(cand_075$n_failed, 1L)
  expect_equal(cand_075$bias, 1)
  expect_equal(cand_075$emp_sd, stats::sd(c(1, 3)))
  expect_equal(cand_075$se_to_sd, mean(c(0.2, 0.4)) / stats::sd(c(1, 3)))
  expect_equal(cand_1$n_reps, 1L)
  expect_true(is.na(cand_1$emp_sd))
  expect_true(is.na(cand_1$se_to_sd))
  expect_equal(bench$method, "ols_interact")
  expect_true(is.na(bench$s))
})

test_that("summarise_sim_results returns NA rather than NaN for all-failed groups", {
  raw <- data.frame(
    method = c("dnn_oracle", "dnn_oracle", "ols_interact"),
    design = "dgp1_linear",
    n = 50L,
    kx = 2L,
    k = 1L,
    tuning_source = "tuning",
    candidate_id = c("kappa_bad", "kappa_bad", NA_character_),
    kappa = c(9, 9, NA_real_),
    s = c(60L, 60L, NA_integer_),
    seed = c(1L, 2L, 1L),
    tau_hat = NA_real_,
    tau0 = 1,
    bias = NA_real_,
    sq_err = NA_real_,
    se = NA_real_,
    cover = NA_real_,
    ci_len = NA_real_,
    neff = NA_real_,
    status = c("failed", "infeasible", "failed"),
    stringsAsFactors = FALSE
  )

  tab <- summarise_sim_results(raw)
  expect_equal(nrow(tab), 2L)
  expect_equal(tab$n_reps, c(2L, 1L))
  expect_equal(tab$n_failed, c(2L, 1L))

  metrics <- c("bias", "rmse", "se_mean", "cover", "ci_len", "emp_sd", "se_to_sd", "neff")
  metric_values <- unlist(tab[metrics], use.names = FALSE)
  expect_true(all(is.na(metric_values)))
  expect_false(any(is.nan(metric_values)))
})

test_that("summarise_sim_results excludes non-ok rows even if metrics are present", {
  raw <- data.frame(
    method = "dnn_oracle",
    design = "dgp1_linear",
    n = 50L,
    kx = 2L,
    k = 1L,
    tuning_source = "tuning",
    candidate_id = "kappa_1",
    kappa = 1,
    s = 6L,
    seed = c(1L, 2L),
    tau_hat = c(1, 100),
    tau0 = 1,
    bias = c(0, 99),
    sq_err = c(0, 9801),
    se = c(0.2, 99),
    cover = c(1, 0),
    ci_len = c(0.4, 198),
    neff = c(10, 1),
    status = c("ok", "failed"),
    stringsAsFactors = FALSE
  )

  tab <- summarise_sim_results(raw)
  expect_equal(tab$n_reps, 2L)
  expect_equal(tab$n_failed, 1L)
  expect_equal(tab$bias, 0)
  expect_equal(tab$rmse, 0)
  expect_equal(tab$se_mean, 0.2)
  expect_equal(tab$cover, 1)
  expect_equal(tab$ci_len, 0.4)
  expect_true(is.na(tab$emp_sd))
  expect_true(is.na(tab$se_to_sd))
  expect_equal(tab$neff, 10)
})

test_that("grid mode requires a tuning object", {
  expect_error(
    sim_cate_one_rep(
      n = 50L, kx = 2L, k = 1L,
      design = "dgp1_linear",
      tuning_mode = "grid"
    ),
    "requires `tuning`"
  )
  expect_error(
    sim_cate_run(
      designs = "dgp1_linear",
      n_vals = 50L,
      kx_vals = 2L,
      n_reps = 1L,
      tuning_mode = "grid",
      verbose = FALSE
    ),
    "requires `tuning`"
  )
})
