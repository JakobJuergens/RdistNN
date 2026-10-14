# CATE simulation designs, truth identities, and input validation.
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
