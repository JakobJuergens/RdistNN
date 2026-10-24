test_that("diagnostics(rdistnn_cate_result) reports overlap when treatment state is available", {
  .test_rng_scope({
    demo <- cate_demo_data(n = 70L, p = 2L, design = "hetero", sigma = 0.2, seed = 20L)
    fit <- cate_dnn(
      data = demo$data,
      x_cols = c("X1", "X2"),
      y_col = "Y",
      t_col = "D",
      m0_col = "m0_oracle", m1_col = "m1_oracle",
      e_col = "e_oracle",
      x0 = demo$x0,
      method = "exact",
      k = 2L,
      s = 4L,
      se = FALSE,
      return = "W"
    )
    d <- diagnostics(fit)

    expect_s3_class(d, "rdistnn_diagnostics")
    expect_true(all(d$overlap_available))
    expect_true(all(d$propensity_available))

    t_vec <- demo$data$D
    e <- demo$data$e_oracle
    expect_equal(d$treated_mass, as.numeric(crossprod(t_vec, fit$W)), tolerance = 1e-12)
    expect_equal(d$control_mass, as.numeric(crossprod(1 - t_vec, fit$W)), tolerance = 1e-12)
    expect_equal(
      d$weighted_mean_propensity,
      as.numeric(crossprod(e, fit$W)) / colSums(fit$W),
      tolerance = 1e-12
    )
    expect_equal(d$treated_share, d$treated_mass / (d$treated_mass + d$control_mass), tolerance = 1e-12)
    expect_equal(
      d$weighted_propensity_q10,
      vapply(seq_len(ncol(fit$W)), function(j) .diag_weighted_quantile(e, fit$W[, j], prob = 0.1), numeric(1L)),
      tolerance = 1e-12
    )
    expect_equal(
      d$weighted_propensity_q90,
      vapply(seq_len(ncol(fit$W)), function(j) .diag_weighted_quantile(e, fit$W[, j], prob = 0.9), numeric(1L)),
      tolerance = 1e-12
    )
  }, seed = 777L)
})

test_that("diagnostics(rdistnn_cate_result) requires retained AIPW propensity state", {
  .test_rng_scope({
    demo <- cate_demo_data(n = 60L, p = 2L, design = "hetero", sigma = 0.2, seed = 22L)
    fit <- cate_dnn(
      data = demo$data,
      x_cols = c("X1", "X2"),
      y_col = "Y",
      t_col = "D",
      m0_col = "m0_oracle", m1_col = "m1_oracle",
      e_col = "e_oracle",
      x0 = demo$x0,
      method = "exact",
      k = 2L,
      s = 4L,
      se = FALSE,
      return = "W"
    )
    fit$.model["e"] <- NULL
    expect_error(diagnostics(fit), "e.*[Rr]efit")
    fit$.model["e"] <- list(NULL)
    expect_error(diagnostics(fit), "e.*[Rr]efit")
  }, seed = 777L)
})

test_that("diagnostics(rdistnn_cate_result) validates retained treatment and propensity state", {
  .test_rng_scope({
    demo <- cate_demo_data(n = 60L, p = 2L, design = "hetero", sigma = 0.2, seed = 23L)
    fit <- cate_dnn(
      data = demo$data,
      x_cols = c("X1", "X2"),
      y_col = "Y",
      t_col = "D",
      m0_col = "m0_oracle", m1_col = "m1_oracle",
      e_col = "e_oracle",
      x0 = demo$x0,
      method = "exact",
      k = 2L,
      s = 4L,
      se = FALSE,
      return = "W"
    )

    bad_length <- fit
    bad_length$.model$t_vec <- bad_length$.model$t_vec[-1L]
    expect_error(diagnostics(bad_length), "t_vec.*[Rr]efit")

    bad_binary <- fit
    bad_binary$.model$t_vec[1L] <- 2
    expect_error(diagnostics(bad_binary), "t_vec.*[Rr]efit")

    bad_e_length <- fit
    bad_e_length$.model$e <- bad_e_length$.model$e[-1L]
    expect_error(diagnostics(bad_e_length), "e.*[Rr]efit")

    bad_e_na <- fit
    bad_e_na$.model$e[1L] <- NA_real_
    expect_error(diagnostics(bad_e_na), "e.*[Rr]efit")
  }, seed = 777L)
})

test_that("diagnostics(rdistnn_cate_result) marks overlap unavailable on pseudo-outcome-only path", {
  .test_rng_scope({
    n <- 45L
    x_mat <- matrix(rnorm(n * 2L), n, 2L)
    fit <- cate_dnn(
      X = x_mat,
      x0 = c(0, 0),
      pseudo_outcome = rnorm(n),
      method = "exact",
      k = 2L,
      s = 4L,
      se = FALSE,
      return = "W"
    )
    d <- diagnostics(fit)

    expect_true(all(!d$overlap_available))
    expect_true(all(!d$propensity_available))
    expect_true(all(is.na(d$treated_mass)))
    expect_true(all(is.na(d$weighted_mean_propensity)))
  }, seed = 777L)
})
