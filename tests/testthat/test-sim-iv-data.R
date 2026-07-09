test_that("sim_iv_data: output schema is correct for each design", {
  designs <- c(
    "dgp1_iv_justid_linear",
    "dgp2_iv_homog_overid",
    "dgp3_iv_heterog_overid",
    "dgp4_iv_weak_local",
    "dgp5_iv_correlated_instruments",
    "dgp6_iv_nonlinear_firststage"
  )

  for (des in designs) {
    dgp <- sim_iv_data(des, n = 80L, kx = 2L, q = 2L, seed = 1L)
    q <- dgp$q
    expected_cols <- c(
      "X1", "X2", "Y", "D",
      paste0("Z", seq_len(q)),
      paste0("piZ", seq_len(q), "_oracle"),
      "muY_oracle", "muD_oracle", "theta_oracle"
    )

    expect_type(dgp, "list")
    expect_named(dgp, c("data", "x0", "q", "theta_x0", "truth"))
    expect_s3_class(dgp$data, "data.frame")
    expect_equal(nrow(dgp$data), 80L)
    expect_true(all(expected_cols %in% names(dgp$data)), info = des)
    expect_length(dgp$x0, 2L)
    expect_true(is.finite(dgp$theta_x0))
    expect_true(all(is.finite(unlist(dgp$data[, c(
      paste0("piZ", seq_len(q), "_oracle"),
      "muY_oracle", "muD_oracle", "theta_oracle"
    )]))))
  }
})

test_that("sim_iv_data: deterministic with fixed seed", {
  dgp1 <- sim_iv_data("dgp2_iv_homog_overid", n = 100L, kx = 2L, q = 2L, seed = 7L)
  dgp2 <- sim_iv_data("dgp2_iv_homog_overid", n = 100L, kx = 2L, q = 2L, seed = 7L)
  expect_equal(dgp1$data, dgp2$data)
  expect_equal(dgp1$theta_x0, dgp2$theta_x0)
})

test_that("sim_iv_data: different seeds produce different data", {
  dgp1 <- sim_iv_data("dgp1_iv_justid_linear", n = 80L, kx = 2L, seed = 1L)
  dgp2 <- sim_iv_data("dgp1_iv_justid_linear", n = 80L, kx = 2L, seed = 2L)
  expect_false(identical(dgp1$data$Y, dgp2$data$Y))
})

test_that("sim_iv_data: just-identified design forces q = 1", {
  dgp <- sim_iv_data("dgp1_iv_justid_linear", n = 60L, kx = 2L, q = 3L, seed = 1L)
  expect_equal(dgp$q, 1L)
  expect_true(all(c("Z1", "piZ1_oracle") %in% names(dgp$data)))
  expect_false("Z2" %in% names(dgp$data))
})

test_that("sim_iv_data: binary instrument designs are binary where expected", {
  for (des in c("dgp1_iv_justid_linear", "dgp2_iv_homog_overid",
                "dgp3_iv_heterog_overid", "dgp4_iv_weak_local")) {
    dgp <- sim_iv_data(des, n = 120L, kx = 2L, q = 2L, seed = 3L)
    for (z in paste0("Z", seq_len(dgp$q))) {
      expect_true(all(dgp$data[[z]] %in% c(0, 1)), info = paste(des, z))
    }
  }
})

test_that("sim_iv_data: mixed design has binary first instrument and continuous second", {
  dgp <- sim_iv_data("dgp6_iv_nonlinear_firststage", n = 200L, kx = 2L, q = 2L, seed = 4L)
  expect_true(all(dgp$data$Z1 %in% c(0, 1)))
  expect_gt(length(unique(dgp$data$Z2)), 10L)
})

test_that("sim_iv_data: correlated instrument design creates strong correlation", {
  dgp <- sim_iv_data(
    "dgp5_iv_correlated_instruments",
    n = 1500L, kx = 2L, q = 2L, seed = 5L, rho = 0.85
  )
  r <- stats::cor(dgp$data$Z1 - dgp$data$piZ1_oracle,
                  dgp$data$Z2 - dgp$data$piZ2_oracle)
  expect_gt(r, 0.75)
  expect_lt(r, 0.95)
})

test_that("sim_iv_data: theta_x0 matches oracle function at x0", {
  dgp <- sim_iv_data("dgp3_iv_heterog_overid", n = 50L, kx = 2L, q = 2L, seed = 6L)
  theta_at_x0 <- dgp$truth$theta_fn(matrix(dgp$x0, nrow = 1L))
  expect_equal(dgp$theta_x0, as.numeric(theta_at_x0), tolerance = 1e-12)
})

test_that("sim_iv_data: truth functions have stable vectorized contracts", {
  designs <- c(
    "dgp1_iv_justid_linear",
    "dgp2_iv_homog_overid",
    "dgp3_iv_heterog_overid",
    "dgp4_iv_weak_local",
    "dgp5_iv_correlated_instruments",
    "dgp6_iv_nonlinear_firststage"
  )

  for (des in designs) {
    dgp <- sim_iv_data(des, n = 90L, kx = 3L, q = 3L, seed = 14L)
    xx <- rbind(dgp$x0, dgp$x0 + c(0.1, -0.1, 0.2))

    pi_x <- dgp$truth$pi_fn(xx)
    muD_x <- dgp$truth$muD_fn(xx)
    muY_x <- dgp$truth$muY_fn(xx)
    theta_x <- dgp$truth$theta_fn(xx)

    expect_equal(dim(pi_x), c(2L, dgp$q), info = des)
    expect_equal(length(muD_x), 2L, info = des)
    expect_equal(length(muY_x), 2L, info = des)
    expect_equal(length(theta_x), 2L, info = des)
    expect_true(all(is.finite(pi_x)), info = des)
    expect_true(all(is.finite(muD_x)), info = des)
    expect_true(all(muD_x >= 0 & muD_x <= 1), info = des)
    expect_true(all(is.finite(muY_x)), info = des)
    expect_true(all(is.finite(theta_x)), info = des)
    expect_equal(length(dgp$truth$instrument_type), dgp$q, info = des)
  }
})

test_that("sim_iv_data: treatment and oracle nuisance columns stay in range", {
  dgp <- sim_iv_data("dgp6_iv_nonlinear_firststage", n = 200L, kx = 2L,
                     q = 3L, seed = 15L)
  pi_cols <- paste0("piZ", seq_len(dgp$q), "_oracle")

  expect_true(all(dgp$data$D %in% c(0L, 1L)))
  expect_true(all(dgp$data$piZ1_oracle > 0 & dgp$data$piZ1_oracle < 1))
  expect_true(all(is.finite(as.matrix(dgp$data[, pi_cols]))))
  expect_true(all(dgp$data$muD_oracle >= 0 & dgp$data$muD_oracle <= 1))
  expect_true(all(is.finite(dgp$data$muY_oracle)))
  expect_true(all(is.finite(dgp$data$Y)))
})

test_that("sim_iv_data: validation is explicit", {
  expect_error(sim_iv_data("bad_design", n = 10L, kx = 2L), "Unknown IV design")
  expect_error(sim_iv_data("dgp2_iv_homog_overid", n = 10L, kx = 1L), "kx")
  expect_error(sim_iv_data("dgp2_iv_homog_overid", n = 10L, kx = 2L, q = 4L), "q")
  expect_error(sim_iv_data("dgp2_iv_homog_overid", n = 10L, kx = 2L, sigma = -1), "sigma")
  expect_error(sim_iv_data("dgp2_iv_homog_overid", n = 10L, kx = 2L, clip = 0.5), "clip")
  expect_error(
    sim_iv_data("dgp2_iv_homog_overid", n = 10L, kx = 2L, x0 = c(0, 0, 0)),
    "x0 length"
  )
  expect_error(
    sim_iv_data("dgp2_iv_homog_overid", n = 10L, kx = 2L, x0 = c(0, NA_real_)),
    "x0.*NA"
  )
  expect_error(
    sim_iv_data("dgp2_iv_homog_overid", n = 10L, kx = 2L, x0 = c(0, Inf)),
    "x0.*non-finite"
  )
})
