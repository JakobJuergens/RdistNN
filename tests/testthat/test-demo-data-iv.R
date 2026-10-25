test_that("iv_demo_data: returns expected schema", {
  demo <- iv_demo_data(n = 80L, p = 2L, q = 2L, design = "homog", seed = 1L)
  expect_type(demo, "list")
  expect_named(demo, c("data", "x0", "q", "theta_x0", "truth"))
  expect_s3_class(demo$data, "data.frame")
  expect_equal(nrow(demo$data), 80L)
  expect_equal(demo$q, 2L)
  expect_true(all(c(
    "X1", "X2", "Y", "D", "Z1", "Z2",
    "piZ1_oracle", "piZ2_oracle",
    "muY_oracle", "muD_oracle", "theta_oracle"
  ) %in% names(demo$data)))
})

test_that("iv_demo_data: deterministic with fixed seed", {
  d1 <- iv_demo_data(n = 80L, p = 2L, q = 2L, design = "heterog", seed = 2L)
  d2 <- iv_demo_data(n = 80L, p = 2L, q = 2L, design = "heterog", seed = 2L)
  expect_equal(d1$data, d2$data)
  expect_equal(d1$theta_x0, d2$theta_x0)
})

test_that("iv_demo_data: injection hooks are honored", {
  d_na <- iv_demo_data(
    n = 40L, p = 2L, q = 2L, design = "homog", seed = 3L,
    inject = "na", inject_col = "Y", idx_inject = 5L
  )
  expect_true(is.na(d_na$data$Y[5L]))

  d_inf <- iv_demo_data(
    n = 40L, p = 2L, q = 2L, design = "homog", seed = 3L,
    inject = "inf", inject_col = "Z1", idx_inject = 4L
  )
  expect_true(is.infinite(d_inf$data$Z1[4L]))

  d_pi <- iv_demo_data(
    n = 40L, p = 2L, q = 2L, design = "homog", seed = 3L,
    inject = "na", inject_col = "piZ2_oracle", idx_inject = 6L
  )
  expect_true(is.na(d_pi$data$piZ2_oracle[6L]))
})

test_that("iv_demo_data: aliases and x0 are forwarded to IV DGPs", {
  x0 <- c(0.25, -0.5, 0.1)
  demo <- iv_demo_data(n = 60L, p = 3L, q = 3L, design = "weak",
                       seed = 4L, x0 = x0)
  expect_equal(demo$x0, x0)
  expect_equal(demo$q, 3L)
  expect_equal(demo$truth$instrument_type, rep("binary", 3L))

  expect_error(iv_demo_data(n = 60L, p = 2L, q = 2L, design = "bad"),
               class = "rdistnn_input_error")
})

test_that("iv_demo_data: rejects invalid injection requests", {
  expect_error(
    iv_demo_data(n = 20L, p = 2L, q = 2L, inject = "bad"),
    "inject"
  )
  expect_error(
    iv_demo_data(n = 20L, p = 2L, q = 2L, inject = "na", inject_col = "bad_col"),
    "inject_col"
  )
  expect_error(
    iv_demo_data(n = 20L, p = 2L, q = 2L, inject = "na", idx_inject = 30L),
    "idx_inject"
  )
})
