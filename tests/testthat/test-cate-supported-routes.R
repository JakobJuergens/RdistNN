test_that("CATE exposes only the supported nuisance construction", {
  expect_false("m_col" %in% names(formals(cate_dnn)))
  expect_false(exists("cate_rho_rlearner", envir = asNamespace("RdistNN"), inherits = FALSE))
  d <- cate_demo_data(n = 20, seed = 9)$data
  expect_false("rho_oracle" %in% names(d))
  expect_true(all(c("m0_oracle", "m1_oracle", "aipw_oracle") %in% names(d)))
  expect_error(cate_dnn(data = d, x_cols = c("X1", "X2"),
    y_col = "Y", t_col = "D", e_col = "e_oracle", m_col = "m_oracle",
    x0 = c(0, 0), s = 3), "unused argument.*m_col")
  sim <- sim_cate_data("dgp1_linear", n = 20, seed = 9)$data
  expect_false("rho_oracle" %in% names(sim))
})
