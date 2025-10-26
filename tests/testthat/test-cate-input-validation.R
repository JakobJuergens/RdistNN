test_that("score_cate_rho validates inputs early and in setup", {
  demo <- cate_demo_data(n = 30, p = 3, seed = 1)
  df <- demo$data
  rho <- df$rho_oracle

  # XOR at factory
  expect_error(score_cate_rho(x0 = rep(0, 3), x_cols = 1:3), "exactly one", ignore.case = TRUE)
  expect_error(score_cate_rho(x0 = rep(0, 3), x_cols = 1:3, rho = rho, rho_col = "rho"),
               "exactly one", ignore.case = TRUE)

  # x0 length mismatch
  expect_error(score_cate_rho(x0 = 0, x_cols = 1:3, rho = rho) |> score_setup(df),
               regexp = "length|ncol", ignore.case = TRUE)

  # rho length mismatch
  expect_error(score_cate_rho(x0 = rep(0, 3), x_cols = 1:3, rho = rho[-1]) |> score_setup(df),
               regexp = "length", ignore.case = TRUE)

  # rho via column vs direct are identical
  df$rho <- rho
  sc_vec <- score_cate_rho(x0 = rep(0, 3), x_cols = 1:3, rho = rho)
  sc_col <- score_cate_rho(x0 = rep(0, 3), x_cols = 1:3, rho_col = "rho")

  ctx_vec <- score_setup(sc_vec, df); ctx_col <- score_setup(sc_col, df)
  subs <- cate_demo_subsets(nrow(df), K = 25, s = 6)
  ctx_vec$ptr_ctx <- NULL; ctx_col$ptr_ctx <- NULL  # force R path
  expect_equal(score_eval_batch(sc_vec, subs, df, ctx_vec),
               score_eval_batch(sc_col, subs, df, ctx_col),
               tolerance = 1e-12)
})
