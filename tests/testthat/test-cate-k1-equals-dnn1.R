test_that("dnnk with k=1 equals dnn1 exactly", {
  demo <- cate_demo_data(n = 80, p = 2, seed = 2)
  df <- demo$data; x0 <- demo$x0; rho <- df$rho_oracle

  sc1 <- score_cate_rho(x0, c("X1","X2"), rho = rho, localizer = "dnn1")
  sck <- score_cate_rho(x0, c("X1","X2"), rho = rho, localizer = "dnnk", k = 1)

  ctx1 <- score_setup(sc1, df); ctxk <- score_setup(sck, df)
  ctx1$ptr_ctx <- NULL; ctxk$ptr_ctx <- NULL  # pure-R to be strict

  subs <- cate_demo_subsets(nrow(df), K = 40, s = 7)
  expect_equal(score_eval_batch(sc1, subs, df, ctx1),
               score_eval_batch(sck, subs, df, ctxk),
               tolerance = 0)  # exact
})
