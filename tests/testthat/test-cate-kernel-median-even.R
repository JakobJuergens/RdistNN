test_that("kernel auto-bandwidth uses true median for even s (matches R::median)", {
  demo <- cate_demo_data(n = 50, p = 2, seed = 5)
  df <- demo$data; x0 <- demo$x0; rho <- df$rho_oracle

  sc  <- score_cate_rho(x0, c("X1","X2"), rho = rho, localizer = "kernel", h = NULL)
  ctx <- score_setup(sc, df)

  # Use even s to exercise median averaging
  set.seed(12)
  subs <- cate_demo_subsets(nrow(df), K = 60, s = 8)

  v_cpp <- score_eval_batch(sc, subs, df, ctx)   # compiled if present
  ctx$ptr_ctx <- NULL
  v_r   <- score_eval_batch(sc, subs, df, ctx)   # R fallback
  expect_equal(v_cpp, v_r, tolerance = 1e-8)
})
