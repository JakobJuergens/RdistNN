test_that("compiled CATE batch matches R batch across localizers", {
  skip_if_not(has_cate_cpp(), "C++ CATE functions not available")

  demo <- cate_demo_data(n = 120, p = 2, design = "hetero", seed = 5)
  df <- demo$data
  x0 <- demo$x0
  rho <- df$rho_oracle

  for (loc in c("dnn1","dnnk","kernel")) {
    sc  <- score_cate_rho(x0 = x0, x_cols = c("X1","X2"), rho = rho, localizer = loc, k = 5)
    ctx <- score_setup(sc, df)     # ctx$ptr_ctx will be non-null if compiled available

    set.seed(99)
    subs <- cate_demo_subsets(n = nrow(df), K = 200, s = 8)

    # Compiled path (score_eval_batch will delegate to rustat_cate_eval_batch)
    vals_cpp <- score_eval_batch(sc, subs, df, ctx)

    # Force R fallback by removing ptr from ctx
    ctx_r <- ctx
    ctx_r$ptr_ctx <- NULL
    vals_r <- score_eval_batch(sc, subs, df, ctx_r)

    # Allow tiny FP differences for kernel bandwidth auto-selection
    expect_equal(vals_cpp, vals_r, tolerance = 1e-8, ignore_attr = TRUE)
  }
})

test_that("compiled CATE returns finite values and respects subset size", {
  skip_if_not(has_cate_cpp(), "C++ CATE functions not available")

  demo <- cate_demo_data(n = 90, p = 2, seed = 10)
  df <- demo$data
  x0 <- demo$x0
  rho <- df$rho_oracle

  sc  <- score_cate_rho(x0 = x0, x_cols = c("X1","X2"), rho = rho, localizer = "dnnk", k = 4)
  ctx <- score_setup(sc, df)

  for (s in c(3L, 5L, 9L)) {
    set.seed(123 + s)
    subs <- cate_demo_subsets(n = nrow(df), K = 50, s = s)
    vals <- score_eval_batch(sc, subs, df, ctx)  # compiled branch
    expect_true(all(is.finite(vals)))
    expect_length(vals, nrow(subs))
  }
})
