test_that("localized value lies within convex hull of subset rhos", {
  demo <- cate_demo_data(n = 60, p = 2, seed = 3)
  df <- demo$data; x0 <- demo$x0; rho <- df$rho_oracle

  for (loc in c("dnn1","dnnk","kernel")) {
    sc  <- score_cate_rho(x0, c("X1","X2"), rho = rho, localizer = loc, k = 3, h = NULL)
    ctx <- score_setup(sc, df); ctx$ptr_ctx <- NULL  # test R path

    set.seed(10)
    subs <- cate_demo_subsets(nrow(df), K = 25, s = 6)
    vals <- score_eval_batch(sc, subs, df, ctx)

    mins <- apply(subs, 1L, function(id) min(rho[id]))
    maxs <- apply(subs, 1L, function(id) max(rho[id]))
    expect_true(all(vals >= mins - 1e-12))
    expect_true(all(vals <= maxs + 1e-12))
  }
})
