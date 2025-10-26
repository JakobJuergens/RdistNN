test_that("ties in distances produce deterministic tie-breaking (lower index first)", {
  # Craft X so rows 1 and 2 are exactly equidistant to x0, others far
  n <- 10
  df <- data.frame(
    X1 = c(0, 0, rep(10, n-2)),
    X2 = c(1, -1, rep(10, n-2)),
    Y = 0, D = 0
  )
  rho <- seq_len(n)   # rho equals index -> reveals which indices are chosen
  x0  <- c(0, 0)

  sc  <- score_cate_rho(x0, c("X1","X2"), rho = rho, localizer = "dnnk", k = 2)
  ctx <- score_setup(sc, df)

  subs <- rbind(c(1,2,3,4,5,6), c(2,1,3,4,5,6)) # both contain the tied pair
  v_cpp <- score_eval_batch(sc, subs, df, ctx)   # compiled if available
  ctx$ptr_ctx <- NULL
  v_r   <- score_eval_batch(sc, subs, df, ctx)   # R fallback

  # With ties, indices (1,2) should be selected (lower indices), so mean = (1+2)/2 = 1.5
  expect_equal(v_cpp[1], 1.5, tolerance = 1e-12)
  expect_equal(v_cpp, v_r, tolerance = 1e-12)  # identical paths
})
