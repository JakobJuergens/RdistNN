test_that("score_cate_rho enforces no internal fitting (rho XOR rho_col)", {
  demo <- cate_demo_data(n = 50, p = 2, seed = 123)
  df <- demo$data
  x0 <- demo$x0
  rho <- df$rho_oracle

  # neither provided -> error
  expect_error(
    score_cate_rho(x0 = x0, x_cols = c("X1","X2")),
    regexp = "exactly one|mutually exclusive",
    ignore.case = TRUE
  )

  # both provided -> error
  expect_error(
    score_cate_rho(x0 = x0, x_cols = c("X1","X2"), rho_col = "rho", rho = df$rho_oracle),
    regexp = "exactly one|mutually exclusive",
    ignore.case = TRUE
  )

  # only rho -> ok
  sc1 <- score_cate_rho(x0 = x0, x_cols = c("X1","X2"), rho = rho)
  expect_true(is_score(sc1))

  # only rho_col -> ok
  sc2 <- score_cate_rho(x0 = x0, x_cols = c("X1","X2"), rho_col = "rho")
  expect_true(is_score(sc2))
})

test_that("score_cate_rho setup builds ctx and returns finite batch outputs", {
  demo <- cate_demo_data(n = 80, p = 2, seed = 1)
  df <- demo$data
  x0 <- demo$x0
  rho <- df$rho_oracle

  for (loc in c("dnn1","dnnk","kernel")) {
    sc <- score_cate_rho(x0 = x0, x_cols = c("X1","X2"), rho = rho, localizer = loc, k = 5)
    ctx <- score_setup(sc, df)
    s   <- score_order(sc, nrow(df), df, ctx)
    expect_type(s, "integer")
    expect_true(s >= 1L)

    set.seed(42)
    subs <- cate_demo_subsets(n = nrow(df), K = 25, s = 8)
    vals <- score_eval_batch(sc, subs, df, ctx)

    expect_length(vals, nrow(subs))
    expect_true(all(is.finite(vals)))
  }
})

test_that("single-subset eval equals batch entry for same subset", {
  demo <- cate_demo_data(n = 60, p = 2, seed = 2)
  df <- demo$data
  x0 <- demo$x0
  rho <- df$rho_oracle

  sc  <- score_cate_rho(x0 = x0, x_cols = c("X1","X2"), rho = rho, localizer = "dnnk", k = 3)
  ctx <- score_setup(sc, df)
  # Force R fallback (no compiled) to test S3 eval vs eval_batch
  ctx$ptr_ctx <- NULL

  set.seed(7)
  subs <- cate_demo_subsets(n = nrow(df), K = 10, s = 7)
  vals_b <- score_eval_batch(sc, subs, df, ctx)

  vals_s <- vapply(
    seq_len(nrow(subs)),
    function(i) score_eval(sc, subs[i, ], df, ctx),
    numeric(1)
  )
  expect_equal(vals_b, vals_s, tolerance = 1e-12)
})

test_that("kernel is permutation-invariant within a subset", {
  demo <- cate_demo_data(n = 75, p = 2, seed = 3)
  df <- demo$data
  x0 <- demo$x0
  rho <- df$rho_oracle

  sc  <- score_cate_rho(x0 = x0, x_cols = c("X1","X2"), rho = rho, localizer = "kernel", h = NULL)
  ctx <- score_setup(sc, df)
  ctx$ptr_ctx <- NULL  # ensure we test the R path deterministically

  set.seed(11)
  subs <- cate_demo_subsets(n = nrow(df), K = 8, s = 6)
  for (i in seq_len(nrow(subs))) {
    idxs <- subs[i, ]
    perm <- sample(idxs)
    v1 <- score_eval(sc, idxs, df, ctx)
    v2 <- score_eval(sc, perm, df, ctx)
    expect_equal(v1, v2, tolerance = 1e-12)
  }
})

test_that("rho must be finite", {
  demo <- cate_demo_data(n = 20, p = 2, seed = 1)
  df <- demo$data
  bad_rho <- demo$data$rho_oracle
  bad_rho[3] <- NA_real_

  sc <- score_cate_rho(x0 = demo$x0, x_cols = c("X1","X2"), rho = bad_rho)
  ctx <- score_setup(sc, df); ctx$ptr_ctx <- NULL

  subs <- matrix(c(1,2,3,4,5,6,7,8), nrow = 1)
  expect_error(score_eval_batch(sc, subs, df, ctx), "non-finite", ignore.case = TRUE)
})

