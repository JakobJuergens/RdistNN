test_that("compiled dispatch can be disabled explicitly", {
  skip_if(
    is.null(.rdistnn_compiled_fun("rustat_dist_l2")),
    "Compiled RdistNN symbols are unavailable."
  )
  old <- getOption("RdistNN.use_compiled")
  on.exit(options(RdistNN.use_compiled = old), add = TRUE)

  X <- matrix(c(
    0, 0,
    1, -1,
    2, 0.5,
    -0.5, 2
  ), ncol = 2, byrow = TRUE)
  X0 <- rbind(c(0.25, -0.5), c(1.5, 1))

  options(RdistNN.use_compiled = TRUE)
  compiled <- list(
    l1 = dist_l1(X, X0[1, ]),
    l2 = dist_l2(X, X0[1, ]),
    l1_batch = distance_l1()$eval_batch(X, X0),
    l2_batch = distance_l2()$eval_batch(X, X0),
    l1_order = order_by_l1(X, X0[1, ]),
    l2_order = order_by_l2(X, X0[1, ])
  )

  options(RdistNN.use_compiled = FALSE)
  expect_false(.rdistnn_compiled_enabled())
  expect_null(.rdistnn_compiled_fun("rustat_dist_l2"))
  fallback <- list(
    l1 = dist_l1(X, X0[1, ]),
    l2 = dist_l2(X, X0[1, ]),
    l1_batch = distance_l1()$eval_batch(X, X0),
    l2_batch = distance_l2()$eval_batch(X, X0),
    l1_order = order_by_l1(X, X0[1, ]),
    l2_order = order_by_l2(X, X0[1, ])
  )

  expect_equal(fallback, compiled, tolerance = 1e-12)
})

test_that("rank jackknife compiled and R engines agree under the option", {
  withr::local_preserve_seed()
  skip_if(
    is.null(.rdistnn_compiled_fun("rustat_jackknife_fast_rank")),
    "Compiled jackknife symbol is unavailable."
  )
  old <- getOption("RdistNN.use_compiled")
  on.exit(options(RdistNN.use_compiled = old), add = TRUE)

  set.seed(20260712)
  n <- 48L
  s <- 6L
  X <- matrix(rnorm(n * 3L), n, 3L)
  ranks <- knn_ranks(X, k = s + 1L, exclude_self = TRUE)
  PHI <- cbind(rnorm(n), rnorm(n))
  a <- seq_len(s) / sum(seq_len(s))

  options(RdistNN.use_compiled = TRUE)
  compiled <- jackknife_fast_rank(ranks, a, PHI, normalize = TRUE)
  options(RdistNN.use_compiled = FALSE)
  fallback <- jackknife_fast_rank(ranks, a, PHI, normalize = TRUE)

  expect_equal(fallback$theta0, compiled$theta0, tolerance = 1e-12)
  expect_equal(fallback$theta_loo, compiled$theta_loo, tolerance = 1e-12)
})

test_that("np_dnn estimates and inference agree across engine modes", {
  withr::local_preserve_seed()
  skip_if(
    is.null(.rdistnn_compiled_fun("rustat_dist_l2")),
    "Compiled RdistNN symbols are unavailable."
  )
  old <- getOption("RdistNN.use_compiled")
  on.exit(options(RdistNN.use_compiled = old), add = TRUE)

  set.seed(20260713)
  X <- matrix(rnorm(120L), 60L, 2L)
  Y <- 1 + X[, 1L] - 0.5 * X[, 2L] + rnorm(60L, sd = 0.2)

  options(RdistNN.use_compiled = TRUE)
  compiled <- np_dnn(
    X = X, x0 = c(0.2, -0.1), Y = Y,
    k = 1L, s = 12L, se = TRUE, jackknife = "fast"
  )
  options(RdistNN.use_compiled = FALSE)
  fallback <- np_dnn(
    X = X, x0 = c(0.2, -0.1), Y = Y,
    k = 1L, s = 12L, se = TRUE, jackknife = "fast"
  )

  expect_equal(fallback$mu, compiled$mu, tolerance = 1e-12)
  expect_equal(fallback$se, compiled$se, tolerance = 1e-12)
  expect_equal(fallback$meta$s, compiled$meta$s)
})

test_that("PLQR score contributions agree across compiled and R engines", {
  skip_if(
    is.null(.rdistnn_compiled_fun("rustat_plqr_ctx_create")) ||
      is.null(.rdistnn_compiled_fun("rustat_subset_row_sums_g")),
    "Compiled PLQR score symbols are unavailable."
  )
  old <- getOption("RdistNN.use_compiled")
  on.exit(options(RdistNN.use_compiled = old), add = TRUE)

  demo <- plqr_demo_data(
    n = 36L, p = 2L, p_d = 2L, tau = 0.4, seed = 75L
  )
  score <- score_plqr_component(
    x0 = demo$x0,
    x_cols = c("X1", "X2"),
    D_cols = c("D1", "D2"),
    g_col = "g_oracle",
    mD_cols = c("mD1", "mD2"),
    y_col = "Y",
    theta = demo$truth$theta_true,
    j = 2L,
    tau = 0.4,
    s = 3L
  )
  subsets <- rbind(
    c(1L, 2L, 3L),
    c(4L, 8L, 12L),
    c(20L, 21L, 36L)
  )

  options(RdistNN.use_compiled = TRUE)
  compiled_ctx <- score_setup(score, demo$data)
  compiled <- list(
    contribution = compiled_ctx$score,
    single = score_eval(score, subsets[1L, ], demo$data, compiled_ctx),
    batch = score_eval_batch(score, subsets, demo$data, compiled_ctx)
  )

  options(RdistNN.use_compiled = FALSE)
  fallback_ctx <- score_setup(score, demo$data)
  fallback <- list(
    contribution = fallback_ctx$score,
    single = score_eval(score, subsets[1L, ], demo$data, fallback_ctx),
    batch = score_eval_batch(score, subsets, demo$data, fallback_ctx)
  )

  expect_equal(fallback, compiled, tolerance = 1e-12)
})
