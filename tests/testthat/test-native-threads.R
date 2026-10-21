test_that("native teams obey package limits without changing the runtime default", {
  old <- options(RdistNN.num_threads = NULL, RdistNN.use_compiled = TRUE)
  on.exit(options(old), add = TRUE)
  info <- .rdistnn_compiled_fun("rustat_thread_info")
  expect_type(info, "closure")
  if (is.null(info)) return(invisible(NULL))
  default <- info()
  expect_lte(default$team_size, 2L)
  expect_lte(default$limit, 2L)
  for (threads in c(1L, 2L)) {
    options(RdistNN.num_threads = threads)
    result <- info()
    expect_identical(result$runtime_max, default$runtime_max)
    expect_identical(result$runtime_limit, default$runtime_limit)
    expect_identical(result$limit, min(threads, result$runtime_max, result$runtime_limit))
    expect_lte(result$team_size, result$limit)
    expect_gte(result$team_size, 1L)
  }
})

test_that("invalid native thread options fail before a parallel calculation", {
  old <- options(RdistNN.num_threads = NULL, RdistNN.use_compiled = TRUE)
  on.exit(options(old), add = TRUE)
  for (bad in list(0, -1, 1.5, NA_integer_, Inf, TRUE, "2", c(1, 2), numeric(), 2^31)) {
    options(RdistNN.num_threads = bad)
    expect_error(rustat_thread_info(), "RdistNN.num_threads")
  }
  options(RdistNN.num_threads = 1L)
  if (rustat_thread_info()$openmp) {
    options(RdistNN.num_threads = 0L)
    expect_error(dist_l2(matrix(0, 1024L, 2L), c(0, 0)), "RdistNN.num_threads")
  }
})

test_that("native parallel kernels retain their one-thread results", {
  old <- options(RdistNN.num_threads = 1L, RdistNN.use_compiled = TRUE)
  on.exit(options(old), add = TRUE)
  required <- c("rustat_dist_l1", "rustat_dist_l2", "rustat_cate_g_indexed",
    "rustat_plqr_g_all", "rustat_plqr_g_indexed", "rustat_subset_row_sums_g",
    "rustat_jackknife_fast_rank")
  for (name in required) expect_type(.rdistnn_compiled_fun(name), "closure")
  n <- 1100L
  X <- cbind(seq_len(n) / n, sin(seq_len(n)))
  X0 <- X[seq_len(8L), , drop = FALSE]
  rho <- cos(seq_len(n))
  cate <- rustat_cate_ctx_create(rho)
  plqr <- rustat_plqr_ctx_create(rho, X, rep(0, n), X / 2, c(0.2, -0.1), 0.5, 1L)
  ranks <- cbind(seq_len(n), c(seq.int(2L, n), 1L))
  calculate <- function() {
    list(
      l1 = dist_l1(X, X0[1L, ]), l2 = dist_l2(X, X0[1L, ]),
      l1_batch = distance_l1()$eval_batch(X, X0),
      l2_batch = distance_l2()$eval_batch(X, X0),
      cate = rustat_cate_g_indexed(seq_len(n), cate),
      plqr = rustat_plqr_g_all(plqr),
      plqr_indexed = rustat_plqr_g_indexed(seq_len(n), plqr),
      subset = rustat_subset_row_sums_g(ranks, rho),
      jackknife = rustat_jackknife_fast_rank(ranks, 1, matrix(rho, ncol = 1L))
    )
  }
  one <- calculate()
  options(RdistNN.num_threads = 2L)
  expect_identical(calculate(), one)
})
