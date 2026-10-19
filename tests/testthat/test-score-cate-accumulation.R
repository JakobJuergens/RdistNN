test_that("pure CATE scores retain cancellation in scalar and batch engines", {
  old <- options(RdistNN.use_compiled = TRUE)
  on.exit(options(old), add = TRUE)
  data <- data.frame(x = 1:3)
  sc <- score_cate_rho(0, "x", 3L, rho = c(2^54, 1, -2^54))
  ctx <- score_setup(sc, data)
  subsets <- rbind(1:3, 3:1, c(1L, 3L, 2L))
  for (native in c(TRUE, FALSE)) {
    options(RdistNN.use_compiled = native)
    if (native) expect_false(is.null(.rdistnn_compiled_fun("rustat_subset_row_sums_g")))
    expect_identical(score_eval_batch(sc, subsets, data, ctx), rep(1, 3))
    for (i in 1:3) expect_identical(score_eval(sc, subsets[i, ], data, ctx), 1)
  }
  options(RdistNN.use_compiled = TRUE)
  cpp <- .rdistnn_compiled_fun("rustat_subset_row_sums_g")
  for (index in list(0L, -1L, 4L, NA_integer_)) {
    expect_error(cpp(matrix(index), ctx$g), "Index out of bounds")
  }
})
