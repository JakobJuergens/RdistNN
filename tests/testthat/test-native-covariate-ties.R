test_that("native metric orders use input index for exact distance ties", {
  old <- options(RdistNN.use_compiled = TRUE)
  on.exit(options(old), add = TRUE)
  X <- rbind(c(1, 0), c(-1, 0), c(0, 1), c(0, -1), c(-1, 0), c(2, 0))
  subsets <- list(1:6, c(6L, 4L, 2L, 1L, 5L, 3L), c(5L, 4L, 1L, 6L))
  for (name in c("rustat_order_by_l1", "rustat_order_by_l1_scaled",
                 "rustat_order_by_l2", "rustat_order_by_l2_scaled")) {
    native <- .rdistnn_compiled_fun(name)
    expect_true(is.function(native))
    for (at in subsets) {
      args <- list(X = X[at, , drop = FALSE], x0 = c(0, 0))
      if (grepl("scaled", name, fixed = TRUE)) args$inv_s <- c(1, 1)
      actual <- do.call(native, args)
      expect_identical(actual$ord, order(ifelse(at == 6L, 2, 1), seq_along(at)))
      expect_identical(as.numeric(actual$dx), ifelse(at == 6L, 2, 1))
    }
    # A represented one-ULP distance difference precedes input index.
    args <- list(X = matrix(c(-1 - 2^-52, 1)), x0 = 0)
    if (grepl("scaled", name, fixed = TRUE)) args$inv_s <- 1
    expect_identical(do.call(native, args)$ord, c(2L, 1L))
  }
})

test_that("distance-only native order and top-k retain their index tie contract", {
  old <- options(RdistNN.use_compiled = TRUE)
  on.exit(options(old), add = TRUE)
  expect_identical(.rdistnn_compiled_fun("rustat_argsort_dist")(c(1, 1, 0)), c(3L, 1L, 2L))
  expect_identical(.rdistnn_compiled_fun("rustat_ranks_from_dist")(c(1, 1, 0)), c(2L, 3L, 1L))
  expect_identical(.rdistnn_compiled_fun("rustat_topk_dist")(c(1, 1, 0), 2L), c(3L, 1L))
})
