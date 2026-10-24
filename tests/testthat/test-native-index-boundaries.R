test_that("native production subset scores reject invalid indices at serial and parallel sizes", {
  old <- options(RdistNN.use_compiled = TRUE, RdistNN.num_threads = 2L)
  on.exit(options(old), add = TRUE)
  native <- RdistNN:::.rdistnn_compiled_fun("rustat_subset_row_sums_g")
  expect_true(is.function(native))
  for (size in c(4L, 1100L)) for (column in 1:3) {
    for (bad in c(-1L, NA_integer_, 0L, 4L, .Machine$integer.max)) {
      subsets <- matrix(rep(c(3L, 1L, 3L), each = size), size, 3L)
      subsets[size, column] <- bad
      expect_error(native(subsets, c(.5, -1, 2)), "Index out of bounds")
    }
  }
})

test_that("native production subset scores preserve repeated and reordered requests", {
  old <- options(RdistNN.use_compiled = TRUE, RdistNN.num_threads = 1L)
  on.exit(options(old), add = TRUE)
  native <- RdistNN:::.rdistnn_compiled_fun("rustat_subset_row_sums_g")
  expect_true(is.function(native))
  for (threads in c(1L, 2L)) {
    options(RdistNN.num_threads = threads)
    for (size in c(4L, 1100L)) {
      indices <- rep(c(3L, 1L, 3L, 2L), length.out = size)
      expect_identical(native(matrix(indices, ncol = 1L), c(.5, -1, 2)), c(.5, -1, 2)[indices])
      subsets <- cbind(indices, rep(3L, size), rep(1L, size))
      expect_identical(native(subsets, c(.5, -1, 2)), c(.5, -1, 2)[indices] + 2.5)
    }
  }
})

test_that("native finite-support ranks validate before converting indices", {
  old <- options(RdistNN.use_compiled = TRUE)
  on.exit(options(old), add = TRUE)
  native <- RdistNN:::.rdistnn_compiled_fun("rustat_jackknife_fast_rank")
  expect_true(is.function(native))
  signal <- cbind(c(1, -3, 5), c(2, 4, -2))
  for (position in 1:2) for (bad in c(NA_integer_, 0L, -1L, 4L, .Machine$integer.max)) {
    ranks <- matrix(c(1L, 2L), nrow = 1L)
    ranks[position] <- bad
    expect_error(native(ranks, 1, signal), "ranks indices out of range")
  }
  ranks <- rbind(c(3L, 1L, 3L), c(1L, 2L, 3L), c(2L, 1L, 2L))
  actual <- native(ranks, c(1, -.5), signal)
  expected <- jackknife_slow_rank(ranks, c(1, -.5), signal)
  expect_identical(actual$theta0, expected$theta0)
  expect_identical(actual$theta_loo, expected$theta_loo)
})

test_that("native full-rank deletions validate every permutation before indexing", {
  old <- options(RdistNN.use_compiled = TRUE)
  on.exit(options(old), add = TRUE)
  native <- RdistNN:::.rdistnn_compiled_fun("rustat_jackknife_rank_full")
  expect_true(is.function(native))
  signal <- cbind(c(1, -3, 5), c(2, 4, -2))
  full <- c(.5, .25, .25)
  table <- c(.75, .25)
  for (target in 1:2) for (bad in c(NA_integer_, 0L, -1L, 4L, .Machine$integer.max, 2L)) {
    ranks <- cbind(1:3, 1:3)
    ranks[1L, target] <- bad
    expect_error(native(ranks, signal, full, table), "permutation")
  }
  ranks <- cbind(c(3L, 1L, 2L), c(2L, 3L, 1L))
  actual <- native(ranks, signal, full, table)
  baseline <- matrix(0, 2L, 2L)
  deleted <- array(0, c(3L, 2L, 2L))
  for (target in 1:2) for (component in 1:2) {
    baseline[component, target] <- sum(full[ranks[, target]] * signal[, component])
    for (i in 1:3) {
      keep <- (1:3) != i
      reduced <- ranks[keep, target] - (ranks[keep, target] > ranks[i, target])
      deleted[i, component, target] <- sum(table[reduced] * signal[keep, component])
    }
  }
  expect_identical(actual$theta0, baseline)
  expect_identical(actual$theta_loo, deleted)
})
