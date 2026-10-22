test_that("native indexed scores reject represented invalid indices", {
  old <- options(RdistNN.use_compiled = TRUE, RdistNN.num_threads = 2L)
  on.exit(options(old), add = TRUE)
  native <- function(name) {
    out <- RdistNN:::.rdistnn_compiled_fun(name)
    expect_true(is.function(out))
    out
  }
  cate_create <- native("rustat_cate_ctx_create")
  plqr_create <- native("rustat_plqr_ctx_create")
  cases <- list(
    list(eval = native("rustat_cate_g_indexed"), ctx = cate_create(c(.5, -1, 2))),
    list(eval = native("rustat_plqr_g_indexed"), ctx = plqr_create(
      c(.5, 1.5, 2.5), matrix(c(1, 2, 3)), rep(0, 3L),
      matrix(0, 3L, 1L), .5, .5, 1L))
  )
  for (case in cases) for (size in c(4L, 1100L)) {
    for (bad in c(-1L, NA_integer_, 0L, 4L, .Machine$integer.max)) {
      indices <- rep(c(3L, 1L, 3L, 2L), length.out = size)
      indices[size] <- bad
      expect_error(case$eval(indices, case$ctx), "Index out of bounds")
    }
  }
})

test_that("native indexed scores preserve repeated and reordered requests", {
  old <- options(RdistNN.use_compiled = TRUE, RdistNN.num_threads = 1L)
  on.exit(options(old), add = TRUE)
  cate <- RdistNN:::rustat_cate_ctx_create(c(.5, -1, 2))
  plqr <- RdistNN:::rustat_plqr_ctx_create(
    c(.5, 1.5, 2.5), matrix(c(1, 2, 3)), rep(0, 3L),
    matrix(0, 3L, 1L), .5, .5, 1L)
  for (threads in c(1L, 2L)) {
    options(RdistNN.num_threads = threads)
    for (size in c(4L, 1100L)) {
      indices <- rep(c(3L, 1L, 3L, 2L), length.out = size)
      expect_identical(RdistNN:::rustat_cate_g_indexed(indices, cate),
        c(.5, -1, 2)[indices])
      expect_identical(RdistNN:::rustat_plqr_g_indexed(indices, plqr),
        c(-.5, 1, 1.5)[indices])
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
