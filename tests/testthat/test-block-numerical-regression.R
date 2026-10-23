test_that("deleting a dominant duplicate retains the actual small survivor", {
  original <- options(RdistNN.use_compiled = TRUE)
  on.exit(options(original), add = TRUE)
  X <- matrix(c(0, 0, 1))
  for (native in c(FALSE, TRUE)) {
    options(RdistNN.use_compiled = native)
    expect_identical(is.null(RdistNN:::.rdistnn_compiled_fun("rustat_jackknife_rank_full")), !native)
    expect_identical(is.null(RdistNN:::.rdistnn_compiled_fun("rustat_block_rank_loo")), !native)
    for (method in c("fast", "loo")) {
      fit <- suppressWarnings(affine_ratio_dnn(X = X, x0 = 0,
        numerator = c(2^54, 1, 0), denominator = rep(1, 3L),
        s = 2L, jackknife = method, return = "deletions"))
      expect_equal(as.numeric(fit$deletions[, 1L, 1L]), c(1, 2^54, 2^53), tolerance = 0)
    }
  }
})

test_that("direct native block ranks reject missing integer values before indexing", {
  old <- options(RdistNN.use_compiled = TRUE)
  on.exit(options(old), add = TRUE)
  native <- RdistNN:::.rdistnn_compiled_fun("rustat_block_rank_loo")
  expect_true(is.function(native))
  expect_error(native(matrix(c(NA_integer_, 2L, 3L)), matrix(1:3), c(1, 0),
    c(1L, 1L, 2L)), "permutation")
  expect_error(native(matrix(c(.Machine$integer.max, 2L, 3L)), matrix(1:3), c(1, 0),
    c(1L, 1L, 2L)), "permutation")
})

test_that("signed block deletions cannot turn exactly represented variation into zero", {
  original <- options(RdistNN.use_compiled = TRUE)
  on.exit(options(original), add = TRUE)
  for (native in c(FALSE, TRUE)) {
    options(RdistNN.use_compiled = native)
    jk <- RdistNN:::.dnn_loo_centered(c(2^54, 1, 2^53, 2^53),
      matrix(1:4), c(0, 1, -1), groups = c(1L, 1L, 2L, 3L))
    expect_equal(as.numeric(jk$loo), c(0, 0, .5, .5), tolerance = 0)
    expect_equal(as.numeric(RdistNN:::.covariance_matrix(jk$loo, 3 / 4)), 3 / 16, tolerance = 0)
  }
})

test_that("actual subnormal duplicate signals follow direct rounded products", {
  original <- options(RdistNN.use_compiled = TRUE)
  on.exit(options(original), add = TRUE)
  eta <- .Machine$double.xmin * .Machine$double.eps
  expect_gt(eta, 0)
  for (native in c(FALSE, TRUE)) {
    options(RdistNN.use_compiled = native)
    for (mass in c(1, 2^1000)) {
      jk <- RdistNN:::.dnn_loo_centered(c(eta, 0, 0), matrix(1:3), c(mass, 0),
        groups = c(1L, 1L, 2L))
      actual <- RdistNN:::.dnn_restore_loo(jk$loo, jk$offset)
      expected <- c(0, eta * mass, eta * (mass / 2))
      expect_equal(as.numeric(actual) / (eta * mass), expected / (eta * mass), tolerance = 0)
      weight_fun <- function(i) {
        matrix(if (i == 3L) rep(mass / 2, 2L) else c(mass, 0), ncol = 1L)
      }
      recomputed <- RdistNN:::.jackknife_recomputed_centered(
        matrix(c(eta, 0, 0)), weight_fun, restore = TRUE)
      expect_equal(as.numeric(recomputed) / (eta * mass), expected / (eta * mass), tolerance = 0)
    }
    rank_only <- jackknife_rank_full(matrix(1:3), matrix(c(eta, 0, 0)),
      c(1, 0, 0), c(.5, .5))$theta_loo
    expect_equal(as.numeric(rank_only), rep(0, 3L), tolerance = 0)
  }
})

test_that("interleaved group masses retain prefix and suffix cancellation residuals", {
  table <- c(2^54, 0, 1, 0, 0, -2^54, 0)
  # Directly enumerate surviving equal-distance group members. The table and signal products
  # are dyadic; the expected values do not rely on the block-mean identity.
  groups <- c(1L, 2L, 1L, 3L, 1L, 4L, 1L, 5L)
  signal <- c(1, 0, 1, 0, 2, 0, 1, 0)
  expected <- vapply(1:8, function(i) {
    kept <- setdiff(1:8, i)
    ranks <- kept - (kept > i)
    weights <- table[ranks]
    for (g in unique(groups[kept])) {
      ids <- which(groups[kept] == g)
      weights[ids] <- sum(weights[ids]) / length(ids)
    }
    sum(weights * signal[kept])
  }, numeric(1L))
  original <- options(RdistNN.use_compiled = TRUE)
  on.exit(options(original), add = TRUE)
  for (native in c(FALSE, TRUE)) {
    options(RdistNN.use_compiled = native)
    jk <- RdistNN:::.dnn_loo_centered(signal, matrix(1:8), table, groups = groups)
    actual <- RdistNN:::.dnn_restore_loo(jk$loo, jk$offset)
    expect_equal(as.numeric(actual), expected, tolerance = 1e-15)
  }
})

test_that("interleaved block deletions retain the public observation-order mass", {
  ranks <- matrix(c(1L, 9L, 3L, 7L, 5L, 2L, 4L, 6L, 8L, 10L))
  groups <- c(rep(1L, 5L), 2:6)
  table <- c(2^112, 0, 2^56, 0, 1, 0, -2^56, 0, -2^112)
  signal <- c(0, 1, 1, 1, 1, 0, 0, 0, 0, 0)
  # Original observations give the block entries A,-A,B,-B,1, whose exact
  # mass is one. Reducing in rank order A,B,1,-B,-A loses that residual.
  X <- matrix(c(rep(0, 5L), 1:5))
  rw <- as_rank_weight(function(r) table[r])
  W <- suppressWarnings(weights_matrix(X[-10L, , drop = FALSE], x0 = 0,
    ranking = list(ranks = ranks[-10L, , drop = FALSE],
      dx = matrix(abs(X[-10L, 1L]), ncol = 1L)), weight = rw))$W
  expect_equal(as.numeric(W[1:5, ]), rep(1 / 5, 5L), tolerance = 0)
  expect_equal(sum(W * signal[-10L]), 4 / 5, tolerance = 0)
  original <- options(RdistNN.use_compiled = TRUE)
  on.exit(options(original), add = TRUE)
  for (native in c(FALSE, TRUE)) {
    options(RdistNN.use_compiled = native)
    expect_identical(is.null(RdistNN:::.rdistnn_compiled_fun("rustat_block_rank_loo")), !native)
    for (direct in c(FALSE, TRUE)) {
      actual <- RdistNN:::.dnn_block_loo(matrix(signal), ranks[, 1L],
        table, groups, direct = direct)
      expect_equal(actual[10L, 1L], 4 / 5, tolerance = 0)
    }
  }
})

test_that("contiguous supplied blocks retain original member-order masses", {
  ranks <- c(1L, 5L, 2L, 4L, 3L, 6L)
  groups <- c(rep(1L, 5L), 2L)
  table <- c(2^112, 2^56, 1, -2^56, -2^112)
  signal <- matrix(c(0, 1, 1, 1, 1, 0))
  original <- options(RdistNN.use_compiled = TRUE)
  on.exit(options(original), add = TRUE)
  for (native in c(FALSE, TRUE)) {
    options(RdistNN.use_compiled = native)
    expect_identical(is.null(RdistNN:::.rdistnn_compiled_fun("rustat_block_rank_loo")), !native)
    for (direct in c(FALSE, TRUE)) {
      actual <- RdistNN:::.dnn_block_loo(signal, ranks, table, groups, direct = direct)
      expect_equal(actual[6L, 1L], 4 / 5, tolerance = 0)
    }
  }
})
