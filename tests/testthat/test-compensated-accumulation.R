test_that("ordered matrix products retain cancellation and row identity", {
  x <- rbind(cancel = c(2^54, 1, -2^54), ordinary = c(1, 2, 3))
  B <- cbind(sum = c(1, 1, 1), first = c(1, 0, 0))
  expected <- matrix(c(1, 6, 2^54, 1), 2, 2,
    dimnames = list(rownames(x), colnames(B)))
  old <- options(RdistNN.use_compiled = TRUE)
  on.exit(options(old), add = TRUE)
  for (native in c(TRUE, FALSE)) {
    options(RdistNN.use_compiled = native)
    expect_identical(is.null(.rdistnn_compiled_fun("rustat_compensated_crossprod")), !native)
    value <- .compensated_matprod(x, B)
    expect_identical(value, expected)
    expect_identical(.compensated_matprod(x[c(2, 1, 2), ], B), value[c(2, 1, 2), ])
    for (i in 1:2) expect_identical(.compensated_matprod(x[i, , drop = FALSE], B),
      value[i, , drop = FALSE])
    expect_identical(dim(.compensated_matprod(x[FALSE, , drop = FALSE], B)), c(0L, 2L))
  }
})

test_that("weighted sums retain small signed contributions in both engines", {
  old <- options(RdistNN.use_compiled = TRUE)
  on.exit(options(old), add = TRUE)
  A <- 2^54
  for (native in c(TRUE, FALSE)) {
    options(RdistNN.use_compiled = native)
    if (native) {
      expect_false(is.null(.rdistnn_compiled_fun("rustat_compensated_sum")))
      expect_false(is.null(.rdistnn_compiled_fun("rustat_compensated_crossprod")))
    }
    expect_identical(.compensated_sum(c(A, 1, -A)), 1)
    expect_identical(.compensated_sum(c(1, 2^200, 1, -2^200)), 2)
    expect_identical(.compensated_sum(numeric()), 0)
    x <- cbind(signal = c(A, 1, -A), constant = rep(1, 3))
    w <- cbind(first = rep(1 / 3, 3), second = c(0, 1, 0))
    expected <- matrix(c(1 / 3, 1, 1, 1), 2, 2,
      dimnames = list(colnames(x), colnames(w)))
    expect_identical(.compensated_crossprod(x, w), expected)
    expect_identical(.dnn_estimate_from_signal(x, w, 3), expected)
    # Compensation addresses accumulation; it does not change product rounding.
    expect_identical(unname(.compensated_crossprod(
      matrix(c(1 + 2^-27, 1)), matrix(c(1 - 2^-27, -1)))), matrix(0))
  }
})

test_that("deletion sums preserve corrections until prefix and suffix merge", {
  old <- options(RdistNN.use_compiled = TRUE)
  on.exit(options(old), add = TRUE)
  A <- 2^54
  ranks <- cbind(1:4, 4:1, 1:4)
  signal <- cbind(c(A, 1, 0, -A), rep(7, 4))
  expected <- cbind(c(-A, 0, 1, A), rep(21, 4))
  results <- list()
  for (native in c(TRUE, FALSE)) {
    options(RdistNN.use_compiled = native)
    if (native) expect_false(is.null(.rdistnn_compiled_fun("rustat_jackknife_rank_full")))
    out <- jackknife_rank_full(ranks, signal, rep(1, 4), rep(1, 3))
    expect_identical(unname(out$theta0), rbind(rep(1, 3), rep(28, 3)))
    expect_identical(out$theta_loo[, , 1], expected)
    expect_identical(out$theta_loo[, , 2], expected)
    expect_identical(out$theta_loo[, , 3], expected)
    results[[length(results) + 1L]] <- out
    expect_identical(.dnn_loo_prefix(signal, ranks, rep(1, 3))[, , 1], expected)
  }
  expect_identical(results[[1]], results[[2]])
})

test_that("compensated deletion edges preserve zero support and small samples", {
  old <- options(RdistNN.use_compiled = TRUE)
  on.exit(options(old), add = TRUE)
  for (native in c(TRUE, FALSE)) {
    options(RdistNN.use_compiled = native)
    out <- jackknife_rank_full(matrix(1L), matrix(9), 1, numeric())
    expect_identical(as.numeric(out$theta_loo), 0)
    out <- jackknife_rank_full(matrix(1:2), matrix(c(2, 3)), c(.5, .5), 1)
    expect_identical(as.numeric(out$theta_loo), c(3, 2))
    out <- jackknife_rank_full(matrix(1:4), matrix(c(1, 2, 3, 2^800)),
      c(1, 0, 0, 0), c(1, 0, 0))
    expect_identical(as.numeric(out$theta_loo), c(2, 1, 1, 1))
  }
})

test_that("finite-support jackknife engines share stable weighted arithmetic", {
  old <- options(RdistNN.use_compiled = TRUE)
  on.exit(options(old), add = TRUE)
  A <- 2^54
  neighbors <- matrix(1:4, nrow = 1L)
  signal <- matrix(c(A, 1, -A, 0))
  for (native in c(TRUE, FALSE)) {
    options(RdistNN.use_compiled = native)
    for (engine in list(jackknife_fast_rank, jackknife_slow_rank)) {
      out <- engine(neighbors, rep(1, 3), signal, normalize = FALSE)
      expect_identical(out$theta0, 1)
      expect_identical(as.numeric(out$theta_loo), c(-A, 0, A, 1))
    }
  }
})

test_that("retained deletion states survive reverse-buffer boundaries", {
  old <- options(RdistNN.use_compiled = TRUE)
  on.exit(options(old), add = TRUE)
  n <- 4101L
  z <- numeric(n)
  z[c(1L, 2048L, 2050L, n)] <- c(2^54, 1, -2^54, 1)
  ranks <- cbind(seq_len(n), rev(seq_len(n)))
  options(RdistNN.use_compiled = TRUE)
  native <- jackknife_rank_full(ranks, cbind(z, 7), rep(1, n), rep(1, n - 1L))
  options(RdistNN.use_compiled = FALSE)
  interpreted <- jackknife_rank_full(ranks, cbind(z, 7), rep(1, n), rep(1, n - 1L))
  expect_identical(native, interpreted)
  expect_identical(native$theta_loo[2049L, 1L, ], c(2, 2))
  expect_identical(native$theta_loo[2049L, 2L, ], rep(7 * (n - 1L), 2))
})

test_that("public means and restored predictions share compensated weighting", {
  .with_distance_tie_warning({
  old <- options(RdistNN.use_compiled = TRUE)
  on.exit(options(old), add = TRUE)
  x <- matrix(1:4, ncol = 1L)
  targets <- matrix(c(0, 2, 0), ncol = 1L)
  signal <- c(2^54, 1, -2^54, 1)
  for (native in c(TRUE, FALSE)) {
    options(RdistNN.use_compiled = native)
    np <- suppressWarnings(np_dnn(X = x, X0 = targets, Y = signal, k = 4, s = 4, se = FALSE))
    cate <- suppressWarnings(cate_dnn(X = x, X0 = targets, pseudo_outcome = signal,
      k = 4, s = 4, se = FALSE))
    ratio <- suppressWarnings(affine_ratio_dnn(X = x, X0 = targets,
      numerator = signal, denominator = rep(2, 4), k = 4, s = 4, se = FALSE))
    expect_identical(as.numeric(np$mu), rep(.5, 3))
    expect_identical(as.numeric(cate$tau), rep(.5, 3))
    expect_identical(as.numeric(ratio$theta), rep(.25, 3))
    expect_identical(as.numeric(predict(unserialize(serialize(np, NULL)), X0 = targets)), rep(.5, 3))
    expect_identical(as.numeric(predict(unserialize(serialize(cate, NULL)), X0 = targets)), rep(.5, 3))
  }
  })
})


test_that("weighted diagnostics and GMM moments share compensated means", {
  old <- options(RdistNN.use_compiled = TRUE)
  on.exit(options(old), add = TRUE)
  x <- matrix(1:4, ncol = 1L)
  signal <- c(2^54, 1, -2^54, 1)
  model <- list(moment_obs = function(sample, theta, ctx) cbind(signal),
    jacobian_obs = function(sample, theta, ctx) cbind(signal))
  for (native in c(TRUE, FALSE)) {
    options(RdistNN.use_compiled = native)
    np <- suppressWarnings(np_dnn(X = x, X0 = matrix(0), Y = signal, k = 4, s = 4, se = FALSE))
    expect_identical(as.numeric(diagnostics(np)$weighted_mean_y), .5)
    expect_identical(.gmm_gbar(model, NULL, 0, NULL, rep(.25, 4), 4), .5)
    expect_identical(.gmm_eval_jacobian(model, NULL, 0, NULL, rep(.25, 4), 4,
      1, 1, "analytic", 1e-5), matrix(.5))
  }
})

test_that("finite-support references preserve repeated-index fast semantics", {
  old <- options(RdistNN.use_compiled = TRUE)
  on.exit(options(old), add = TRUE)
  ranks <- matrix(c(1L, 1L, 2L), nrow = 1L)
  for (native in c(TRUE, FALSE)) {
    options(RdistNN.use_compiled = native)
    for (engine in list(jackknife_fast_rank, jackknife_slow_rank)) {
      out <- engine(ranks, c(1, 1), matrix(c(2, 6)), normalize = FALSE)
      expect_identical(out$theta0, 4)
      # Each occurrence adds its single-position replacement change: 4+4+4.
      expect_identical(as.numeric(out$theta_loo), c(12, 4))
    }
  }
})
