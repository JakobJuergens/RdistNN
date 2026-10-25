# These regressions catch cancellation in Beta CDF increments and in the
# two-scale coefficient denominator. Decimal literals below were calculated
# independently with Python Decimal at 100-digit precision.

test_that("selecting every subset neighbor is exactly uniform under either law", {
  # Summing all s order-statistic densities gives s on [0,1]. Thus every
  # rank receives 1/n, independently of s, without subtracting rounded CDFs.
  for (n in c(1L, 7L, 31L, 257L)) {
    for (s in unique(c(1L, max(1L, n %/% 2L), n))) {
      ranks <- rev(seq_len(n))
      uniform <- rep(1 / n, n)
      for (normalize in c(FALSE, TRUE)) {
        direct <- dnn_weights_asymp_full_from_ranks(ranks, s, s, normalize)
        expect_identical(direct$w_rank, uniform)
        expect_identical(direct$w, uniform)
        rw <- rank_weight_dnn_asymp(s, s, normalize)
        expected <- if (normalize) uniform / .compensated_sum(uniform) else uniform
        expect_identical(rank_weight_eval(rw, ranks), expected)
        expect_identical(rank_weight_eval_batch(rw, cbind(ranks, rev(ranks))),
          matrix(expected, nrow = n, ncol = 2L))
      }
    }
    for (normalize in c(FALSE, TRUE)) {
      expect_identical(.dnn_rank_table_asymp_k1(n, 1L, normalize), rep(1 / n, n))
    }
    if (n >= 2L) {
      # dim=2 and scales (1,2) have coefficients (-1,2). The first component
      # must retain the same exact uniform identity before signed blending.
      j <- seq_len(n)
      second <- .beta_interval_mass((j - 1) / n, j / n, 1L, 2L)
      expect_identical(tdnn_rank_table_asymp(n, 1L, 2L, 2L, FALSE),
        -rep(1 / n, n) + 2 * second)
    }
  }
})

test_that("asymptotic DNN and TDNN retain represented right-tail mass", {
  want <- 1.6542999783940550651810050759970065291673e-68
  for (normalize in c(FALSE, TRUE)) {
    out <- dnn_weights_asymp_full_from_ranks(100:1, 1L, 100L, normalize)
    expect_gt(out$w_rank[80L], 0)
    expect_equal(out$w_rank[80L] / want, 1, tolerance = 5e-14)
    expect_identical(out$w, out$w_rank[100:1])
    single <- RdistNN:::.dnn_rank_table_asymp_k1(100L, 100L, normalize)
    expect_equal(single[80L] / want, 1, tolerance = 5e-14)
    expect_true(all(single > 0))
    expect_equal(sum(single), 1, tolerance = 2e-15)
  }
})

test_that("adjacent TDNN orders use accurate finite bias coefficients", {
  fixtures <- list(
    c(1000000000, 1000000001, 2, -1000000000),
    c(1000000000, 1000000001, 3, -1500000000.2499999999305555555903),
    c(1000, 1001, 3, -1500.2499305902563162827243569130),
    c(2147483646, 2147483647, 2147483647, -2305843006529339392.2083333333)
  )
  for (f in fixtures) {
    sc <- RdistNN:::.tdnn_scalars(f[1], f[2], f[3])
    expect_true(all(is.finite(sc)))
    expect_equal(unname(sc[1]) / f[4], 1, tolerance = 5e-15)
    expect_identical(unname(sc[2]), 1 - unname(sc[1]))
    # A unit difference is representable only below 2^53.
    if (abs(sc[1]) < 2^53) expect_equal(sum(sc), 1, tolerance = 0)
  }
})

test_that("Beta interval masses retain both tails, endpoints and tiny probabilities", {
  mass <- RdistNN:::.beta_interval_mass
  want <- 1.6542999783940550651810050759970065291673e-68
  got <- mass(c(.79, .2), c(.8, .21), c(1, 100), c(100, 1))
  expect_true(all(got > 0))
  expect_equal(got / want, c(1, 1), tolerance = 5e-14)
  expect_equal(mass(c(0, 0, 1), c(0, 1, 1), 7, 9), c(0, 1, 0), tolerance = 0)
  expect_equal(mass(0, 1e-200, 1, 1) / 1e-200, 1, tolerance = 5e-14)
  # Beta(1,1) is uniform, including intervals crossing the crossover.
  expect_equal(mass(c(0, .49, .5, .75), c(.25, .51, .75, 1), 1, 1),
               c(.25, .02, .25, .25), tolerance = 1e-15)
  # Every order-statistic shape partitions the whole probability space.
  for (a in 1:25) {
    p <- mass((0:99) / 100, (1:100) / 100, a, 26 - a)
    expect_true(all(p >= 0))
    expect_equal(sum(p), 1, tolerance = 2e-15)
  }
})

test_that("all k and normalization choices preserve the Beta mixture", {
  # For integer shapes F_Beta(a,s-a+1)(x) = P(Binomial(s,x) >= a).
  # Integrate its polynomial by evaluating its independent binomial sum.
  n <- 9L
  s <- 7L
  cdf <- function(x, a) {
    sum(vapply(a:s, function(i) choose(s, i) * x^i * (1 - x)^(s - i), 0.0))
  }
  for (k in seq_len(s)) {
    want <- vapply(seq_len(n), function(j) {
      mean(vapply(seq_len(k), function(a) cdf(j / n, a) - cdf((j - 1) / n, a), 0.0))
    }, 0.0)
    for (normalize in c(FALSE, TRUE)) {
      got <- dnn_weights_asymp_full_from_ranks(seq_len(n), k, s, normalize)$w_rank
      expect_true(all(got > 0))
      expect_equal(got, want, tolerance = 2e-15)
      expect_equal(sum(got), 1, tolerance = 2e-15)
    }
  }
  # Decimal binomial-sum fixture for k > 1 retains a far-right mixture tail.
  want <- 4.9760024804384519642649082785548295536855e-52
  for (normalize in c(FALSE, TRUE)) {
    rw <- rank_weight_dnn_asymp(10L, 100L, normalize)
    R <- cbind(1:100, 100:1, 1:100)
    W <- rank_weight_eval_batch(rw, R)
    expect_equal(W[80, 1] / want, 1, tolerance = 1e-13)
    expect_equal(W[, 2], rev(W[, 1]), tolerance = 0)
    expect_identical(W[, 1], W[, 3])
    expect_equal(rank_weight_eval(rw, 1:100), W[, 1], tolerance = 0)
    expect_equal(colSums(W), rep(1, 3), tolerance = 2e-15)
    uniform <- dnn_weights_asymp_full_from_ranks(1:100, 100L, 100L, normalize)$w_rank
    expect_equal(uniform / .01, rep(1, 100), tolerance = 2e-14)
  }
})

test_that("TDNN exact and asymptotic tables use stable coefficients and signed tails", {
  sc <- -1500.2499305902563162827243569130
  # For n=s2=s1+1 only the first two exact ranks have mass.
  want <- c(1 - sc / 1001, sc / 1001, rep(0, 999))
  exact <- rank_weight_eval(rank_weight_tdnn_exact(1000L, 1001L, 3L), 1:1001)
  expect_equal(exact, want, tolerance = 2e-9)
  expect_equal(unname(RdistNN:::.tdnn_scalars(1000L, 1001L, 3L)[1]) / sc, 1, tolerance = 5e-15)
  # dim=2, s2=2*s1 gives coefficients (-1,2). The tail is negative:
  # -(.21^50-.2^50) + 2*(.21^100-.2^100).
  tail <- -1.1785244443664566263864564756466283356941e-34
  for (normalize in c(FALSE, TRUE)) {
    out <- rank_weight_eval(rank_weight_tdnn_asymp(50L, 100L, 2L, normalize = normalize), 100:1)
    table <- tdnn_rank_table_asymp(100L, 50L, 100L, 2L, normalize)
    expect_lt(table[80], 0)
    expect_equal(table[80] / tail, 1, tolerance = 5e-14)
    expect_identical(out, table[100:1])
    rw <- rank_weight_tdnn_asymp(50L, 100L, 2L, normalize)
    W <- rank_weight_eval_batch(rw, cbind(1:100, 100:1))
    expect_equal(W[, 1], table, tolerance = 0)
    expect_equal(W[, 2], out, tolerance = 0)
    expect_equal(colSums(W), c(1, 1), tolerance = 2e-15)
  }
  # An ordinary small table has the independent polynomial reference.
  j <- 1:6
  p2 <- ((6 - j + 1) / 6)^2 - ((6 - j) / 6)^2
  p4 <- ((6 - j + 1) / 6)^4 - ((6 - j) / 6)^4
  expect_equal(tdnn_rank_table_asymp(6L, 2L, 4L, 2L), -p2 + 2 * p4,
               tolerance = 2e-15)
})
