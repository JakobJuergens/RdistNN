test_that("CATE concentration recovers finite weighted residuals after subtraction overflow", {
  old <- options(RdistNN.use_compiled = TRUE)
  on.exit(options(old), add = TRUE)
  for (native in c(FALSE, TRUE)) {
    options(RdistNN.use_compiled = native)
    for (amplitude in c(6e307, 1.6e308)) {
      signal <- rep(-amplitude, 12L)
      signal[11L] <- if (amplitude == 6e307) 8e307 else amplitude
      fit <- cate_dnn(X = matrix(seq_len(12L)), x0 = 0, pseudo_outcome = signal,
        s = 2L, se = FALSE, return = "W")
      # w_i=(12-i)/66. Concentration coefficients are 65 at row11 and
      # 11,10,...,2 at the other active rows: total130, largest five103.
      for (object in list(fit, unserialize(serialize(fit, NULL)))) {
        d <- diagnostics(object)
        expect_true(d$score_contribution_available)
        expect_equal(d$top5_score_share, 103 / 130, tolerance = 1e-15)
      }
    }
  }
})

test_that("ordinary score concentration preserves its existing arithmetic", {
  w <- c(.2, .1, .3, .15, .05, .2)
  left <- c(2, -1, 7, 0, 4, -3)
  right <- rep(.75, 6L)
  contribution <- abs(w * (left - right))
  expected <- sum(sort(contribution, decreasing = TRUE)[1:5]) / sum(contribution)
  share <- getFromNamespace(".dnn_score_contributions", "RdistNN")(w, left, right)
  expect_true(share$available)
  expect_identical(share$share, expected)
  expect_identical(share$contribution, w * (left - right))
})

test_that("score concentration normalizes overflowing totals and tiny weighted residuals", {
  share <- getFromNamespace(".dnn_score_contributions", "RdistNN")
  total <- share(rep(1, 6L), rep(.Machine$double.xmax / 2, 6L), rep(0, 6L))
  expect_true(total$available)
  expect_equal(total$share, 5 / 6, tolerance = 1e-15)
  tiny <- share(c(2^-1074, rep(1, 6L)), c(-1.6e308, rep(2^-50, 6L)),
    c(8e307, rep(0, 6L)))
  expect_true(tiny$available)
  # Fraction.from_float products/differences independently give this ratio;
  # normalizing the tiny weight before the difference would erase its mass.
  expect_equal(tiny$share, 0.7273363439834166, tolerance = 1e-15)
  expect_equal(tiny$contribution[1L], -1.1857575500189916e-15, tolerance = 1e-15)
  unrepresentable <- share(rep(1, 6L), rep(-1.6e308, 6L), rep(1.6e308, 6L))
  expect_true(unrepresentable$available)
  expect_equal(unrepresentable$share, 5 / 6, tolerance = 1e-15)
  expect_true(all(is.na(unrepresentable$contribution)))
})

test_that("score concentration excludes zero support before checking residuals", {
  share <- getFromNamespace(".dnn_score_contributions", "RdistNN")
  inactive <- share(c(1, 0), c(1, Inf), c(1, NA_real_))
  expect_true(inactive$available)
  expect_true(is.na(inactive$share))
  expect_identical(inactive$contribution, c(0, 0))
  absent <- share(rep(0, 2L), c(Inf, NA_real_), c(NA_real_, Inf))
  expect_true(absent$available)
  expect_true(is.na(absent$share))
  invalid <- share(c(1, 1), c(1, 1), c(0, Inf))
  expect_false(invalid$available)
  expect_true(is.na(invalid$share))
  expect_identical(invalid$contribution, c(1, NA_real_))
})

test_that("influence plot data retain finite signed CATE and CLATE contributions", {
  skip_if_not_installed("ggplot2")
  old <- options(RdistNN.use_compiled = TRUE)
  on.exit(options(old), add = TRUE)
  signal <- rep(-1.6e308, 12L)
  signal[11L] <- 1.6e308
  for (native in c(FALSE, TRUE)) {
    options(RdistNN.use_compiled = native)
    cate <- cate_dnn(X = matrix(seq_len(12L)), x0 = 0, pseudo_outcome = signal,
      s = 2L, se = FALSE, return = "W")
    A <- rep(0:1, 6L)
    d1 <- rep(1, 12L)
    d1[12L] <- 0
    clate <- clate_dnn(X = matrix(seq_len(12L)), x0 = 0,
      Y = A * signal, D = A, A = A, e = rep(.5, 12L),
      m_y0 = rep(0, 12L), m_y1 = signal, m_d0 = rep(0, 12L), m_d1 = d1,
      s = 2L, se = FALSE, return = "W")
    for (fit in list(cate, clate)) {
      center <- if (inherits(fit, "rdistnn_cate_result")) fit$tau else fit$theta
      # Binary-normalized subtraction is an independent dimensional reference.
      expected <- as.numeric(fit$W) * (signal / 2^1022 - center / 2^1022) * 2^1022
      expected[12L] <- 0
      for (object in list(fit, unserialize(serialize(fit, NULL)))) {
        p <- plot_influence(object)
        expect_true(all(is.finite(p$data$contribution)))
        expect_equal(p$data$contribution, expected[p$data$obs_id], tolerance = 1e-14)
        expect_identical(p$data$contribution[p$data$obs_id == 12L], 0)
      }
    }
  }
})
