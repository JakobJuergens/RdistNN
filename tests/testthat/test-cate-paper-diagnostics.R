test_that("CATE paper diagnostics use physical Euclidean radii and retained scores", {
  .with_distance_tie_warning({
  X <- rbind(c(-1, 0), c(1, 0), c(0, -1), c(0, 1))
  fit <- suppressWarnings(
    cate_dnn(
      X = X, x0 = c(0, 0), pseudo_outcome = 1:4,
      k = 1L, s = 1L, se = FALSE, return = "W"
    )
  )
  diag <- diagnostics(fit)

  expect_equal(diag$physical_radius_90, 1, tolerance = 1e-12)
  expect_false(diag$arm_ess_available)
  expect_true(is.na(diag$treated_ess))
  expect_true(diag$score_contribution_available)
  expect_true(is.finite(diag$top5_score_share))
  })
})

test_that("CATE arm ESS distinguishes tiny positive and zero local mass", {
  dat <- data.frame(
    x = 1:10,
    y = (1:10) / 10,
    D = as.integer(1:10 == 6L),
    e = rep(0.5, 10),
    m0 = 0,
    m1 = 1
  )
  fit <- cate_dnn(
    data = dat, x_cols = "x", y_col = "y", t_col = "D",
    e_col = "e", m0_col = "m0", m1_col = "m1", x0 = 0,
    k = 1L, s = 5L, se = FALSE, return = "W"
  )
  diag <- diagnostics(fit)
  expect_true(diag$arm_ess_available)
  expect_gt(diag$treated_mass, 0)
  expect_lt(diag$treated_mass, 0.01)
  expect_equal(diag$treated_ess, 1, tolerance = 1e-12)

  dat$D <- as.integer(1:10 == 10L)
  fit_zero <- cate_dnn(
    data = dat, x_cols = "x", y_col = "y", t_col = "D",
    e_col = "e", m0_col = "m0", m1_col = "m1", x0 = 0,
    k = 1L, s = 5L, se = FALSE, return = "W"
  )
  diag_zero <- diagnostics(fit_zero)
  expect_equal(diag_zero$treated_mass, 0, tolerance = 1e-12)
  expect_true(is.na(diag_zero$treated_ess))
})

test_that("CATE top-five score share detects concentrated contributions", {
  X <- matrix(seq_len(20), ncol = 1L)
  signal <- c(1e6, rep(0, 19))
  fit <- cate_dnn(
    X = X, x0 = 0, pseudo_outcome = signal,
    k = 1L, s = 5L, se = FALSE, return = "W"
  )
  diag <- diagnostics(fit)

  expect_gt(diag$top5_score_share, 0.85)
  expect_gt(abs(diag$top5_score_share - diag$top5_weight_share), 0.01)
})

test_that("CATE within-arm ESS recovers zero and subnormal squared mass", {
  old <- options(RdistNN.use_compiled = TRUE)
  on.exit(options(old), add = TRUE)
  for (compiled in c(TRUE, FALSE)) {
    options(RdistNN.use_compiled = compiled)
    for (power in c(537, 600)) for (flip in c(FALSE, TRUE)) {
      # Independent binary fixture: the arm ratio is 2:1, hence ESS is 9/5.
      # At exponent 537 the first square is positive but the second is lost.
      w <- c(1, 2^-power, 2^(-power - 1))
      arm <- if (flip) c(1, 0, 0) else c(0, 1, 1)
      d <- RdistNN:::.cate_paper_diagnostics(matrix(w), matrix(1:3),
        matrix(0), rep(0, 3), 0, arm)
      value <- if (flip) d$control_ess else d$treated_ess
      expect_equal(value, 9/5, tolerance = 1e-14)
      normal <- RdistNN:::.cate_paper_diagnostics(matrix(c(1, .5, .25)),
        matrix(1:3), matrix(0), rep(0, 3), 0, arm)
      expect_equal(value, if (flip) normal$control_ess else normal$treated_ess,
        tolerance = 1e-14)
      singleton <- RdistNN:::.cate_paper_diagnostics(matrix(c(1, w[2L], 0)),
        matrix(1:3), matrix(0), rep(0, 3), 0, arm)
      expect_identical(if (flip) singleton$control_ess else singleton$treated_ess, 1)
    }
    d <- RdistNN:::.cate_paper_diagnostics(matrix(c(1, 0, 0)),
      matrix(1:3), matrix(0), rep(0, 3), 0, c(0, 1, 1))
    expect_true(is.na(d$treated_ess))
    expect_identical(d$control_ess, 1)
    w <- c(.5, .25, .125)
    # The ordinary branch retains the pre-existing expression, bit for bit.
    d <- RdistNN:::.cate_paper_diagnostics(matrix(w), matrix(1:3),
      matrix(0), rep(0, 3), 0, c(0, 1, 1))
    normalized <- w / sum(w)
    expect_identical(d$treated_ess,
      sum(normalized[2:3])^2 / sum(normalized[2:3]^2))
  }
})

test_that("CATE singleton tail-arm ESS survives public fit and current reuse", {
  old <- options(RdistNN.use_compiled = TRUE)
  on.exit(options(old), add = TRUE)
  tuning <- new_tuning("arm_ess_boundary",
    eval = function(...) list(s = 300L),
    eval_grid = function(...) data.frame(s = c(299L, 300L)))
  for (compiled in c(TRUE, FALSE)) for (flip in c(FALSE, TRUE)) {
    options(RdistNN.use_compiled = compiled)
    arm <- as.integer(seq_len(600L) == 301L)
    if (flip) arm <- 1L - arm
    fit <- cate_dnn(data = data.frame(x = seq_len(600L), y = 0,
        arm = arm, e = .5, m0 = 0, m1 = 0), x_cols = "x", x0 = 0,
      y_col = "y", t_col = "arm", e_col = "e", m0_col = "m0", m1_col = "m1",
      k = 1L, tuning = tuning, tuning_mode = "grid", se = FALSE, return = "W")
    restored <- unserialize(serialize(fit, NULL))
    for (view in list(fit, restored, candidate_fit(restored, 2L))) {
      d <- diagnostics(view)
      value <- if (flip) d$control_ess else d$treated_ess
      mass <- if (flip) d$control_mass else d$treated_mass
      expect_gt(mass, 0)
      expect_equal(value, 1, tolerance = 0)
      expect_equal(predict(view, x0 = 0), view$tau, tolerance = 0)
      expect_equal(d, diagnostics(fit), tolerance = 0)
    }
  }
})
