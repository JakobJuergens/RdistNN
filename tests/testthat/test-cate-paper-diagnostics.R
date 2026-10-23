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
