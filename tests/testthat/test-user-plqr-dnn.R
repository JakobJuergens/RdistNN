test_that("plqr_dnn matrix API matches direct g-weighted moment identity", {
  d <- .make_plqr_data()
  out <- plqr_dnn(
    X = d$X, D = d$D, g = d$g, mD = d$mD, Y = d$Y,
    theta_eval = d$theta_eval, tau = d$tau,
    X0 = d$X0, distance = "l2",
    method = "exact", k = 2L, s = 3L, se = FALSE,
    return = "W"
  )

  score <- (d$tau - as.numeric(d$Y <= as.vector(d$D %*% d$theta_eval) + d$g)) * (d$D - d$mD)
  manual <- crossprod(score, out$W)

  expect_s3_class(out, "rdistnn_plqr_result")
  expect_equal(dim(out$moment), c(ncol(d$D), nrow(d$X0)))
  expect_equal(unname(out$moment), unname(manual), tolerance = 1e-12)
})

test_that("plqr_dnn returns finite jackknife SEs with expected shape", {
  d <- .make_plqr_data(n = 50L, p = 2L, p_d = 2L, m = 2L, tau = 0.5)
  out <- plqr_dnn(
    X = d$X, D = d$D, g = d$g, mD = d$mD, Y = d$Y,
    theta_eval = d$theta_eval, tau = d$tau,
    X0 = d$X0,
    k = 2L, s = 3L, method = "exact", se = TRUE, jackknife = "loo"
  )

  expect_true(is.matrix(out$moment_se))
  expect_equal(dim(out$moment_se), dim(out$moment))
  expect_true(all(is.finite(out$moment_se)))
  expect_true(all(out$moment_se >= 0))
})

test_that("default jackknife mode uses fast path and matches loo reference", {
  d <- .make_plqr_data(n = 45L, p = 2L, p_d = 2L, m = 2L, tau = 0.5)

  out_default <- plqr_dnn(
    X = d$X, D = d$D, g = d$g, mD = d$mD, Y = d$Y,
    theta_eval = d$theta_eval, tau = d$tau, X0 = d$X0,
    k = 2L, s = 3L, method = "exact", se = TRUE
  )
  out_fast <- plqr_dnn(
    X = d$X, D = d$D, g = d$g, mD = d$mD, Y = d$Y,
    theta_eval = d$theta_eval, tau = d$tau, X0 = d$X0,
    k = 2L, s = 3L, method = "exact", se = TRUE, jackknife = "fast"
  )
  out_loo <- plqr_dnn(
    X = d$X, D = d$D, g = d$g, mD = d$mD, Y = d$Y,
    theta_eval = d$theta_eval, tau = d$tau, X0 = d$X0,
    k = 2L, s = 3L, method = "exact", se = TRUE, jackknife = "loo"
  )

  expect_equal(out_default$moment, out_fast$moment, tolerance = 1e-12)
  expect_equal(out_default$moment_se, out_fast$moment_se, tolerance = 1e-12)
  expect_equal(out_fast$moment_se, out_loo$moment_se, tolerance = 1e-12)
})

test_that("asymp method fast and loo jackknife paths agree", {
  d <- .make_plqr_data(n = 48L, p = 2L, p_d = 2L, m = 2L, tau = 0.5)

  out_fast <- plqr_dnn(
    X = d$X, D = d$D, g = d$g, mD = d$mD, Y = d$Y,
    theta_eval = d$theta_eval, tau = d$tau, X0 = d$X0,
    method = "asymp", k = 2L, s = 3L, se = TRUE, jackknife = "fast"
  )
  out_loo <- plqr_dnn(
    X = d$X, D = d$D, g = d$g, mD = d$mD, Y = d$Y,
    theta_eval = d$theta_eval, tau = d$tau, X0 = d$X0,
    method = "asymp", k = 2L, s = 3L, se = TRUE, jackknife = "loo"
  )

  expect_equal(out_fast$moment, out_loo$moment, tolerance = 1e-12)
  expect_equal(out_fast$moment_se, out_loo$moment_se, tolerance = 1e-12)
})

test_that("return extras are populated consistently for plqr_dnn", {
  d <- .make_plqr_data(n = 55L, p = 2L, p_d = 2L, m = 3L, tau = 0.5)
  out <- plqr_dnn(
    X = d$X, D = d$D, g = d$g, mD = d$mD, Y = d$Y,
    theta_eval = d$theta_eval, tau = d$tau, X0 = d$X0,
    method = "exact", k = 3L, s = 4L, se = FALSE,
    return = c("W", "ranks", "ord", "dx", "score", "ess", "neighbors")
  )

  expect_equal(dim(out$W), c(nrow(d$X), nrow(d$X0)))
  if (!is.null(out$ranks)) expect_equal(dim(out$ranks), c(nrow(d$X), nrow(d$X0)))
  if (!is.null(out$ord)) expect_equal(dim(out$ord), c(nrow(d$X), nrow(d$X0)))
  if (!is.null(out$dx)) expect_equal(dim(out$dx), c(nrow(d$X), nrow(d$X0)))
  expect_equal(dim(out$score), dim(d$D))
  expect_equal(unname(out$moment), unname(crossprod(out$score, out$W)), tolerance = 1e-12)
  expect_equal(out$ess, 1 / colSums(out$W^2), tolerance = 1e-12)
  expect_equal(dim(out$neighbors), c(3L, nrow(d$X0)))
  expect_true(all(out$neighbors >= 1 & out$neighbors <= nrow(d$X)))
  support_edge <- nrow(d$X) - out$meta$s + out$meta$k
  for (j in seq_len(nrow(d$X0))) {
    expect_equal(out$neighbors[, j], order(out$ranks[, j])[seq_len(out$meta$k)])
    expect_equal(sum(out$W[, j] > 0), support_edge)
    expect_gt(sum(out$W[, j] > 0), nrow(out$neighbors))
  }
})

test_that("x0 path matches corresponding single target from X0", {
  d <- .make_plqr_data(n = 44L, p = 2L, p_d = 2L, m = 4L, tau = 0.5)
  j <- 2L

  out_X0 <- plqr_dnn(
    X = d$X, D = d$D, g = d$g, mD = d$mD, Y = d$Y,
    theta_eval = d$theta_eval, tau = d$tau, X0 = d$X0,
    k = 2L, s = 3L, method = "exact", se = FALSE
  )
  out_x0 <- plqr_dnn(
    X = d$X, D = d$D, g = d$g, mD = d$mD, Y = d$Y,
    theta_eval = d$theta_eval, tau = d$tau, x0 = d$X0[j, ],
    k = 2L, s = 3L, method = "exact", se = FALSE
  )

  expect_equal(
    unname(out_x0$moment),
    unname(out_X0$moment[, j, drop = FALSE]),
    tolerance = 1e-12
  )
})

test_that("print.rdistnn_plqr_result runs without error", {
  d <- .make_plqr_data(n = 35L, p = 2L, p_d = 2L, m = 1L, tau = 0.5)
  out <- plqr_dnn(
    X = d$X, D = d$D, g = d$g, mD = d$mD, Y = d$Y,
    theta_eval = d$theta_eval, tau = d$tau, X0 = d$X0,
    k = 1L, s = 2L, se = FALSE
  )
  txt <- capture.output(print(out))
  expect_true(any(grepl("<rdistnn_plqr_result>", txt, fixed = TRUE)))
})

test_that("se = FALSE bypasses delete-1 jackknife feasibility check for plqr_dnn", {
  d <- .make_plqr_data(n = 10L, p = 2L, p_d = 1L, m = 2L, tau = 0.5)
  expect_no_error(
    plqr_dnn(
      X = d$X, D = d$D, g = d$g, mD = d$mD, Y = d$Y,
      theta_eval = d$theta_eval, tau = d$tau, X0 = d$X0,
      distance = "l2", k = 2L, s = nrow(d$X), method = "exact",
      se = FALSE, jackknife = "loo"
    )
  )
})

test_that("constant PLQR moment component implies zero jackknife SE", {
  d <- .make_plqr_data(n = 18L, p = 2L, p_d = 1L, m = 2L, tau = 0.5)
  # Force a degenerate case where D == mD everywhere, so the score is
  # identically 0 regardless of Y/g/theta, and the fast-jackknife SE on
  # `moment` must be 0.
  d$mD <- d$D
  out <- plqr_dnn(
    X = d$X, D = d$D, g = d$g, mD = d$mD, Y = d$Y,
    theta_eval = d$theta_eval, tau = d$tau, X0 = d$X0,
    distance = "l2", k = 2L, s = 4L, method = "exact",
    se = TRUE, jackknife = "fast"
  )
  expect_true(all(abs(out$moment) < 1e-12))
  expect_true(all(abs(out$moment_se) < 1e-12))
})
