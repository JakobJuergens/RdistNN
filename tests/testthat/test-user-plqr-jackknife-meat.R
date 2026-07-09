skip_if_not(
  all(vapply(
    c("plqr_dnn", "plqr_demo_data"),
    exists, logical(1),
    mode = "function"
  )),
  "plqr_dnn() or required helpers not found; skipping file"
)

.fit_plqr_meat <- function(d, theta_meat, jackknife = "fast", ...) {
  plqr_dnn(
    X = d$X, D = d$D, g = d$g, mD = d$mD, Y = d$Y,
    theta_eval = d$theta_eval, tau = d$tau, density = d$density,
    X0 = d$X0, method = "exact", k = 2L, s = 8L,
    se = TRUE, jackknife = jackknife, theta_meat = theta_meat,
    return = "W", ...
  )
}

test_that("jackknife meat returns full moment covariance with matching shapes", {
  d <- .make_plqr_inference_data(n = 70L, m = 2L)
  fit <- .fit_plqr_meat(d, theta_meat = "jackknife")
  p_d <- length(d$theta_eval)
  m <- nrow(d$X0)

  expect_true(is.list(fit$moment_vcov))
  expect_length(fit$moment_vcov, m)
  expect_named(fit$moment_vcov, paste0("target_", seq_len(m)))
  for (j in seq_len(m)) {
    v <- fit$moment_vcov[[j]]
    expect_true(is.matrix(v))
    expect_equal(dim(v), c(p_d, p_d))
    expect_equal(v, t(v), tolerance = 1e-12)
  }

  expect_equal(dim(fit$moment_se), c(p_d, m))
  for (j in seq_len(m)) {
    expect_equal(
      unname(fit$moment_se[, j]^2),
      unname(diag(fit$moment_vcov[[j]])),
      tolerance = 1e-12
    )
  }

  expect_identical(fit$meta$theta_meat_method, "jackknife")
  expect_identical(fit$meta$moment_se_method, "jackknife")
  expect_identical(fit$meta$theta_se_method, "direct_sandwich")
  expect_null(fit$score_covariance)
})

test_that("jackknife moment covariance matches a manual leave-one-out computation", {
  d <- .make_plqr_inference_data(n = 24L, m = 1L)
  fit <- .fit_plqr_meat(d, theta_meat = "jackknife", jackknife = "loo")

  n <- nrow(d$X)
  score <- (d$tau - as.numeric(d$Y <= as.vector(d$D %*% d$theta_eval) + d$g)) *
    (d$D - d$mD)

  loo <- matrix(NA_real_, nrow = n, ncol = length(d$theta_eval))
  for (i in seq_len(n)) {
    keep <- seq_len(n) != i
    wt_i <- weights_dnn(
      X = d$X[keep, , drop = FALSE], X0 = d$X0,
      distance = "l2", method = "exact", k = 2L, s = 8L,
      normalize = TRUE, return = "W"
    )
    loo[i, ] <- as.numeric(crossprod(score[keep, , drop = FALSE], wt_i$W))
  }
  centered <- sweep(loo, 2L, colMeans(loo), FUN = "-")
  manual_vcov <- ((n - 1) / n) * crossprod(centered)

  expect_equal(unname(fit$moment_vcov[[1L]]), unname(manual_vcov), tolerance = 1e-10)
  expect_equal(
    unname(fit$moment_se[, 1L]),
    unname(sqrt(diag(manual_vcov))),
    tolerance = 1e-10
  )
})

test_that("jackknife-meat coefficient covariance is the J-sandwiched moment covariance", {
  d <- .make_plqr_inference_data(n = 80L, m = 2L)
  fit <- .fit_plqr_meat(d, theta_meat = "jackknife")

  for (j in seq_len(nrow(d$X0))) {
    J_inv <- solve(fit$jacobian[[j]])
    manual <- J_inv %*% fit$moment_vcov[[j]] %*% t(J_inv)
    expect_equal(unname(fit$theta_vcov[[j]]), unname(manual), tolerance = 1e-10)
    expect_equal(
      unname(fit$theta_se[, j]),
      unname(sqrt(pmax(diag(manual), 0))),
      tolerance = 1e-10
    )
    expect_equal(
      unname(fit$theta_omega[[j]]),
      unname(manual * fit$neff[[j]]),
      tolerance = 1e-10
    )
  }
})

test_that("correlated scores produce nonzero off-diagonal moment covariance", {
  d <- .make_plqr_inference_data(n = 80L, m = 1L)
  # Make the second regressor a noisy copy of the first so the two score
  # components are strongly correlated.
  d$D[, 2L] <- d$D[, 1L] + 0.05 * d$D[, 2L]
  d$mD[, 2L] <- d$mD[, 1L] + 0.05 * d$mD[, 2L]
  fit <- .fit_plqr_meat(d, theta_meat = "jackknife")

  off_diag <- fit$moment_vcov[[1L]][1L, 2L]
  expect_true(is.finite(off_diag))
  expect_gt(abs(off_diag), 0)
})

test_that("analytic and jackknife meat share the localized moment and bread", {
  d <- .make_plqr_inference_data(n = 70L, m = 2L)
  fit_analytic <- .fit_plqr_meat(d, theta_meat = "analytic")
  fit_jk <- .fit_plqr_meat(d, theta_meat = "jackknife")

  expect_equal(fit_jk$moment, fit_analytic$moment, tolerance = 1e-12)
  expect_equal(fit_jk$theta, fit_analytic$theta, tolerance = 1e-12)
  expect_equal(fit_jk$jacobian, fit_analytic$jacobian, tolerance = 1e-12)
  expect_equal(fit_jk$neff, fit_analytic$neff, tolerance = 1e-12)

  expect_identical(fit_analytic$meta$theta_meat_method, "analytic")
  expect_null(fit_analytic$moment_vcov)
  expect_null(fit_analytic$moment_se)
  expect_true(is.list(fit_analytic$score_covariance))
})

test_that("jackknife-meat moment_se matches the no-density moment diagnostic", {
  d <- .make_plqr_inference_data(n = 60L, m = 2L)
  fit_jk <- .fit_plqr_meat(d, theta_meat = "jackknife", jackknife = "fast")
  fit_diag <- plqr_dnn(
    X = d$X, D = d$D, g = d$g, mD = d$mD, Y = d$Y,
    theta_eval = d$theta_eval, tau = d$tau,
    X0 = d$X0, method = "exact", k = 2L, s = 8L,
    se = TRUE, jackknife = "fast"
  )
  expect_equal(
    unname(fit_jk$moment_se),
    unname(fit_diag$moment_se),
    tolerance = 1e-12
  )
})

test_that("theta_meat = 'jackknife' validates its prerequisites", {
  d <- .make_plqr_inference_data(n = 40L, m = 1L)

  expect_error(
    plqr_dnn(
      X = d$X, D = d$D, g = d$g, mD = d$mD, Y = d$Y,
      theta_eval = d$theta_eval, tau = d$tau,
      X0 = d$X0, method = "exact", k = 2L, s = 6L,
      se = TRUE, theta_meat = "jackknife"
    ),
    "theta_meat.*density",
    ignore.case = TRUE
  )
  expect_error(
    plqr_dnn(
      X = d$X, D = d$D, g = d$g, mD = d$mD, Y = d$Y,
      theta_eval = d$theta_eval, tau = d$tau, density = d$density,
      X0 = d$X0, method = "exact", k = 2L, s = 6L,
      se = FALSE, theta_meat = "jackknife"
    ),
    "theta_meat.*se = TRUE",
    ignore.case = TRUE
  )
  expect_error(
    plqr_dnn(
      X = d$X, D = d$D, g = d$g, mD = d$mD, Y = d$Y,
      theta_eval = d$theta_eval, tau = d$tau, density = d$density,
      X0 = d$X0, method = "exact", k = 2L, s = 6L,
      se = TRUE, jackknife = "none", theta_meat = "jackknife"
    ),
    "theta_meat.*jackknife",
    ignore.case = TRUE
  )

  fit_none <- plqr_dnn(
    X = d$X, D = d$D, g = d$g, mD = d$mD, Y = d$Y,
    theta_eval = d$theta_eval, tau = d$tau,
    X0 = d$X0, method = "exact", k = 2L, s = 6L,
    se = FALSE
  )
  expect_identical(fit_none$meta$theta_meat_method, "none")
})

test_that("fast and loo jackknife meat agree for exact weights", {
  d <- .make_plqr_inference_data(n = 40L, m = 2L)
  fit_fast <- .fit_plqr_meat(d, theta_meat = "jackknife", jackknife = "fast")
  fit_loo <- .fit_plqr_meat(d, theta_meat = "jackknife", jackknife = "loo")

  expect_equal(fit_fast$moment_vcov, fit_loo$moment_vcov, tolerance = 1e-10)
  expect_equal(fit_fast$theta_vcov, fit_loo$theta_vcov, tolerance = 1e-10)
})
