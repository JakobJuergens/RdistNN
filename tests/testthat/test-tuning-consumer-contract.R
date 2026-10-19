.consumer_tuning <- function(s = 3L, grid = 2:4) new_tuning(
  name = "consumer_contract", eval = function(...) list(s = s),
  eval_grid = function(...) data.frame(s = grid))

test_that("experimental moment plots preserve the shared retained-tuning contract", {
  skip_if_not_installed("ggplot2")
  demo <- plqr_demo_data(n = 24L, p = 1L, p_d = 1L, seed = 9103L)
  d <- demo$data
  fit <- plqr_dnn(X = as.matrix(d["X1"]), Y = d$Y, D = as.matrix(d["D1"]),
    g = d$g_oracle, mD = as.matrix(d["mD1"]), theta_eval = demo$truth$theta_true,
    x0 = 0, tuning = .consumer_tuning(), se = FALSE)
  expect_setequal(plot_tuning(fit, se = FALSE)$data$s, 2:4)
  expect_setequal(plot_tuning(fit, tuning = .consumer_tuning(5L, c(3L, 5L)),
    se = FALSE)$data$s, c(3L, 5L))
})

test_that("manual-order sensitivity sweeps are independent of optional tuning returns", {
  skip_if_not_installed("ggplot2")
  X <- matrix(seq_len(100L))
  Y <- seq_len(100L)^2
  for (family in c("np", "cate", "ratio", "clate")) {
    fit <- function(extra) switch(family,
      np = np_dnn(X = X, Y = Y, x0 = 0, s = 10, se = FALSE, return = extra),
      cate = cate_dnn(X = X, pseudo_outcome = Y, x0 = 0, s = 10,
        se = FALSE, return = extra),
      ratio = affine_ratio_dnn(X = X, numerator = Y, denominator = rep(1, 100),
        x0 = 0, s = 10, se = FALSE, return = extra),
      clate = clate_dnn(X = X, Y = Y, D = rep(0:1, 50), A = rep(0:1, 50),
        e = rep(.5, 100), m_y0 = rep(0, 100), m_y1 = rep(0, 100),
        m_d0 = rep(0, 100), m_d1 = rep(0, 100),
        x0 = 0, s = 10, se = FALSE, return = extra))
    without_extra <- plot_tuning(fit(character()), se = FALSE)$data
    with_extra <- plot_tuning(fit("tuning"), se = FALSE)$data
    # n=100, kx=1: floor(c(.3,.4,.75)*100^(2/3)), plus the fitted order.
    expect_setequal(without_extra$s, c(6L, 8L, 16L, 10L))
    expect_setequal(with_extra$s, c(6L, 8L, 16L, 10L))
    expect_identical(with_extra, without_extra, info = family)
  }
})

test_that("point-only tuning grids reject s above n before evaluation", {
  tun <- new_tuning("core_feasibility", eval = function(...) list(s = 3L),
    eval_grid = function(...) data.frame(s = c(3L, 6L, 99L)))
  expect_warning(fit <- np_dnn(X = matrix(1:6), Y = 1:6, x0 = 0, tuning = tun,
    tuning_mode = "grid", se = FALSE), "Skipped infeasible")
  expect_identical(fit$tuning_grid$candidates$feasible, c(TRUE, TRUE, FALSE))
  expect_identical(fit$tuning_grid$candidates$status, c("ok", "ok", "infeasible"))
  expect_match(fit$tuning_grid$candidates$reason[3], "s > n")
})

test_that("scalar sensitivity honors tuning and a single uncertainty policy", {
  for (family in c("np", "cate")) for (fit_se in c(FALSE, TRUE)) {
    fit <- if (family == "np") np_dnn(X = matrix(1:8), Y = 1:8, x0 = 0,
      tuning = .consumer_tuning(), se = fit_se) else cate_dnn(X = matrix(1:8),
      pseudo_outcome = 1:8, x0 = 0, tuning = .consumer_tuning(), se = fit_se)
    for (se in c(FALSE, TRUE)) {
      grid <- if (family == "np") .plot_tuning_grid_np(fit, .consumer_tuning(grid = c(3L, 5L)), se) else
        .plot_tuning_grid_cate(fit, .consumer_tuning(grid = c(3L, 5L)), se)
      expect_setequal(grid$candidates$s, c(3L, 5L))
      expect_true(all(is.na(grid$estimates$se) == !se))
      if (se) expect_equal(grid$estimates$se[grid$estimates$s == 3],
        np_dnn(X = matrix(1:8), Y = 1:8, x0 = 0, s = 3)$se, tolerance = 1e-14)
    }
    expect_error(plot_tuning(fit, se = "yes"), "se.*TRUE or FALSE")
  }
})

test_that("NP and CATE tidy extractors select stored candidate rows", {
  for (family in c("np", "cate")) {
    fit <- if (family == "np") np_dnn(X = matrix(1:8), Y = 1:8, x0 = 0,
      tuning = .consumer_tuning(), tuning_mode = "grid", se = FALSE) else
      cate_dnn(X = matrix(1:8), pseudo_outcome = 1:8, x0 = 0,
        tuning = .consumer_tuning(), tuning_mode = "grid", se = FALSE)
    rows <- plot_data_estimates(fit, candidates = "all")
    expect_setequal(rows$s, 2:4)
    expect_equal(rows$estimate[rows$s == 4], 1.8, tolerance = 1e-14)
  }
})

test_that("scalar plots revalidate point-only candidate feasibility for SEs", {
  tuning <- .consumer_tuning(8L, c(3L, 8L))
  for (family in c("np", "cate")) {
    f <- if (family == "np") np_dnn(X = matrix(1:8), Y = 1:8, x0 = 0,
      tuning = tuning, se = FALSE) else cate_dnn(X = matrix(1:8), pseudo_outcome = 1:8,
        x0 = 0, tuning = tuning, se = FALSE)
    grid <- suppressWarnings(if (family == "np") .plot_tuning_grid_np(f, NULL, TRUE) else
      .plot_tuning_grid_cate(f, NULL, TRUE))
    expect_false(grid$candidates$feasible[grid$candidates$s == 8])
    expect_identical(grid$candidates$status[grid$candidates$s == 8], "infeasible")
    expect_true(is.finite(grid$estimates$se[grid$estimates$s == 3]))
  }
})

test_that("explicit grids skip an inserted default with infeasible uncertainty", {
  for (family in c("np", "cate", "ratio")) {
    f <- switch(family,
      np = np_dnn(X = matrix(1:8), Y = 1:8, x0 = 0, s = 8, se = FALSE),
      cate = cate_dnn(X = matrix(1:8), pseudo_outcome = 1:8, x0 = 0, s = 8, se = FALSE),
      ratio = affine_ratio_dnn(X = matrix(1:8), numerator = 1:8,
        denominator = rep(1, 8), x0 = 0, s = 8, se = FALSE))
    p <- suppressWarnings(plot_tuning(f, tuning = .consumer_tuning(3L, 3L), se = TRUE))
    expect_setequal(p$data$s, 3L)
    expect_true(all(is.finite(p$data$se)))
  }
})
