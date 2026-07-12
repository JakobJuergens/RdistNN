test_that("np_dnn tuning path matches manual s and exposes tuning payload", {
  demo <- np_demo_data(n = 80L, p = 2L, seed = 1L)
  df <- demo$data
  x_mat <- as.matrix(df[, c("X1", "X2")])
  y_vec <- df$Y
  x0 <- demo$x0

  tun_fixed4 <- as_tuning(function(...) 4L, name = "fixed_4")

  fit_manual <- np_dnn(
    X = x_mat, Y = y_vec, x0 = x0,
    method = "exact", k = 2L, s = 4L, se = TRUE, jackknife = "fast"
  )
  fit_tuned <- np_dnn(
    X = x_mat, Y = y_vec, x0 = x0,
    method = "exact", k = 2L, tuning = tun_fixed4,
    se = TRUE, jackknife = "fast",
    return = "tuning"
  )

  expect_equal(fit_tuned$mu, fit_manual$mu, tolerance = 1e-12)
  expect_equal(fit_tuned$se, fit_manual$se, tolerance = 1e-12)
  expect_equal(fit_tuned$meta$s, 4L)
  expect_true(is.list(fit_tuned$meta$tuning))
  expect_equal(fit_tuned$meta$tuning$source, "tuning")
  expect_true(is.list(fit_tuned$tuning$default))
  expect_s3_class(fit_tuned$tuning$grid$candidates, "data.frame")
  expect_true(all(c("lambda", "s", "feasible", "reason") %in% names(fit_tuned$tuning$grid$candidates)))

  fit_df <- np_dnn(
    data = df, x_cols = c("X1", "X2"), y_col = "Y", x0 = x0,
    method = "exact", k = 2L, tuning = tun_fixed4, se = FALSE
  )
  fit_mat <- np_dnn(
    X = x_mat, Y = y_vec, x0 = x0,
    method = "exact", k = 2L, tuning = tun_fixed4, se = FALSE
  )
  expect_equal(fit_df$mu, fit_mat$mu, tolerance = 1e-12)
})

test_that("explicit s and tuning conflict is rejected across DNN wrappers", {
  x_mat <- matrix(rnorm(80), nrow = 40L, ncol = 2L)
  x0 <- c(0, 0)
  y_vec <- rnorm(40)
  psi <- rnorm(40)
  tun <- as_tuning(function(...) 4L, name = "fixed_4_conflict")

  expect_error(
    np_dnn(X = x_mat, Y = y_vec, x0 = x0, k = 2L, s = 4L, tuning = tun, se = FALSE),
    "cannot provide both explicit `s` and `tuning`"
  )
  expect_error(
    cate_dnn(X = x_mat, x0 = x0, pseudo_outcome = psi, k = 2L, s = 4L, tuning = tun, se = FALSE),
    "cannot provide both explicit `s` and `tuning`"
  )

  demo_plqr <- plqr_demo_data(n = 80L, p = 2L, p_d = 2L, tau = 0.5, seed = 1L)
  expect_error(
    plqr_dnn(
      data = demo_plqr$data,
      x_cols = c("X1", "X2"),
      D_cols = c("D1", "D2"),
      g_col = "g_oracle",
      mD_cols = c("mD1", "mD2"),
      y_col = "Y",
      theta_eval = demo_plqr$truth$theta_true,
      tau = 0.5,
      x0 = demo_plqr$x0,
      k = 2L,
      s = 4L,
      tuning = tun,
      se = FALSE
    ),
    "cannot provide both explicit `s` and `tuning`"
  )
})

test_that("cate_dnn and plqr_dnn tuning path matches manual s", {
  set.seed(123)
  n <- 70L
  x_mat <- matrix(rnorm(n * 2L), nrow = n, ncol = 2L)
  x0 <- c(0, 0)
  psi <- rnorm(n)
  tun <- as_tuning(function(...) 4L, name = "fixed_4_match")

  cate_manual <- cate_dnn(
    X = x_mat, x0 = x0, pseudo_outcome = psi,
    method = "exact", k = 2L, s = 4L, se = FALSE
  )
  cate_tuned <- cate_dnn(
    X = x_mat, x0 = x0, pseudo_outcome = psi,
    method = "exact", k = 2L, tuning = tun, se = FALSE
  )
  expect_equal(cate_tuned$tau, cate_manual$tau, tolerance = 1e-12)
  expect_equal(cate_tuned$meta$s, 4L)
  expect_true(is.list(cate_tuned$meta$tuning))

  demo_plqr <- plqr_demo_data(n = 90L, p = 2L, p_d = 2L, tau = 0.5, seed = 2L)
  plqr_manual <- plqr_dnn(
    data = demo_plqr$data,
    x_cols = c("X1", "X2"),
    D_cols = c("D1", "D2"),
    g_col = "g_oracle",
    mD_cols = c("mD1", "mD2"),
    y_col = "Y",
    theta_eval = demo_plqr$truth$theta_true,
    tau = 0.5,
    x0 = demo_plqr$x0,
    method = "exact",
    k = 2L,
    s = 4L,
    se = FALSE
  )
  plqr_tuned <- plqr_dnn(
    data = demo_plqr$data,
    x_cols = c("X1", "X2"),
    D_cols = c("D1", "D2"),
    g_col = "g_oracle",
    mD_cols = c("mD1", "mD2"),
    y_col = "Y",
    theta_eval = demo_plqr$truth$theta_true,
    tau = 0.5,
    x0 = demo_plqr$x0,
    method = "exact",
    k = 2L,
    tuning = tun,
    se = FALSE
  )
  expect_equal(plqr_tuned$moment, plqr_manual$moment, tolerance = 1e-12)
  expect_equal(plqr_tuned$meta$s, 4L)
  expect_true(is.list(plqr_tuned$meta$tuning))
})

test_that("omitted s resolves to the paper default; explicit s keeps the manual payload", {
  set.seed(10)
  x_mat <- matrix(rnorm(80), nrow = 40L, ncol = 2L)
  y_vec <- rnorm(40L)

  # Omitted s and tuning: the paper undersmoothing rule is the default.
  # (The historical silent default s = 1 produced uniform weights and is
  # deliberately no longer reachable without an explicit `s = 1L`.)
  fit <- np_dnn(
    X = x_mat, Y = y_vec, x0 = c(0, 0),
    method = "exact", k = 1L, se = FALSE,
    return = "tuning"
  )

  s_paper <- as.integer(floor(0.4 * 40^(3 / 4)))
  expect_equal(fit$meta$s, s_paper)
  expect_true(is.list(fit$meta$tuning))
  expect_equal(fit$tuning$source, "paper_default")
  expect_equal(fit$tuning$default$s, s_paper)
  expect_equal(nrow(fit$tuning$grid$candidates), 4L)

  # Explicit s: manual payload, no tuning metadata attached to meta.
  fit_manual <- np_dnn(
    X = x_mat, Y = y_vec, x0 = c(0, 0),
    method = "exact", k = 1L, s = 5L, se = FALSE,
    return = "tuning"
  )

  expect_equal(fit_manual$meta$s, 5L)
  expect_null(fit_manual$meta$tuning)
  expect_equal(fit_manual$tuning$source, "manual")
  expect_equal(fit_manual$tuning$default$s, 5L)
  expect_equal(fit_manual$tuning$default$meta$name, "manual")
  expect_equal(fit_manual$tuning$grid$candidates$s, 5L)
  expect_true(isTRUE(fit_manual$tuning$grid$candidates$feasible[1L]))
})

test_that("tuning works via registered name and bare function in wrappers", {
  demo <- np_demo_data(n = 90L, p = 2L, seed = 2L)
  x_mat <- as.matrix(demo$data[, c("X1", "X2")])
  y_vec <- demo$data$Y
  x0 <- demo$x0

  nm <- paste0("np_fixed_", sample.int(1e6, 1L))
  register_tuning(nm, function(...) 4L)

  fit_named <- np_dnn(
    X = x_mat, Y = y_vec, x0 = x0,
    method = "exact", k = 2L, tuning = nm,
    se = FALSE, return = "tuning"
  )
  fit_fun <- np_dnn(
    X = x_mat, Y = y_vec, x0 = x0,
    method = "exact", k = 2L, tuning = function(...) 4L,
    se = FALSE, return = "tuning"
  )
  fit_manual <- np_dnn(
    X = x_mat, Y = y_vec, x0 = x0,
    method = "exact", k = 2L, s = 4L, se = FALSE
  )

  expect_equal(fit_named$mu, fit_manual$mu, tolerance = 1e-12)
  expect_equal(fit_fun$mu, fit_manual$mu, tolerance = 1e-12)
  expect_equal(fit_named$meta$s, 4L)
  expect_equal(fit_fun$meta$s, 4L)
  expect_true(is.na(fit_fun$tuning$grid$candidates$lambda[1L]))
  expect_equal(fit_fun$tuning$grid$candidates$s, 4L)
})

test_that("tuning metadata is attached even when `return` omits tuning payload", {
  demo <- np_demo_data(n = 100L, p = 3L, seed = 3L)
  x_mat <- as.matrix(demo$data[, c("X1", "X2", "X3")])
  y_vec <- demo$data$Y
  tun <- tuning_dnn_paper(lambda = 1.25)
  expected_s <- tuning_eval(tun, n = nrow(x_mat), kx = ncol(x_mat), k = 2L, mode = "default")$s

  fit <- np_dnn(
    X = x_mat, Y = y_vec, x0 = demo$x0,
    method = "exact", k = 2L, tuning = tun, se = FALSE
  )

  expect_equal(fit$meta$s, expected_s)
  expect_true(is.list(fit$meta$tuning))
  expect_equal(fit$meta$tuning$source, "tuning")
  expect_null(fit$tuning)
})

test_that("tuning infeasibility k > s errors are propagated across wrappers", {
  set.seed(21)
  x_mat <- matrix(rnorm(80), nrow = 40L, ncol = 2L)
  x0 <- c(0, 0)
  y_vec <- rnorm(40L)
  psi <- rnorm(40L)
  tun_small <- as_tuning(function(...) 2L, name = "too_small_s")

  expect_error(
    np_dnn(X = x_mat, Y = y_vec, x0 = x0, method = "exact", k = 3L, tuning = tun_small, se = FALSE),
    "k <= s"
  )
  expect_error(
    cate_dnn(X = x_mat, x0 = x0, pseudo_outcome = psi, method = "exact", k = 3L, tuning = tun_small, se = FALSE),
    "k <= s"
  )

  demo_plqr <- plqr_demo_data(n = 85L, p = 2L, p_d = 2L, tau = 0.5, seed = 8L)
  expect_error(
    plqr_dnn(
      data = demo_plqr$data,
      x_cols = c("X1", "X2"),
      D_cols = c("D1", "D2"),
      g_col = "g_oracle",
      mD_cols = c("mD1", "mD2"),
      y_col = "Y",
      theta_eval = demo_plqr$truth$theta_true,
      tau = 0.5,
      x0 = demo_plqr$x0,
      method = "exact",
      k = 3L,
      tuning = tun_small,
      se = FALSE
    ),
    "k <= s"
  )
})

test_that("jackknife feasibility failures from tuning (s >= n) error clearly", {
  set.seed(22)
  x_mat <- matrix(rnorm(80), nrow = 40L, ncol = 2L)
  x0 <- c(0, 0)
  y_vec <- rnorm(40L)
  psi <- rnorm(40L)
  tun_jk_bad <- as_tuning(function(n, ...) as.integer(n), name = "jk_bad")

  expect_error(
    np_dnn(
      X = x_mat, Y = y_vec, x0 = x0,
      method = "exact", k = 1L, tuning = tun_jk_bad,
      se = TRUE, jackknife = "fast"
    ),
    "s < n"
  )
  expect_error(
    cate_dnn(
      X = x_mat, x0 = x0, pseudo_outcome = psi,
      method = "exact", k = 1L, tuning = tun_jk_bad,
      se = TRUE, jackknife = "fast"
    ),
    "s < n"
  )

  demo_plqr <- plqr_demo_data(n = 80L, p = 2L, p_d = 2L, tau = 0.5, seed = 9L)
  expect_error(
    plqr_dnn(
      data = demo_plqr$data,
      x_cols = c("X1", "X2"),
      D_cols = c("D1", "D2"),
      g_col = "g_oracle",
      mD_cols = c("mD1", "mD2"),
      y_col = "Y",
      theta_eval = demo_plqr$truth$theta_true,
      tau = 0.5,
      x0 = demo_plqr$x0,
      method = "exact",
      k = 1L,
      tuning = tun_jk_bad,
      se = TRUE,
      jackknife = "fast"
    ),
    "s < n"
  )
})

test_that("cate_dnn and plqr_dnn expose stable tuning payload with `return = 'tuning'`", {
  set.seed(23)
  n <- 70L
  x_mat <- matrix(rnorm(n * 2L), nrow = n, ncol = 2L)
  psi <- rnorm(n)
  tun <- tuning_dnn_paper(lambda = 1.25)

  cate_fit <- cate_dnn(
    X = x_mat, x0 = c(0, 0), pseudo_outcome = psi,
    method = "exact", k = 2L, tuning = tun, se = FALSE,
    return = "tuning"
  )
  expect_true(is.list(cate_fit$meta$tuning))
  expect_equal(cate_fit$meta$tuning$source, "tuning")
  expect_s3_class(cate_fit$tuning$grid$candidates, "data.frame")
  expect_true(all(c("lambda", "s", "feasible", "reason") %in% names(cate_fit$tuning$grid$candidates)))

  demo_plqr <- plqr_demo_data(n = 95L, p = 2L, p_d = 2L, tau = 0.5, seed = 10L)
  plqr_fit <- plqr_dnn(
    data = demo_plqr$data,
    x_cols = c("X1", "X2"),
    D_cols = c("D1", "D2"),
    g_col = "g_oracle",
    mD_cols = c("mD1", "mD2"),
    y_col = "Y",
    theta_eval = demo_plqr$truth$theta_true,
    tau = 0.5,
    x0 = demo_plqr$x0,
    method = "exact",
    k = 2L,
    tuning = tun,
    se = FALSE,
    return = "tuning"
  )
  expect_true(is.list(plqr_fit$meta$tuning))
  expect_equal(plqr_fit$meta$tuning$source, "tuning")
  expect_s3_class(plqr_fit$tuning$grid$candidates, "data.frame")
  expect_true(all(c("lambda", "s", "feasible", "reason") %in% names(plqr_fit$tuning$grid$candidates)))
})

test_that("matrix/data-frame interface parity is preserved with tuning enabled", {
  set.seed(24)
  cate_demo <- cate_demo_data(n = 90L, p = 2L, design = "hetero", sigma = 0.3, seed = 24L)
  df <- cate_demo$data
  x_cols <- c("X1", "X2")
  x_mat <- as.matrix(df[, x_cols, drop = FALSE])
  x0_mat <- matrix(rnorm(6L), nrow = 3L, ncol = 2L)
  den <- pmax(df$e_oracle * (1 - df$e_oracle), 1e-8)
  psi <- (df$D * (df$Y - (df$m0_oracle + df$tau_oracle)) - (1 - df$D) * (df$Y - df$m0_oracle)) / den + df$tau_oracle

  tun_fixed4 <- as_tuning(function(...) 4L, name = "fixed4_iface")
  out_cate_mat <- cate_dnn(
    X = x_mat, X0 = x0_mat, pseudo_outcome = psi,
    method = "exact", k = 2L, tuning = tun_fixed4, se = FALSE
  )
  out_cate_df <- cate_dnn(
    data = df, x_cols = x_cols, X0 = x0_mat, pseudo_outcome = psi,
    method = "exact", k = 2L, tuning = tun_fixed4, se = FALSE
  )
  expect_equal(out_cate_df$tau, out_cate_mat$tau, tolerance = 1e-12)

  plqr_demo <- plqr_demo_data(n = 90L, p = 2L, p_d = 2L, tau = 0.5, seed = 11L)
  out_plqr_mat <- plqr_dnn(
    X = as.matrix(plqr_demo$data[, c("X1", "X2")]),
    D = as.matrix(plqr_demo$data[, c("D1", "D2")]),
    g = plqr_demo$data$g_oracle,
    mD = as.matrix(plqr_demo$data[, c("mD1", "mD2")]),
    Y = plqr_demo$data$Y,
    theta_eval = plqr_demo$truth$theta_true,
    tau = 0.5,
    x0 = plqr_demo$x0,
    method = "exact",
    k = 2L,
    tuning = tun_fixed4,
    se = FALSE
  )
  out_plqr_df <- plqr_dnn(
    data = plqr_demo$data,
    x_cols = c("X1", "X2"),
    D_cols = c("D1", "D2"),
    g_col = "g_oracle",
    mD_cols = c("mD1", "mD2"),
    y_col = "Y",
    theta_eval = plqr_demo$truth$theta_true,
    tau = 0.5,
    x0 = plqr_demo$x0,
    method = "exact",
    k = 2L,
    tuning = tun_fixed4,
    se = FALSE
  )
  expect_equal(out_plqr_df$moment, out_plqr_mat$moment, tolerance = 1e-12)
})
