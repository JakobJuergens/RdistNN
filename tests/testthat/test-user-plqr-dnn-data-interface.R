# PLQR matrix/data-frame equivalence and data-frame input contracts.
test_that("data-frame and matrix APIs agree for plqr_dnn", {
  d <- .make_plqr_data(n = 70L, p = 2L, p_d = 2L, m = 3L, tau = 0.5)
  nms <- names(d$df)

  out_df <- plqr_dnn(
    data = d$df,
    x_cols = d$x_cols,
    D_cols = d$D_cols,
    g_col = "g_oracle",
    mD_cols = d$mD_cols,
    y_col = "Y",
    theta_eval = d$theta_eval,
    tau = d$tau,
    X0 = d$X0,
    k = 2L, s = 3L, method = "exact", se = FALSE
  )
  out_df_idx <- plqr_dnn(
    data = d$df,
    x_cols = match(d$x_cols, nms),
    D_cols = match(d$D_cols, nms),
    g_col = match("g_oracle", nms),
    mD_cols = match(d$mD_cols, nms),
    y_col = match("Y", nms),
    theta_eval = d$theta_eval,
    tau = d$tau,
    X0 = d$X0,
    k = 2L, s = 3L, method = "exact", se = FALSE
  )
  out_m <- plqr_dnn(
    X = d$X, D = d$D, g = d$g, mD = d$mD, Y = d$Y,
    theta_eval = d$theta_eval, tau = d$tau,
    X0 = d$X0,
    k = 2L, s = 3L, method = "exact", se = FALSE
  )

  expect_equal(out_df$moment, out_m$moment, tolerance = 1e-12)
  expect_equal(out_df$moment, out_df_idx$moment, tolerance = 1e-12)
})

test_that("plqr_dnn interfaces agree for exact method with fast jackknife", {
  d <- .make_plqr_data(n = 54L, p = 2L, p_d = 2L, m = 3L, tau = 0.5)
  nms <- names(d$df)

  out_mat <- plqr_dnn(
    X = d$X, D = d$D, g = d$g, mD = d$mD, Y = d$Y,
    theta_eval = d$theta_eval, tau = d$tau, X0 = d$X0,
    method = "exact", k = 2L, s = 3L, se = TRUE, jackknife = "fast"
  )
  out_df_name <- plqr_dnn(
    data = d$df, x_cols = d$x_cols, D_cols = d$D_cols, g_col = "g_oracle",
    mD_cols = d$mD_cols, y_col = "Y",
    theta_eval = d$theta_eval, tau = d$tau, X0 = d$X0,
    method = "exact", k = 2L, s = 3L, se = TRUE, jackknife = "fast"
  )
  out_df_idx <- plqr_dnn(
    data = d$df,
    x_cols = match(d$x_cols, nms),
    D_cols = match(d$D_cols, nms),
    g_col = match("g_oracle", nms),
    mD_cols = match(d$mD_cols, nms),
    y_col = match("Y", nms),
    theta_eval = d$theta_eval, tau = d$tau, X0 = d$X0,
    method = "exact", k = 2L, s = 3L, se = TRUE, jackknife = "fast"
  )

  expect_equal(out_df_name$moment, out_mat$moment, tolerance = 1e-12)
  expect_equal(out_df_name$moment_se, out_mat$moment_se, tolerance = 1e-12)
  expect_equal(out_df_idx$moment, out_mat$moment, tolerance = 1e-12)
  expect_equal(out_df_idx$moment_se, out_mat$moment_se, tolerance = 1e-12)
})

test_that("plqr_dnn interfaces agree for exact method with loo jackknife", {
  d <- .make_plqr_data(n = 40L, p = 2L, p_d = 2L, m = 2L, tau = 0.5)

  out_mat <- plqr_dnn(
    X = d$X, D = d$D, g = d$g, mD = d$mD, Y = d$Y,
    theta_eval = d$theta_eval, tau = d$tau, X0 = d$X0,
    method = "exact", k = 1L, s = 2L, se = TRUE, jackknife = "loo"
  )
  out_df <- plqr_dnn(
    data = d$df, x_cols = d$x_cols, D_cols = d$D_cols, g_col = "g_oracle",
    mD_cols = d$mD_cols, y_col = "Y",
    theta_eval = d$theta_eval, tau = d$tau, X0 = d$X0,
    method = "exact", k = 1L, s = 2L, se = TRUE, jackknife = "loo"
  )

  expect_equal(out_df$moment, out_mat$moment, tolerance = 1e-12)
  expect_equal(out_df$moment_se, out_mat$moment_se, tolerance = 1e-12)
})

test_that("plqr_dnn interfaces agree for asymp method with fast jackknife", {
  d <- .make_plqr_data(n = 52L, p = 2L, p_d = 2L, m = 3L, tau = 0.5)
  nms <- names(d$df)

  out_mat <- plqr_dnn(
    X = d$X, D = d$D, g = d$g, mD = d$mD, Y = d$Y,
    theta_eval = d$theta_eval, tau = d$tau, X0 = d$X0,
    method = "asymp", k = 2L, s = 4L, se = TRUE, jackknife = "fast"
  )
  out_df_name <- plqr_dnn(
    data = d$df, x_cols = d$x_cols, D_cols = d$D_cols, g_col = "g_oracle",
    mD_cols = d$mD_cols, y_col = "Y",
    theta_eval = d$theta_eval, tau = d$tau, X0 = d$X0,
    method = "asymp", k = 2L, s = 4L, se = TRUE, jackknife = "fast"
  )
  out_df_idx <- plqr_dnn(
    data = d$df,
    x_cols = match(d$x_cols, nms),
    D_cols = match(d$D_cols, nms),
    g_col = match("g_oracle", nms),
    mD_cols = match(d$mD_cols, nms),
    y_col = match("Y", nms),
    theta_eval = d$theta_eval, tau = d$tau, X0 = d$X0,
    method = "asymp", k = 2L, s = 4L, se = TRUE, jackknife = "fast"
  )

  expect_equal(out_df_name$moment, out_mat$moment, tolerance = 1e-12)
  expect_equal(out_df_name$moment_se, out_mat$moment_se, tolerance = 1e-12)
  expect_equal(out_df_idx$moment, out_mat$moment, tolerance = 1e-12)
  expect_equal(out_df_idx$moment_se, out_mat$moment_se, tolerance = 1e-12)
})

test_that("plqr_dnn x0 path agrees across matrix and data interfaces", {
  d <- .make_plqr_data(n = 44L, p = 2L, p_d = 2L, m = 4L, tau = 0.5)
  j <- 3L

  out_mat <- plqr_dnn(
    X = d$X, D = d$D, g = d$g, mD = d$mD, Y = d$Y,
    theta_eval = d$theta_eval, tau = d$tau, x0 = d$X0[j, ],
    method = "exact", k = 2L, s = 3L, se = FALSE
  )
  out_df <- plqr_dnn(
    data = d$df, x_cols = d$x_cols, D_cols = d$D_cols, g_col = "g_oracle",
    mD_cols = d$mD_cols, y_col = "Y",
    theta_eval = d$theta_eval, tau = d$tau, x0 = d$X0[j, ],
    method = "exact", k = 2L, s = 3L, se = FALSE
  )

  expect_equal(out_df$moment, out_mat$moment, tolerance = 1e-12)
})

test_that("plqr_dnn return extras and meta are identical across interfaces", {
  d <- .make_plqr_data(n = 56L, p = 2L, p_d = 2L, m = 3L, tau = 0.5)

  out_mat <- plqr_dnn(
    X = d$X, D = d$D, g = d$g, mD = d$mD, Y = d$Y,
    theta_eval = d$theta_eval, tau = d$tau, X0 = d$X0,
    method = "exact", k = 3L, s = 4L, se = FALSE,
    return = c("W", "ranks", "ord", "dx", "score", "ess", "neighbors")
  )
  out_df <- plqr_dnn(
    data = d$df, x_cols = d$x_cols, D_cols = d$D_cols, g_col = "g_oracle",
    mD_cols = d$mD_cols, y_col = "Y",
    theta_eval = d$theta_eval, tau = d$tau, X0 = d$X0,
    method = "exact", k = 3L, s = 4L, se = FALSE,
    return = c("W", "ranks", "ord", "dx", "score", "ess", "neighbors")
  )

  expect_equal(out_df$moment, out_mat$moment, tolerance = 1e-12)
  expect_equal(out_df$W, out_mat$W, tolerance = 1e-12)
  expect_equal(out_df$ranks, out_mat$ranks)
  expect_equal(out_df$ord, out_mat$ord)
  expect_equal(out_df$dx, out_mat$dx, tolerance = 1e-12)
  expect_equal(out_df$score, out_mat$score, tolerance = 1e-12)
  expect_equal(out_df$ess, out_mat$ess, tolerance = 1e-12)
  expect_equal(out_df$neighbors, out_mat$neighbors)

  expect_identical(out_df$meta$n, out_mat$meta$n)
  expect_identical(out_df$meta$m, out_mat$meta$m)
  expect_identical(out_df$meta$p, out_mat$meta$p)
  expect_identical(out_df$meta$p_d, out_mat$meta$p_d)
  expect_identical(out_df$meta$k, out_mat$meta$k)
  expect_identical(out_df$meta$s, out_mat$meta$s)
  expect_identical(out_df$meta$method, out_mat$meta$method)
  expect_identical(out_df$meta$distance, out_mat$meta$distance)
})

test_that("plqr_dnn preserves interface equivalence with distance spec object", {
  d <- .make_plqr_data(n = 46L, p = 2L, p_d = 2L, m = 3L, tau = 0.5)
  dist_spec <- distance_l2()

  out_mat <- plqr_dnn(
    X = d$X, D = d$D, g = d$g, mD = d$mD, Y = d$Y,
    theta_eval = d$theta_eval, tau = d$tau, X0 = d$X0,
    distance = dist_spec, method = "exact", k = 2L, s = 3L, se = FALSE
  )
  out_df <- plqr_dnn(
    data = d$df, x_cols = d$x_cols, D_cols = d$D_cols, g_col = "g_oracle",
    mD_cols = d$mD_cols, y_col = "Y",
    theta_eval = d$theta_eval, tau = d$tau, X0 = d$X0,
    distance = dist_spec, method = "exact", k = 2L, s = 3L, se = FALSE
  )

  expect_equal(out_df$moment, out_mat$moment, tolerance = 1e-12)
})

test_that("plqr_dnn tibble input matches data.frame input", {
  skip_if_not_installed("tibble")
  d <- .make_plqr_data(n = 42L, p = 2L, p_d = 2L, m = 3L, tau = 0.5)

  out_df <- plqr_dnn(
    data = d$df, x_cols = d$x_cols, D_cols = d$D_cols, g_col = "g_oracle",
    mD_cols = d$mD_cols, y_col = "Y",
    theta_eval = d$theta_eval, tau = d$tau, X0 = d$X0,
    method = "exact", k = 2L, s = 3L, se = FALSE
  )
  out_tbl <- plqr_dnn(
    data = tibble::as_tibble(d$df), x_cols = d$x_cols, D_cols = d$D_cols, g_col = "g_oracle",
    mD_cols = d$mD_cols, y_col = "Y",
    theta_eval = d$theta_eval, tau = d$tau, X0 = d$X0,
    method = "exact", k = 2L, s = 3L, se = FALSE
  )

  expect_equal(out_tbl$moment, out_df$moment, tolerance = 1e-12)
})

test_that("plqr_dnn rejects ambiguous data-frame column-name schemas", {
  d <- .make_plqr_data(n = 42L, p = 2L, p_d = 2L, m = 2L, tau = 0.5)
  nms <- names(d$df)
  df_dup <- d$df
  names(df_dup)[match(d$D_cols[2L], nms)] <- d$D_cols[1L]

  expect_error(
    plqr_dnn(
      data = df_dup,
      x_cols = match(d$x_cols, nms),
      D_cols = match(d$D_cols, nms),
      g_col = match("g_oracle", nms),
      mD_cols = match(d$mD_cols, nms),
      y_col = match("Y", nms),
      theta_eval = d$theta_eval,
      tau = d$tau,
      X0 = d$X0,
      k = 1L, s = 2L, se = FALSE
    ),
    "unique column names",
    ignore.case = TRUE
  )
})

test_that("plqr_dnn interface equivalence holds in p=1, p_d=1, m=1 edge case", {
  d <- .make_plqr_data(n = 38L, p = 1L, p_d = 1L, m = 1L, tau = 0.5)

  out_mat <- plqr_dnn(
    X = d$X, D = d$D, g = d$g, mD = d$mD, Y = d$Y,
    theta_eval = d$theta_eval, tau = d$tau, X0 = d$X0,
    method = "exact", k = 1L, s = 2L, se = FALSE,
    return = c("W", "score")
  )
  out_df <- plqr_dnn(
    data = d$df, x_cols = d$x_cols, D_cols = d$D_cols, g_col = "g_oracle",
    mD_cols = d$mD_cols, y_col = "Y",
    theta_eval = d$theta_eval, tau = d$tau, X0 = d$X0,
    method = "exact", k = 1L, s = 2L, se = FALSE,
    return = c("W", "score")
  )

  expect_equal(out_df$moment, out_mat$moment, tolerance = 1e-12)
  expect_equal(out_df$W, out_mat$W, tolerance = 1e-12)
  expect_equal(out_df$score, out_mat$score, tolerance = 1e-12)
})

test_that("data-frame API distinguishes NA vs non-finite x_cols index diagnostics for plqr_dnn", {
  d <- .make_plqr_data(n = 25L, p = 2L, p_d = 1L, m = 2L, tau = 0.5)

  expect_error(
    plqr_dnn(
      data = d$df, x_cols = c(1L, NA_integer_), D_cols = d$D_cols,
      g_col = "g_oracle", mD_cols = d$mD_cols, y_col = "Y",
      theta_eval = d$theta_eval, tau = d$tau, X0 = d$X0, k = 2L, s = 3L, se = FALSE
    ),
    regexp = "x_cols.*NA", ignore.case = TRUE
  )

  expect_error(
    plqr_dnn(
      data = d$df, x_cols = c(1L, Inf), D_cols = d$D_cols,
      g_col = "g_oracle", mD_cols = d$mD_cols, y_col = "Y",
      theta_eval = d$theta_eval, tau = d$tau, X0 = d$X0, k = 2L, s = 3L, se = FALSE
    ),
    regexp = "x_cols.*non-finite", ignore.case = TRUE
  )
})

test_that("data-frame API distinguishes NA vs non-finite covariate diagnostics for plqr_dnn", {
  d <- .make_plqr_data(n = 25L, p = 2L, p_d = 1L, m = 2L, tau = 0.5)
  base <- function(df) {
    tryCatch(
      plqr_dnn(
        data = df, x_cols = d$x_cols, D_cols = d$D_cols,
        g_col = "g_oracle", mD_cols = d$mD_cols, y_col = "Y",
        theta_eval = d$theta_eval, tau = d$tau, X0 = d$X0, k = 2L, s = 3L, se = FALSE
      ),
      error = conditionMessage
    )
  }

  df_na <- d$df
  df_na[[d$x_cols[[1L]]]][2L] <- NA_real_
  df_inf <- d$df
  df_inf[[d$x_cols[[1L]]]][2L] <- Inf

  msg_na <- base(df_na)
  msg_inf <- base(df_inf)

  expect_match(msg_na, "X.*NA", ignore.case = TRUE)
  expect_match(msg_inf, "X.*non-finite", ignore.case = TRUE)
  expect_false(identical(msg_na, msg_inf))
})

test_that("data-frame API reports y_col as missing when omitted", {
  d <- .make_plqr_data(n = 20L, p = 2L, p_d = 1L, m = 2L, tau = 0.5)
  expect_error(
    plqr_dnn(
      data = d$df, x_cols = d$x_cols, D_cols = d$D_cols,
      g_col = "g_oracle", mD_cols = d$mD_cols,
      theta_eval = d$theta_eval, tau = d$tau, X0 = d$X0, k = 2L, s = 3L, se = FALSE
    ),
    regexp = "missing: y_col", ignore.case = TRUE
  )
})
