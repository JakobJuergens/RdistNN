# PLQR matrix inputs, wrapper controls, and removed argument names.
test_that("plqr_dnn validation catches shape and theta_eval issues", {
  d <- .make_plqr_data(n = 30L, p = 2L, p_d = 2L, m = 2L, tau = 0.5)

  expect_error(
    plqr_dnn(
      X = d$X, D = d$D, g = d$g, mD = d$mD, Y = d$Y,
      theta_eval = d$theta_eval[-1], tau = d$tau, X0 = d$X0,
      k = 1L, s = 2L, se = FALSE
    ),
    "theta_eval.*ncol\\(D\\)",
    ignore.case = TRUE
  )

  expect_error(
    plqr_dnn(
      X = d$X, D = d$D, g = d$g, mD = d$mD[, 1, drop = FALSE], Y = d$Y,
      theta_eval = d$theta_eval, tau = d$tau, X0 = d$X0,
      k = 1L, s = 2L, se = FALSE
    ),
    "match ncol\\(D\\)|same number of columns",
    ignore.case = TRUE
  )

  expect_error(
    plqr_dnn(
      X = d$X, D = d$D, g = d$g, mD = d$mD, Y = d$Y,
      theta_eval = d$theta_eval, tau = d$tau, X0 = d$X0,
      k = 1L, s = nrow(d$X), se = TRUE, jackknife = "loo"
    ),
    "s < n",
    ignore.case = TRUE
  )
  expect_error(
    plqr_dnn(
      X = d$X, D = d$D, g = d$g, mD = d$mD, Y = d$Y,
      theta_eval = d$theta_eval, tau = 1, X0 = d$X0,
      k = 1L, s = 2L, se = FALSE
    ),
    "numeric scalar in \\(0,1\\)",
    ignore.case = TRUE
  )
  expect_error(
    plqr_dnn(
      X = d$X, D = d$D, g = d$g, mD = d$mD, Y = d$Y,
      theta_eval = d$theta_eval, tau = d$tau, X0 = d$X0, x0 = d$X0[1, ],
      k = 1L, s = 2L, se = FALSE
    ),
    "either `X0` or `x0`",
    ignore.case = TRUE
  )
  expect_error(
    plqr_dnn(
      X = d$X, D = d$D, g = d$g, mD = d$mD, Y = d$Y,
      theta_eval = d$theta_eval, tau = d$tau, x0 = c(NA_real_, 0),
      k = 1L, s = 2L, se = FALSE
    ),
    "x0.*NA",
    ignore.case = TRUE
  )
  bad_X0 <- d$X0
  bad_X0[1L, 1L] <- Inf
  expect_error(
    plqr_dnn(
      X = d$X, D = d$D, g = d$g, mD = d$mD, Y = d$Y,
      theta_eval = d$theta_eval, tau = d$tau, X0 = bad_X0,
      k = 1L, s = 2L, se = FALSE
    ),
    "X0.*non-finite",
    ignore.case = TRUE
  )
  expect_error(
    plqr_dnn(
      X = d$X, D = d$D, g = d$g, mD = d$mD, Y = d$Y,
      theta_eval = d$theta_eval, tau = d$tau,
      X0 = matrix(numeric(0), nrow = 0L, ncol = ncol(d$X)),
      k = 1L, s = 2L, se = FALSE
    ),
    "X0.*at least 1 row",
    ignore.case = TRUE
  )
  expect_error(
    plqr_dnn(
      data = d$df,
      x_cols = d$x_cols,
      D_cols = d$D_cols,
      theta_eval = d$theta_eval,
      tau = d$tau,
      X0 = d$X0,
      k = 1L, s = 2L, se = FALSE
    ),
    "x_cols.*D_cols.*g_col.*mD_cols.*y_col",
    ignore.case = TRUE
  )
  expect_error(
    plqr_dnn(
      data = d$df,
      x_cols = d$x_cols,
      D_cols = c(d$D_cols[1], "D_missing"),
      g_col = "g_oracle",
      mD_cols = d$mD_cols,
      y_col = "Y",
      theta_eval = d$theta_eval,
      tau = d$tau,
      X0 = d$X0,
      k = 1L, s = 2L, se = FALSE
    ),
    "specified by `D_cols`",
    ignore.case = TRUE
  )
  df_bad_mD <- d$df
  df_bad_mD[[d$mD_cols[1]]] <- as.character(df_bad_mD[[d$mD_cols[1]]])
  expect_error(
    plqr_dnn(
      data = df_bad_mD,
      x_cols = d$x_cols,
      D_cols = d$D_cols,
      g_col = "g_oracle",
      mD_cols = d$mD_cols,
      y_col = "Y",
      theta_eval = d$theta_eval,
      tau = d$tau,
      X0 = d$X0,
      k = 1L, s = 2L, se = FALSE
    ),
    "mD.*numeric",
    ignore.case = TRUE
  )
  # k and s validation at entry point
  expect_error(
    plqr_dnn(
      X = d$X, D = d$D, g = d$g, mD = d$mD, Y = d$Y,
      theta_eval = d$theta_eval, tau = d$tau, X0 = d$X0,
      k = 0L, s = 3L, se = FALSE
    ),
    "positive integer"
  )
  expect_error(
    plqr_dnn(
      X = d$X, D = d$D, g = d$g, mD = d$mD, Y = d$Y,
      theta_eval = d$theta_eval, tau = d$tau, X0 = d$X0,
      k = 3L, s = 2L, se = FALSE
    ),
    "k <= s"
  )
  expect_error(
    plqr_dnn(
      X = d$X, D = d$D, g = d$g, mD = d$mD, Y = d$Y,
      theta_eval = d$theta_eval, tau = d$tau, X0 = d$X0,
      k = 1.5, s = 3L, se = FALSE
    ),
    "integer scalar"
  )
  expect_error(
    plqr_dnn(
      X = d$X, D = d$D, g = d$g, mD = d$mD, Y = d$Y,
      theta_eval = d$theta_eval, tau = d$tau, X0 = d$X0,
      k = "foo", s = 3L, se = FALSE
    ),
    "integer scalar"
  )
  expect_error(
    plqr_dnn(
      X = d$X, D = d$D, g = d$g, mD = d$mD, Y = d$Y,
      theta_eval = d$theta_eval, tau = d$tau, X0 = d$X0,
      k = 1L, s = 3L, se = NA
    ),
    "se.*TRUE or FALSE"
  )
  expect_error(
    plqr_dnn(
      X = d$X, D = d$D, g = d$g, mD = d$mD, Y = d$Y,
      theta_eval = d$theta_eval, tau = d$tau, X0 = d$X0,
      k = 1L, s = 3L, se = "yes"
    ),
    "se.*TRUE or FALSE"
  )
  expect_error(
    plqr_dnn(
      X = d$X, D = d$D, g = d$g, mD = d$mD, Y = d$Y,
      theta_eval = d$theta_eval, tau = d$tau, X0 = d$X0,
      method = "asymp", k = 1L, s = 3L, se = FALSE,
      normalize = NA
    ),
    "normalize.*TRUE or FALSE"
  )
  expect_error(
    plqr_dnn(
      X = d$X, D = d$D, g = d$g, mD = d$mD, Y = d$Y,
      theta_eval = d$theta_eval, tau = d$tau, X0 = d$X0,
      method = "asymp", k = 1L, s = 3L, se = FALSE,
      normalize = "yes"
    ),
    "normalize.*TRUE or FALSE"
  )
})

test_that("plqr_dnn shape diagnostics include expected and actual values", {
  d <- .make_plqr_data(n = 28L, p = 2L, p_d = 2L, m = 1L, tau = 0.5)
  msg <- tryCatch(
    {
      plqr_dnn(
        X = d$X, D = d$D[-1, , drop = FALSE], g = d$g, mD = d$mD, Y = d$Y,
        theta_eval = d$theta_eval, tau = d$tau, X0 = d$X0,
        k = 1L, s = 2L, se = FALSE
      )
      ""
    },
    error = function(e) conditionMessage(e)
  )
  expect_match(msg, "Expected:", ignore.case = FALSE)
  expect_match(msg, "Got:", ignore.case = FALSE)
})

test_that("plqr_dnn matrix API distinguishes NA vs non-finite diagnostics", {
  d <- .make_plqr_data(n = 28L, p = 2L, p_d = 2L, m = 1L, tau = 0.5)

  D_na <- d$D
  D_na[1, 1] <- NA_real_
  expect_error(
    plqr_dnn(
      X = d$X, D = D_na, g = d$g, mD = d$mD, Y = d$Y,
      theta_eval = d$theta_eval, tau = d$tau, X0 = d$X0,
      k = 1L, s = 2L, se = FALSE
    ),
    "D.*NA",
    ignore.case = TRUE
  )

  D_inf <- d$D
  D_inf[1, 1] <- Inf
  expect_error(
    plqr_dnn(
      X = d$X, D = D_inf, g = d$g, mD = d$mD, Y = d$Y,
      theta_eval = d$theta_eval, tau = d$tau, X0 = d$X0,
      k = 1L, s = 2L, se = FALSE
    ),
    "D.*non-finite",
    ignore.case = TRUE
  )

  g_na <- d$g
  g_na[2] <- NA_real_
  expect_error(
    plqr_dnn(
      X = d$X, D = d$D, g = g_na, mD = d$mD, Y = d$Y,
      theta_eval = d$theta_eval, tau = d$tau, X0 = d$X0,
      k = 1L, s = 2L, se = FALSE
    ),
    "g.*NA",
    ignore.case = TRUE
  )

  mD_inf <- d$mD
  mD_inf[2, 1] <- Inf
  expect_error(
    plqr_dnn(
      X = d$X, D = d$D, g = d$g, mD = mD_inf, Y = d$Y,
      theta_eval = d$theta_eval, tau = d$tau, X0 = d$X0,
      k = 1L, s = 2L, se = FALSE
    ),
    "mD.*non-finite",
    ignore.case = TRUE
  )

  y_inf <- d$Y
  y_inf[3] <- Inf
  expect_error(
    plqr_dnn(
      X = d$X, D = d$D, g = d$g, mD = d$mD, Y = y_inf,
      theta_eval = d$theta_eval, tau = d$tau, X0 = d$X0,
      k = 1L, s = 2L, se = FALSE
    ),
    "Y.*non-finite",
    ignore.case = TRUE
  )
})

test_that("plqr_dnn warns on unknown return values", {
  d <- .make_plqr_data(n = 35L, p = 2L, p_d = 2L, m = 2L, tau = 0.5)
  expect_warning(
    plqr_dnn(
      X = d$X, D = d$D, g = d$g, mD = d$mD, Y = d$Y,
      theta_eval = d$theta_eval, tau = d$tau, X0 = d$X0,
      k = 1L, s = 2L, se = FALSE, return = c("W", "typo")
    ),
    "Unknown `return` value"
  )
  expect_no_warning(
    plqr_dnn(
      X = d$X, D = d$D, g = d$g, mD = d$mD, Y = d$Y,
      theta_eval = d$theta_eval, tau = d$tau, X0 = d$X0,
      k = 1L, s = 2L, se = FALSE,
      return = c("W", "ranks", "ord", "dx", "score", "ess", "neighbors")
    )
  )
})

test_that("X0 ncol mismatch vs X is rejected", {
  d <- .make_plqr_data(n = 30L, p = 2L, p_d = 2L, m = 2L, tau = 0.5)
  expect_error(
    plqr_dnn(
      X = d$X, D = d$D, g = d$g, mD = d$mD, Y = d$Y,
      theta_eval = d$theta_eval, tau = d$tau,
      X0 = matrix(0, 2, ncol(d$X) + 1),
      distance = "l2", k = 2L, s = 3L, se = FALSE
    ),
    regexp = "ncol\\(X0\\) .* ncol\\(X\\)", ignore.case = TRUE
  )
})

test_that("invalid jackknife argument errors clearly for plqr_dnn", {
  d <- .make_plqr_data(n = 20L, p = 2L, p_d = 1L, m = 3L, tau = 0.5)
  expect_error(
    plqr_dnn(
      X = d$X, D = d$D, g = d$g, mD = d$mD, Y = d$Y,
      theta_eval = d$theta_eval, tau = d$tau, X0 = d$X0,
      distance = "l2", k = 2L, s = 3L,
      se = TRUE, jackknife = "bad_option"
    ),
    regexp = "one of"
  )
})

test_that("matrix interface reports each missing D/g/mD/Y argument individually", {
  d <- .make_plqr_data(n = 15L, p = 2L, p_d = 1L, m = 2L, tau = 0.5)

  expect_error(
    plqr_dnn(X = d$X, g = d$g, mD = d$mD, Y = d$Y, theta_eval = d$theta_eval, tau = d$tau, X0 = d$X0, se = FALSE),
    regexp = "missing: D", ignore.case = TRUE
  )
  expect_error(
    plqr_dnn(X = d$X, D = d$D, mD = d$mD, Y = d$Y, theta_eval = d$theta_eval, tau = d$tau, X0 = d$X0, se = FALSE),
    regexp = "missing: g", ignore.case = TRUE
  )
  expect_error(
    plqr_dnn(X = d$X, D = d$D, g = d$g, Y = d$Y, theta_eval = d$theta_eval, tau = d$tau, X0 = d$X0, se = FALSE),
    regexp = "missing: mD", ignore.case = TRUE
  )
  expect_error(
    plqr_dnn(X = d$X, D = d$D, g = d$g, mD = d$mD, theta_eval = d$theta_eval, tau = d$tau, X0 = d$X0, se = FALSE),
    regexp = "missing: Y", ignore.case = TRUE
  )
})

test_that("matrix API rejects non-finite X for plqr_dnn", {
  d <- .make_plqr_data(n = 20L, p = 2L, p_d = 1L, m = 2L, tau = 0.5)
  X_bad <- d$X
  X_bad[3L, 1L] <- Inf
  expect_error(
    plqr_dnn(
      X = X_bad, D = d$D, g = d$g, mD = d$mD, Y = d$Y,
      theta_eval = d$theta_eval, tau = d$tau, X0 = d$X0, k = 2L, s = 3L, se = FALSE
    ),
    regexp = "X.*non-finite", ignore.case = TRUE
  )
})

test_that("removed legacy PLQR argument names fail instead of aliasing", {
  d <- .make_plqr_data(n = 20L, p = 2L, p_d = 1L, m = 1L, tau = 0.5)
  expect_error(
    plqr_dnn(
      X = d$X, D = d$D, q = d$g, mD = d$mD, Y = d$Y,
      beta = d$theta_eval, tau = d$tau, X0 = d$X0,
      k = 1L, s = 2L, se = FALSE
    ),
    "unused argument"
  )
  expect_error(
    plqr_dnn(
      data = d$df, x_cols = d$x_cols, D_cols = d$D_cols,
      q_col = "g_oracle", mD_cols = d$mD_cols, y_col = "Y",
      theta_eval = d$theta_eval, tau = d$tau, X0 = d$X0,
      k = 1L, s = 2L, se = FALSE
    ),
    "unused argument"
  )
})

test_that("plqr_dnn validates direct sandwich density inputs", {
  .test_rng_scope({
    d <- .make_plqr_inference_data(n = 50L, m = 1L)

    expect_error(
      plqr_dnn(
        X = d$X, D = d$D, g = d$g, mD = d$mD, Y = d$Y,
        theta_eval = d$theta_eval, tau = d$tau, density = d$density[-1L],
        X0 = d$X0, method = "exact", k = 1L, s = 6L, se = FALSE
      ),
      "density.*length",
      ignore.case = TRUE
    )

    density_zero <- d$density
    density_zero[1L] <- 0
    expect_error(
      plqr_dnn(
        X = d$X, D = d$D, g = d$g, mD = d$mD, Y = d$Y,
        theta_eval = d$theta_eval, tau = d$tau, density = density_zero,
        X0 = d$X0, method = "exact", k = 1L, s = 6L, se = FALSE
      ),
      "density.*positive",
      ignore.case = TRUE
    )

    expect_error(
      plqr_dnn(
        data = d$df,
        x_cols = d$x_cols,
        D_cols = d$D_cols,
        g_col = "g_oracle",
        mD_cols = d$mD_cols,
        y_col = "Y",
        density_col = "missing_density",
        theta_eval = d$theta_eval,
        tau = d$tau,
        X0 = d$X0,
        method = "exact", k = 1L, s = 6L, se = FALSE
      ),
      "density_col",
      ignore.case = TRUE
    )

    expect_error(
      plqr_dnn(
        X = d$X, D = d$D, g = d$g, mD = d$mD, Y = d$Y,
        theta_eval = d$theta_eval, tau = d$tau, density = d$density,
        X0 = d$X0, method = "exact", k = 1L,
        tuning_mode = "grid",
        se = TRUE
      ),
      "tuning_mode.*density",
      ignore.case = TRUE
    )

    expect_error(
      plqr_dnn(
        X = d$X, D = d$D, g = d$g, mD = d$mD, Y = d$Y,
        theta_eval = d$theta_eval, tau = d$tau, density = d$density,
        X0 = d$X0, method = "asymp", k = 1L, s = 6L,
        normalize = FALSE, se = TRUE
      ),
      "normalize.*direct PLQR coefficient inference",
      ignore.case = TRUE
    )
  }, seed = 818L)
})
