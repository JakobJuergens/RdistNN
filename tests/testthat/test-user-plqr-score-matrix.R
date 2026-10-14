# Public PLQR score-matrix formula and validation.
test_that("plqr_score_matrix matches the orthogonal score formula", {
  d <- .make_plqr_data(n = 32L, p = 2L, p_d = 2L, m = 1L, tau = 0.35)
  score_names <- plqr_score_matrix(
    data = d$df, y_col = "Y", D_cols = d$D_cols,
    g_col = "g_oracle", mD_cols = d$mD_cols,
    theta = d$theta_eval, tau = d$tau
  )

  selected <- match(
    c("Y", d$D_cols, "g_oracle", d$mD_cols),
    names(d$df)
  )
  score_indices <- plqr_score_matrix(
    data = d$df, y_col = selected[1L], D_cols = selected[2:3],
    g_col = selected[4L], mD_cols = selected[5:6],
    theta = d$theta_eval, tau = d$tau
  )

  threshold <- as.vector(d$D %*% d$theta_eval) + d$g
  expected <- (d$tau - as.numeric(d$Y <= threshold)) * (d$D - d$mD)

  expect_equal(dim(score_names), c(32L, 2L))
  expect_equal(unname(score_names), unname(expected), tolerance = 1e-12)
  expect_equal(score_indices, score_names)
})

test_that("plqr_score_matrix validates its public data contract", {
  d <- .make_plqr_data(n = 30L, p = 2L, p_d = 2L, m = 1L, tau = 0.5)
  call_score <- function(data = d$df, D_cols = d$D_cols,
                         mD_cols = d$mD_cols, theta = d$theta_eval,
                         tau = d$tau) {
    plqr_score_matrix(
      data = data, y_col = "Y", D_cols = D_cols,
      g_col = "g_oracle", mD_cols = mD_cols,
      theta = theta, tau = tau
    )
  }

  expect_error(call_score(data = as.matrix(d$df)), "data")
  expect_error(call_score(tau = 1), "tau")
  expect_error(call_score(mD_cols = d$mD_cols[1L]), "mD_cols")
  expect_error(call_score(theta = d$theta_eval[1L]), "theta")

  bad_g <- d$df
  bad_g$g_oracle[2L] <- Inf
  expect_error(call_score(data = bad_g), "g_col")
})

test_that("plqr_score_matrix distinguishes NA vs non-finite theta diagnostics", {
  d <- .make_plqr_data(n = 30L, p = 2L, p_d = 2L, m = 1L, tau = 0.5)

  theta_na <- d$theta_eval
  theta_na[1] <- NA_real_
  expect_error(
    plqr_score_matrix(
      data = d$df, y_col = "Y", D_cols = d$D_cols, g_col = "g_oracle",
      mD_cols = d$mD_cols, theta = theta_na, tau = d$tau
    ),
    "theta.*NA",
    ignore.case = TRUE
  )

  theta_inf <- d$theta_eval
  theta_inf[1] <- Inf
  expect_error(
    plqr_score_matrix(
      data = d$df, y_col = "Y", D_cols = d$D_cols, g_col = "g_oracle",
      mD_cols = d$mD_cols, theta = theta_inf, tau = d$tau
    ),
    "theta.*non-finite",
    ignore.case = TRUE
  )
})
