# tests/testthat/test-score-plqr.R
# -------------------------------------------------------------------
# PLQR component score — PURE per-observation score only (no weights/USTAT)
# For any subset S, the kernel returns sum_{i in S} psi_i where
#   psi_i = (tau - 1{Y_i <= D_i' theta + g(U_i)}) * (D_{i,j} - mD_{i,j}).
# This suite matches the implementation in R/score-plqr.R.
# -------------------------------------------------------------------

testthat::local_edition(3)

# Score scenarios belong to this suite and use the public demo generator.
# Keep only the four designs actually consumed below, with their original seeds.
.get_plqr_fx <- function(n = 120L, p = 2L, p_d = 2L) {
  make_case <- function(design, sigma, seed) {
    demo <- plqr_demo_data(n, p, p_d, tau = 0.5, sigma = sigma,
                           design = design, seed = seed)
    list(data = demo$data, x0 = demo$x0,
      x_cols = paste0("X", seq_len(p)), d_cols = paste0("D", seq_len(p_d)),
      theta_true = demo$truth$theta_true, tau = 0.5,
      cols = list(Y = "Y", g = "g_oracle", mD = paste0("mD", seq_len(p_d)),
                  rD = paste0("rD", seq_len(p_d))),
      subs_empty = matrix(integer(), ncol = 1L),
      subs_2col = cbind(seq_len(n), rev(seq_len(n))))
  }
  list(ok_linear = make_case("linear", 0.5, 1L),
       ok_sinusoid = make_case("sinusoidal", 0.4, 2L),
       ok_bump = make_case("bump", 0.6, 3L),
       ok_clean = make_case("linear", 1e-8, 4L))
}

psi_tau <- function(u, tau) tau - as.numeric(u <= 0)

# --- tests -----------------------------------------------------------

test_that("[factory] validates tau, theta length, j range, D/mD lengths, and s", {
  fx <- .get_plqr_fx(n = 80, p = 2, p_d = 2)
  scn <- fx$ok_linear

  x0 <- scn$x0
  x_cols <- scn$x_cols
  D_cols <- scn$d_cols
  mD_cols <- scn$cols$mD
  g_col <- scn$cols$g
  tau <- scn$tau

  # tau outside (0,1)
  expect_error(
    score_plqr_component(x0, x_cols, D_cols, g_col, mD_cols, y_col = "Y",
      theta = rep(0, length(D_cols)), j = 1L, tau = 1, s = 1
    ),
    "`tau` must be a numeric scalar in \\(0,1\\)\\.",
    ignore.case = FALSE
  )
  expect_error(
    score_plqr_component(x0, x_cols, D_cols, g_col, mD_cols, y_col = "Y",
      theta = rep(0, length(D_cols)), j = 1L, tau = -1e-3, s = 1
    ),
    "`tau` must be a numeric scalar in \\(0,1\\)\\.",
    ignore.case = FALSE
  )

  # s must be a positive integer
  expect_error(
    score_plqr_component(x0, x_cols, D_cols, g_col, mD_cols, y_col = "Y",
      theta = rep(0, length(D_cols)), j = 1L, tau = tau, s = 0
    ),
    "`s` must be a positive integer\\.",
    ignore.case = FALSE
  )


  # D_cols non-empty
  expect_error(
    score_plqr_component(x0, x_cols,
      D_cols = character(0), g_col, mD_cols, y_col = "Y",
      theta = numeric(0), j = 1L, tau = tau, s = 1
    ),
    "`D_cols` must be a non-empty vector of columns\\.",
    ignore.case = FALSE
  )

  # mD length mismatch
  expect_error(
    score_plqr_component(x0, x_cols, D_cols, g_col, mD_cols[-1], y_col = "Y",
      theta = rep(0, length(D_cols)), j = 1L, tau = tau, s = 1
    ),
    "`mD_cols` must have the same length and order as `D_cols`\\.",
    ignore.case = FALSE
  )

  # theta length mismatch
  expect_error(
    score_plqr_component(x0, x_cols, D_cols, g_col, mD_cols, y_col = "Y",
      theta = numeric(length(D_cols) + 1L), j = 1L, tau = tau, s = 1
    ),
    "theta.*length.*must equal length\\(D_cols\\)",
    ignore.case = TRUE
  )

  # j out of range
  expect_error(
    score_plqr_component(x0, x_cols, D_cols, g_col, mD_cols, y_col = "Y",
      theta = rep(0, length(D_cols)), j = length(D_cols) + 1L,
      tau = tau, s = 1
    ),
    "`j` must be an integer in 1..length\\(D_cols\\)",
    ignore.case = TRUE
  )
})

test_that("[order] reported order equals chosen s", {
  fx <- .get_plqr_fx(n = 60, p = 2, p_d = 2)
  scn <- fx$ok_sinusoid

  for (ss in c(1L, 3L, 7L)) {
    sc <- score_plqr_component(scn$x0, scn$x_cols, scn$d_cols, scn$cols$g, scn$cols$mD,
      y_col = scn$cols$Y, theta = scn$theta_true, j = 1L, tau = scn$tau, s = ss
    )
    ctx <- score_setup(sc, scn$data)
    expect_identical(sc$order(nrow(scn$data), scn$data, ctx), as.integer(ss))
  }
})

test_that("[eval] matches direct formula for singletons (character x_cols)", {
  fx <- .get_plqr_fx(n = 120, p = 2, p_d = 2)
  scn <- fx$ok_linear

  df <- scn$data
  sc <- score_plqr_component(scn$x0, scn$x_cols, scn$d_cols, scn$cols$g, scn$cols$mD,
    y_col = scn$cols$Y, theta = scn$theta_true, j = 1L, tau = scn$tau, s = 1
  )
  ctx <- score_setup(sc, df)

  subs1 <- matrix(seq_len(nrow(df)), ncol = 1L)
  vals <- score_eval_batch(sc, subs1, df, ctx)

  Dmat <- as.matrix(df[scn$d_cols])
  thr <- as.vector(Dmat %*% scn$theta_true) + df[[scn$cols$g]]
  u <- df[[scn$cols$Y]] - thr
  rDj <- df[[scn$cols$rD[1]]]
  gR <- psi_tau(u, scn$tau) * rDj

  expect_equal(vals, gR, tolerance = 1e-12)
  expect_equal(score_eval(sc, 17L, df, ctx), gR[17], tolerance = 1e-12)

  # Empty batch
  expect_identical(score_eval_batch(sc, scn$subs_empty, df, ctx), numeric(0L))
})

test_that("[eval] multi-column subsets sum g across subset rows", {
  fx <- .get_plqr_fx(n = 90, p = 3, p_d = 2)
  scn <- fx$ok_bump

  df <- scn$data
  sc <- score_plqr_component(scn$x0, seq_along(scn$x_cols), scn$d_cols, scn$cols$g, scn$cols$mD,
    y_col = scn$cols$Y, theta = rep(0.1, length(scn$d_cols)), j = 2L, tau = scn$tau, s = 2
  )
  ctx <- score_setup(sc, df)

  subs <- scn$subs_2col
  vals <- score_eval_batch(sc, subs, df, ctx)

  # Direct row sums using the g-definition
  Dmat <- as.matrix(df[scn$d_cols])
  thr <- as.vector(Dmat %*% rep(0.1, length(scn$d_cols))) + df[[scn$cols$g]]
  u <- df[[scn$cols$Y]] - thr
  rDj <- df[[scn$cols$rD[2]]]
  g <- psi_tau(u, scn$tau) * rDj
  gR <- rowSums(matrix(g[subs], nrow = nrow(subs), ncol = ncol(subs)))

  expect_equal(vals, gR, tolerance = 1e-12)
})

test_that("[inputs/setup] finiteness guards fire at setup", {
  fx <- .get_plqr_fx(n = 40, p = 2, p_d = 2)
  scn <- fx$ok_clean
  df <- scn$data

  # Y non-finite
  df_bad <- df
  df_bad[[scn$cols$Y]][3] <- NA_real_
  sc <- score_plqr_component(scn$x0, scn$x_cols, scn$d_cols, scn$cols$g, scn$cols$mD,
    y_col = scn$cols$Y, theta = scn$theta_true, j = 1L, tau = scn$tau, s = 1
  )
  expect_error(score_setup(sc, df_bad), "Y cannot contain `NA` values\\.", ignore.case = FALSE)

  # g non-finite
  df_bad <- df
  df_bad[[scn$cols$g]][2] <- Inf
  expect_error(score_setup(sc, df_bad), "g cannot contain non-finite values", ignore.case = FALSE)

  # D non-finite
  df_bad <- df
  df_bad[[scn$d_cols[1]]][5] <- NaN
  expect_error(score_setup(sc, df_bad), "D cannot contain non-finite values", ignore.case = FALSE)

  # mD non-finite (need to add mD columns to df for this check)
  # If the demo generator didn't include mD columns in df, emulate via join:
  if (!all(scn$cols$mD %in% names(df))) {
    for (j in seq_along(scn$cols$mD)) df[[scn$cols$mD[j]]] <- df[[scn$d_cols[j]]] - df[[scn$cols$rD[j]]]
  }
  df_bad <- df
  df_bad[[scn$cols$mD[2]]][7] <- -Inf
  expect_error(score_setup(sc, df_bad), "mD cannot contain non-finite values", ignore.case = FALSE)

  # X non-finite
  df_bad <- df
  df_bad[[scn$x_cols[1]]][4] <- Inf
  expect_error(score_setup(sc, df_bad), "X cannot contain non-finite values", ignore.case = FALSE)

  # X NA
  df_bad <- df
  df_bad[[scn$x_cols[1]]][6] <- NA_real_
  expect_error(score_setup(sc, df_bad), "X cannot contain `NA` values\\.", ignore.case = FALSE)
})

test_that("[ratio] where rD_ij != 0, psi_i / rD_ij ∈ {tau - 1, tau}", {
  fx <- .get_plqr_fx(n = 110, p = 2, p_d = 2)
  scn <- fx$ok_linear

  df <- scn$data
  sc <- score_plqr_component(scn$x0, scn$x_cols, scn$d_cols, scn$cols$g, scn$cols$mD,
    y_col = scn$cols$Y, theta = scn$theta_true, j = 1L, tau = scn$tau, s = 1
  )
  ctx <- score_setup(sc, df)
  vals <- score_eval_batch(sc, matrix(seq_len(nrow(df)), ncol = 1L), df, ctx)

  rDj <- df[[scn$cols$rD[1]]]
  mask <- rDj != 0
  ratio <- vals[mask] / rDj[mask]
  ok <- (abs(ratio - scn$tau) <= 1e-12) | (abs(ratio - (scn$tau - 1)) <= 1e-12)
  expect_true(all(ok))
})

test_that("[object] score object shape and flags", {
  fx <- .get_plqr_fx(n = 50, p = 2, p_d = 1)
  scn <- fx$ok_linear

  sc <- score_plqr_component(scn$x0, scn$x_cols, scn$d_cols, scn$cols$g, scn$cols$mD,
    y_col = scn$cols$Y, theta = scn$theta_true, j = 1L, tau = scn$tau, s = 5L
  )
  expect_true(inherits(sc, "rdistnn_ustat_score"))
  expect_type(sc$eval, "closure")
  expect_type(sc$eval_batch, "closure")
  expect_false(isTRUE(sc$randomized))
  expect_true(isTRUE(sc$symmetric))
})

test_that("[parity] compiled batch equals direct row sums", {
  .test_rng_scope({
    fx <- .get_plqr_fx(n = 70, p = 2, p_d = 2)
    scn <- fx$ok_linear

    df <- scn$data
    sc <- score_plqr_component(scn$x0, scn$x_cols, scn$d_cols, scn$cols$g, scn$cols$mD,
      y_col = scn$cols$Y, theta = scn$theta_true, j = 2L, tau = scn$tau, s = 3
    )
    ctx <- score_setup(sc, df)

    subs <- cbind(
      seq_len(nrow(df)),
      rev(seq_len(nrow(df))),
      sample(seq_len(nrow(df)))
    )
    withr::local_options(RdistNN.use_compiled = TRUE)
    resolver <- RdistNN:::.rdistnn_compiled_fun
    cpp <- resolver("rustat_subset_row_sums_g")
    expect_true(is.function(cpp))
    calls <- 0L
    testthat::local_mocked_bindings(
      .rdistnn_compiled_fun = function(name) {
        if (!identical(name, "rustat_subset_row_sums_g")) return(resolver(name))
        function(...) {
          calls <<- calls + 1L
          cpp(...)
        }
      }, .package = "RdistNN"
    )
    out_compiled <- score_eval_batch(sc, subs, df, ctx)
    expect_identical(calls, 1L)
    # Direct row sums of the score using ctx$score (already computed by setup)
    gR <- rowSums(matrix(ctx$score[subs], nrow = nrow(subs), ncol = ncol(subs)))
    expect_equal(out_compiled, gR, tolerance = 1e-12)
  }, seed = 123L)
})
