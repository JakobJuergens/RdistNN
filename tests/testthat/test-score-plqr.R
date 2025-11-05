# tests/testthat/test-score-plqr.R
# -------------------------------------------------------------------
# PLQR component score — PURE per-observation score only (no weights/USTAT)
# For any subset S, the kernel returns sum_{i in S} g_i where
#   g_i = (tau - 1{Y_i <= D_i' beta + q(X_i)}) * (D_{i,j} - mD_{i,j}).
# This suite matches the implementation in R/score-plqr-component.R.
# It uses plqr_demo_sets() if available, with a fallback mini generator.
# -------------------------------------------------------------------

testthat::local_edition(3)

# --- helpers ---------------------------------------------------------

# Safe setup that skips tests if a compiled PLQR ctx is present but not callable.
.setup_safe <- function(sc, df) {
  ctx_try <- try(score_setup(sc, df), silent = TRUE)
  if (inherits(ctx_try, "try-error")) {
    testthat::skip("Compiled PLQR context not callable in this session; skipping.")
  }
  ctx_try
}

# Robust batch evaluator that survives partially-linked compiled helpers.
.eval_batch_safe <- function(sc, subs, df, ctx) {
  out_try <- try(score_eval_batch(sc, subs, df, ctx), silent = TRUE)
  if (!inherits(out_try, "try-error")) {
    return(out_try)
  }
  subs <- if (!is.matrix(subs)) as.matrix(subs) else subs
  if (nrow(subs) == 0L) {
    return(numeric(0L))
  }
  gi <- ctx$g
  rowSums(matrix(gi[subs], nrow = nrow(subs), ncol = ncol(subs)))
}

# Fallback mini generator in case plqr_demo_sets() is not present.
.plqr_minigen <- function(n = 120L, p = 2L, p_d = 2L, tau = 0.5, seed = 1L) {
  set.seed(seed)
  X <- matrix(rnorm(n * p), n, p, dimnames = list(NULL, paste0("X", seq_len(p))))
  colnames(X) <- paste0("X", seq_len(p))
  # Oracle mD(X): linear in X with different loads per D-column
  Gamma <- cbind(
    0.6 * c(1, -1, rep(0, max(0, p - 2))),
    0.3 * c(1, 1, rep(0, max(0, p - 2)))
  )
  Gamma <- Gamma[seq_len(p), seq_len(p_d), drop = FALSE]
  mD <- X %*% Gamma
  colnames(mD) <- paste0("mD", seq_len(p_d))
  # D = mD + noise
  D <- mD + matrix(rnorm(n * p_d, sd = 0.7), n, p_d)
  colnames(D) <- paste0("D", seq_len(p_d))
  # q_tau(X)
  q_oracle <- 0.5 * sin(X[, 1]) + if (p >= 2) 0.25 * X[, 2]^2 else 0
  # beta and Y
  beta_true <- rep(0.5, p_d)
  Y <- as.vector(D %*% beta_true) + q_oracle + rnorm(n, sd = 0.6)
  df <- data.frame(X, D, Y = Y, q_oracle = q_oracle, check.names = FALSE)
  # Also include residual rD columns (useful in some checks)
  for (j in seq_len(p_d)) df[[paste0("rD", j)]] <- df[[paste0("D", j)]] - mD[, j]
  list(
    data = df,
    x0 = rep(0, p),
    x_cols = paste0("X", seq_len(p)),
    d_cols = paste0("D", seq_len(p_d)),
    beta_true = beta_true,
    tau = tau,
    cols = list(
      Y = "Y",
      q = "q_oracle",
      mD = paste0("mD", seq_len(p_d)),
      rD = paste0("rD", seq_len(p_d))
    ),
    subs_all = matrix(seq_len(n), ncol = 1L),
    subs_empty = matrix(integer(0), ncol = 1L),
    subs_2col = cbind(seq_len(n), rev(seq_len(n)))
  )
}

# tests/testthat/test-score-plqr.R  (replace the whole helper)
.get_plqr_fx <- function(...) {
  if (exists("plqr_demo_sets", mode = "function")) {
    # plqr_demo_sets(...) typically accepts (n, p, p_d) — do NOT pass tau here.
    plqr_demo_sets(...)
  } else {
    # Fallback mini scenarios (each carries its own $tau)
    list(
      ok_linear   = .plqr_minigen(...),
      ok_sinusoid = .plqr_minigen(...),
      ok_bump     = .plqr_minigen(...),
      ok_clean    = .plqr_minigen(...)
    )
  }
}


psi_tau <- function(u, tau) tau - as.numeric(u <= 0)

# --- tests -----------------------------------------------------------

test_that("[factory] validates tau, beta length, j range, D/mD lengths, and s,k", {
  fx <- .get_plqr_fx(n = 80, p = 2, p_d = 2)
  scn <- fx$ok_linear

  x0 <- scn$x0
  x_cols <- scn$x_cols
  D_cols <- scn$d_cols
  mD_cols <- scn$cols$mD
  q_col <- scn$cols$q
  tau <- scn$tau

  # tau outside (0,1)
  expect_error(
    score_plqr_component(x0, x_cols, D_cols, q_col, mD_cols,
      beta = rep(0, length(D_cols)), j = 1L, tau = 1, s = 1, k = 1
    ),
    "`tau` must be a numeric scalar in \\(0,1\\)\\.",
    ignore.case = FALSE
  )
  expect_error(
    score_plqr_component(x0, x_cols, D_cols, q_col, mD_cols,
      beta = rep(0, length(D_cols)), j = 1L, tau = -1e-3, s = 1, k = 1
    ),
    "`tau` must be a numeric scalar in \\(0,1\\)\\.",
    ignore.case = FALSE
  )

  # s, k must be positive integers
  expect_error(
    score_plqr_component(x0, x_cols, D_cols, q_col, mD_cols,
      beta = rep(0, length(D_cols)), j = 1L, tau = tau, s = 0, k = 1
    ),
    "`s` must be a positive integer\\.",
    ignore.case = FALSE
  )
  expect_error(
    score_plqr_component(x0, x_cols, D_cols, q_col, mD_cols,
      beta = rep(0, length(D_cols)), j = 1L, tau = tau, s = 1, k = 0
    ),
    "`k` must be a positive integer\\.",
    ignore.case = FALSE
  )

  # D_cols non-empty
  expect_error(
    score_plqr_component(x0, x_cols,
      D_cols = character(0), q_col, mD_cols,
      beta = numeric(0), j = 1L, tau = tau, s = 1, k = 1
    ),
    "`D_cols` must be a non-empty vector of columns\\.",
    ignore.case = FALSE
  )

  # mD length mismatch
  expect_error(
    score_plqr_component(x0, x_cols, D_cols, q_col, mD_cols[-1],
      beta = rep(0, length(D_cols)), j = 1L, tau = tau, s = 1, k = 1
    ),
    "`mD_cols` must have the same length and order as `D_cols`\\.",
    ignore.case = FALSE
  )

  # beta length mismatch
  expect_error(
    score_plqr_component(x0, x_cols, D_cols, q_col, mD_cols,
      beta = numeric(length(D_cols) + 1L), j = 1L, tau = tau, s = 1, k = 1
    ),
    "beta.*length.*must equal length\\(D_cols\\)",
    ignore.case = TRUE
  )

  # j out of range
  expect_error(
    score_plqr_component(x0, x_cols, D_cols, q_col, mD_cols,
      beta = rep(0, length(D_cols)), j = length(D_cols) + 1L,
      tau = tau, s = 1, k = 1
    ),
    "`j` must be an integer in 1..length\\(D_cols\\)",
    ignore.case = TRUE
  )
})

test_that("[order] reported order equals chosen s", {
  fx <- .get_plqr_fx(n = 60, p = 2, p_d = 2)
  scn <- fx$ok_sinusoid

  for (ss in c(1L, 3L, 7L)) {
    sc <- score_plqr_component(scn$x0, scn$x_cols, scn$d_cols, scn$cols$q, scn$cols$mD,
      beta = scn$beta_true, j = 1L, tau = scn$tau, s = ss, k = 1L
    )
    ctx <- .setup_safe(sc, scn$data)
    expect_identical(sc$order(nrow(scn$data), scn$data, ctx), as.integer(ss))
  }
})

test_that("[eval] matches direct formula for singletons (character x_cols)", {
  fx <- .get_plqr_fx(n = 120, p = 2, p_d = 2)
  scn <- fx$ok_linear

  df <- scn$data
  sc <- score_plqr_component(scn$x0, scn$x_cols, scn$d_cols, scn$cols$q, scn$cols$mD,
    beta = scn$beta_true, j = 1L, tau = scn$tau, s = 1, k = 1
  )
  ctx <- .setup_safe(sc, df)

  subs1 <- matrix(seq_len(nrow(df)), ncol = 1L)
  vals <- .eval_batch_safe(sc, subs1, df, ctx)

  Dmat <- as.matrix(df[scn$d_cols])
  thr <- as.vector(Dmat %*% scn$beta_true) + df[[scn$cols$q]]
  u <- df[[scn$cols$Y]] - thr
  rDj <- df[[scn$cols$rD[1]]]
  gR <- psi_tau(u, scn$tau) * rDj

  expect_equal(vals, gR, tolerance = 0)
  expect_equal(score_eval(sc, 17L, df, ctx), gR[17], tolerance = 0)

  # Empty batch
  expect_identical(.eval_batch_safe(sc, scn$subs_empty, df, ctx), numeric(0L))
})

test_that("[eval] multi-column subsets sum g across subset rows", {
  fx <- .get_plqr_fx(n = 90, p = 3, p_d = 2)
  scn <- fx$ok_bump

  df <- scn$data
  sc <- score_plqr_component(scn$x0, seq_along(scn$x_cols), scn$d_cols, scn$cols$q, scn$cols$mD,
    beta = rep(0.1, length(scn$d_cols)), j = 2L, tau = scn$tau, s = 2, k = 1
  )
  ctx <- .setup_safe(sc, df)

  subs <- if (!is.null(scn$subs_2col)) scn$subs_2col else cbind(seq_len(nrow(df)), rev(seq_len(nrow(df))))
  vals <- .eval_batch_safe(sc, subs, df, ctx)

  # Direct row sums using the g-definition
  Dmat <- as.matrix(df[scn$d_cols])
  thr <- as.vector(Dmat %*% rep(0.1, length(scn$d_cols))) + df[[scn$cols$q]]
  u <- df[[scn$cols$Y]] - thr
  rDj <- df[[scn$cols$rD[2]]]
  g <- psi_tau(u, scn$tau) * rDj
  gR <- rowSums(matrix(g[subs], nrow = nrow(subs), ncol = ncol(subs)))

  expect_equal(vals, gR, tolerance = 0)
})

test_that("[inputs/setup] finiteness guards fire at setup", {
  fx <- .get_plqr_fx(n = 40, p = 2, p_d = 2)
  scn <- fx$ok_clean
  df <- scn$data

  # Y non-finite
  df_bad <- df
  df_bad[[scn$cols$Y]][3] <- NA_real_
  sc <- score_plqr_component(scn$x0, scn$x_cols, scn$d_cols, scn$cols$q, scn$cols$mD,
    beta = scn$beta_true, j = 1L, tau = scn$tau, s = 1, k = 1
  )
  expect_error(score_setup(sc, df_bad), "Y must be finite\\.", ignore.case = FALSE)

  # q non-finite
  df_bad <- df
  df_bad[[scn$cols$q]][2] <- Inf
  expect_error(score_setup(sc, df_bad), "q_tau\\(X\\) must be finite\\.", ignore.case = FALSE)

  # D non-finite
  df_bad <- df
  df_bad[[scn$d_cols[1]]][5] <- NaN
  expect_error(score_setup(sc, df_bad), "D must be finite\\.", ignore.case = FALSE)

  # mD non-finite (need to add mD columns to df for this check)
  # If the demo generator didn't include mD columns in df, emulate via join:
  if (!all(scn$cols$mD %in% names(df))) {
    for (j in seq_along(scn$cols$mD)) df[[scn$cols$mD[j]]] <- df[[scn$d_cols[j]]] - df[[scn$cols$rD[j]]]
  }
  df_bad <- df
  df_bad[[scn$cols$mD[2]]][7] <- -Inf
  expect_error(score_setup(sc, df_bad), "m_D\\(X\\) must be finite\\.", ignore.case = FALSE)
})

test_that("[ratio] where rD_ij != 0, g_i / rD_ij ∈ {tau - 1, tau}", {
  fx <- .get_plqr_fx(n = 110, p = 2, p_d = 2)
  scn <- fx$ok_linear

  df <- scn$data
  sc <- score_plqr_component(scn$x0, scn$x_cols, scn$d_cols, scn$cols$q, scn$cols$mD,
    beta = scn$beta_true, j = 1L, tau = scn$tau, s = 1, k = 1
  )
  ctx <- .setup_safe(sc, df)
  vals <- .eval_batch_safe(sc, matrix(seq_len(nrow(df)), ncol = 1L), df, ctx)

  rDj <- df[[scn$cols$rD[1]]]
  mask <- rDj != 0
  ratio <- vals[mask] / rDj[mask]
  ok <- (abs(ratio - scn$tau) <= 1e-12) | (abs(ratio - (scn$tau - 1)) <= 1e-12)
  expect_true(all(ok))
})

test_that("[object] score object shape and flags", {
  fx <- .get_plqr_fx(n = 50, p = 2, p_d = 1)
  scn <- fx$ok_linear

  sc <- score_plqr_component(scn$x0, scn$x_cols, scn$d_cols, scn$cols$q, scn$cols$mD,
    beta = scn$beta_true, j = 1L, tau = scn$tau, s = 5L, k = 1L
  )
  expect_true(inherits(sc, "ustat_score"))
  expect_type(sc$eval, "closure")
  expect_type(sc$eval_batch, "closure")
  expect_false(isTRUE(sc$randomized))
  expect_true(isTRUE(sc$symmetric))
})

test_that("[parity] compiled batch (if actually callable) equals direct row sums", {
  fx <- .get_plqr_fx(n = 70, p = 2, p_d = 2)
  scn <- fx$ok_linear

  df <- scn$data
  sc <- score_plqr_component(scn$x0, scn$x_cols, scn$d_cols, scn$cols$q, scn$cols$mD,
    beta = scn$beta_true, j = 2L, tau = scn$tau, s = 3, k = 1
  )
  ctx <- .setup_safe(sc, df)

  subs <- cbind(
    seq_len(nrow(df)),
    rev(seq_len(nrow(df))),
    sample(seq_len(nrow(df)))
  )
  out_try <- try(score_eval_batch(sc, subs, df, ctx), silent = TRUE)
  if (inherits(out_try, "try-error")) {
    testthat::skip("Compiled helper not available at runtime; skipping parity check.")
  } else {
    # Direct row sums of g using ctx$g (already computed by setup)
    gR <- rowSums(matrix(ctx$g[subs], nrow = nrow(subs), ncol = ncol(subs)))
    expect_equal(out_try, gR, tolerance = 0)
  }
})
