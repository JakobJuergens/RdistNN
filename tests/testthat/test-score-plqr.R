# tests/testthat/test-score-plqr.R
# -------------------------------------------------------------------
# PLQR component score — PURE per-observation score only (no weights/USTAT)
# We evaluate singleton subsets (s=1,k=1), so the score equals
#   g_i = (tau - 1{Y_i <= D_i' beta + q(X_i)}) * (D_ij - mD_ij).
# -------------------------------------------------------------------

testthat::local_edition(3)

# --- helpers ---------------------------------------------------------
plqr_demo <- function(n = 200L, px = 2L, tau = 0.5, seed = 1L) {
  set.seed(seed)
  X  <- matrix(runif(n * px), n, px, dimnames = list(NULL, paste0("X", seq_len(px))))
  u  <- rnorm(n)
  qx <- as.numeric(X[, 1] - 0.5 * (px >= 2) * X[, 2])
  Y  <- u + qx
  D1  <- rep(1, n)     # makes rDj = 1 - 0 = 1 in j=1 tests
  mD1 <- rep(0, n)
  data.frame(X, D1 = D1, mD1 = mD1, q = qx, Y = Y, check.names = FALSE)
}

.singletons <- function(n) matrix(seq_len(n), ncol = 1L)

.skip_if_no_factory <- function() {
  testthat::skip_if_not(
    exists("score_plqr_component", mode = "function"),
    "score_plqr_component() not found."
  )
}

.has_plqr_cpp <- function() {
  exists("rustat_plqr_g_all", mode = "function")
}

psi_tau <- function(u, tau) tau - as.numeric(u <= 0)

# --- tests -----------------------------------------------------------

# [factory] validates tau, beta length, j range, and shapes
test_that("[factory] validates tau, beta length, j range, and shapes", {
  .skip_if_no_factory()
  demo <- plqr_demo(n = 60, px = 2, tau = 0.5, seed = 123)
  df <- demo$data; x0 <- demo$x0
  XC <- demo$cols$XC; D_cols <- demo$cols$D; mD_cols <- demo$cols$mD; q_col <- demo$cols$q

  # tau invalid
  expect_error(
    score_plqr_component(x0, XC, D_cols, q_col, mD_cols,
                         beta = 0, j = 1L, tau = 1, s = 5, k = 1),
    regexp = "tau.*\\(0,\\s*1\\)", ignore.case = TRUE
  )
  expect_error(
    score_plqr_component(x0, XC, D_cols, q_col, mD_cols,
                         beta = 0, j = 1L, tau = -0.1, s = 5, k = 1),
    regexp = "tau.*\\(0,\\s*1\\)", ignore.case = TRUE
  )

  # β length invalid (make j valid so only β triggers)
  expect_error(
    score_plqr_component(x0, XC, D_cols, q_col, mD_cols,
                         beta = c(0, 0), j = 1L, tau = 0.5, s = 5, k = 1),
    regexp = "beta.*length\\s*\\(D_cols\\)", ignore.case = TRUE
  )

  # j out of range (make β valid so only j triggers)
  expect_error(
    score_plqr_component(x0, XC, D_cols, q_col, mD_cols,
                         beta = 0, j = 2L, tau = 0.5, s = 5, k = 1),
    regexp = "`j`.*1..length\\(D_cols\\)", ignore.case = TRUE
  )

  # mD length mismatch
  expect_error(
    score_plqr_component(x0, XC, D_cols, q_col, c(mD_cols, "bogus"),
                         beta = 0, j = 1L, tau = 0.5, s = 5, k = 1),
    regexp = "mD_cols.*same length.*D_cols", ignore.case = TRUE
  )
})


test_that("[g] score equals ψ_tau(u_i) * (D_ij - mD_ij) for singletons (R path)", {
  df <- plqr_demo(n = 80, px = 2, tau = 0.4, seed = 5)
  x0 <- colMeans(df[, 1:2, drop = FALSE])
  XC <- colnames(df)[1:2]; j <- 1L; beta <- 0; tau <- 0.4

  sc  <- score_plqr_component(x0, XC, D_cols = "D1", q_col = "q", mD_cols = "mD1",
                              beta = beta, j = j, tau = tau, s = 1, k = 1)
  ctx <- score_setup(sc, df); ctx$ptr_ctx <- NULL

  subs <- .singletons(nrow(df))
  vals <- score_eval_batch(sc, subs, df, ctx)

  # direct R formula
  thr  <- as.vector(df$D1 * beta + df$q)
  u    <- df$Y - thr
  rDj  <- df$D1 - df$mD1
  gR   <- psi_tau(u, tau) * rDj

  expect_equal(vals, gR, tolerance = 0)
  expect_equal(score_eval(sc, 17L, df, ctx), gR[17], tolerance = 0)
})

test_that("[compiled] compiled equals R for singletons (if available)", {
  df <- plqr_demo(n = 90, px = 2, tau = 0.5, seed = 7)
  x0 <- colMeans(df[, 1:2, drop = FALSE])
  XC <- colnames(df)[1:2]; j <- 1L; beta <- 0; tau <- 0.5

  sc   <- score_plqr_component(x0, XC, "D1", "q", "mD1",
                               beta = beta, j = j, tau = tau, s = 1, k = 1)
  ctxC <- score_setup(sc, df)
  ctxR <- ctxC; ctxR$ptr_ctx <- NULL

  subs   <- .singletons(nrow(df))
  vals_C <- score_eval_batch(sc, subs, df, ctxC)
  vals_R <- score_eval_batch(sc, subs, df, ctxR)

  expect_equal(vals_C, vals_R, tolerance = 0)
})

test_that("[inputs] non-finite inputs error on setup/eval", {
  df <- plqr_demo(n = 30, px = 2, tau = 0.5, seed = 9)
  x0 <- colMeans(df[, 1:2, drop = FALSE])

  # Non-finite Y
  df_bad <- df; df_bad$Y[3] <- NA_real_
  sc <- score_plqr_component(x0, colnames(df)[1:2], "D1", "q", "mD1",
                             beta = 0, j = 1L, tau = 0.5, s = 1, k = 1)
  expect_error(score_setup(sc, df_bad), "Y.*finite", ignore.case = TRUE)
})

test_that("[range] with rDj ≡ 1, per-observation g_i ∈ [tau - 1, tau]", {
  df <- plqr_demo(n = 60, px = 2, tau = 0.5, seed = 10)
  x0 <- colMeans(df[, 1:2, drop = FALSE])
  XC <- colnames(df)[1:2]; j <- 1L; beta <- 0; tau <- 0.5

  sc  <- score_plqr_component(x0, XC, "D1", "q", "mD1",
                              beta = beta, j = j, tau = tau, s = 1, k = 1)
  ctx <- score_setup(sc, df); ctx$ptr_ctx <- NULL

  subs <- .singletons(nrow(df))
  vals <- score_eval_batch(sc, subs, df, ctx)

  expect_true(all(vals >= (tau - 1) - 1e-12))
  expect_true(all(vals <=  tau      + 1e-12))
})
