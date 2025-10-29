# tests/testthat/test-score-cate.R
# -------------------------------------------------------------------
# CATE score tests — PURE per-observation score only (no weights/USTAT)
# We evaluate only singleton subsets (s=1,k=1), so the score equals the
# per-row quantity g_i = rho_i.
# -------------------------------------------------------------------

testthat::local_edition(3)

# --- helpers ---------------------------------------------------------
cate_demo <- function(n = 120, p = 2, seed = 1L) {
  set.seed(seed)
  X <- matrix(rnorm(n * p), n, p, dimnames = list(NULL, paste0("X", 1:p)))
  rho <- rnorm(n)
  data.frame(X, rho = rho)
}

.singletons <- function(n) matrix(seq_len(n), ncol = 1L)

# --- tests -----------------------------------------------------------

test_that("[factory] XOR (rho vs rho_col) and basic dimension checks", {
  df <- cate_demo(n = 30, p = 3, seed = 1)
  x0 <- colMeans(df[, 1:3, drop = FALSE])

  # Must supply exactly one of rho or rho_col
  expect_error(score_cate_rho(x0 = x0, x_cols = 1:3, s = 1, k = 1),
    "exactly one|mutually exclusive",
    ignore.case = TRUE
  )
  expect_error(
    score_cate_rho(
      x0 = x0, x_cols = 1:3, s = 1, k = 1,
      rho = df$rho, rho_col = "rho"
    ),
    "exactly one|mutually exclusive",
    ignore.case = TRUE
  )

  # x0 length mismatch (setup-time)
  expect_error(
    score_cate_rho(x0 = rep(0, 4), x_cols = 1:3, s = 1, k = 1, rho = df$rho) |> score_setup(df),
    "x0 length|ncol",
    ignore.case = TRUE
  )

  # rho length mismatch (setup-time)
  expect_error(
    score_cate_rho(x0 = x0, x_cols = 1:3, s = 1, k = 1, rho = df$rho[-1]) |> score_setup(df),
    "rho length|nrow",
    ignore.case = TRUE
  )
})

test_that("[g] score equals rho_i for singleton subsets (R path)", {
  df <- cate_demo(n = 50, p = 2, seed = 2)
  x0 <- colMeans(df[, 1:2, drop = FALSE])

  sc <- score_cate_rho(x0 = x0, x_cols = 1:2, s = 1, k = 1, rho = df$rho)
  ctx <- score_setup(sc, df)
  ctx$ptr_ctx <- NULL # force R fallback

  subs <- .singletons(nrow(df))
  vals <- score_eval_batch(sc, subs, df, ctx)

  expect_equal(vals, df$rho, tolerance = 0)
  # Also check single-index eval == batch row
  one <- score_eval(sc, 17L, df, ctx)
  expect_equal(one, df$rho[17], tolerance = 0)
})

test_that("[compiled] compiled equals R for singleton subsets (if available)", {
  df <- cate_demo(n = 60, p = 2, seed = 3)
  x0 <- colMeans(df[, 1:2, drop = FALSE])

  sc <- score_cate_rho(x0 = x0, x_cols = 1:2, s = 1, k = 1, rho = df$rho)
  ctxC <- score_setup(sc, df) # compiled (if ptr available)
  ctxR <- ctxC
  ctxR$ptr_ctx <- NULL # R path

  subs <- .singletons(nrow(df))
  vals_C <- score_eval_batch(sc, subs, df, ctxC)
  vals_R <- score_eval_batch(sc, subs, df, ctxR)

  expect_equal(vals_C, vals_R, tolerance = 0)
})

test_that("[inputs] non-finite rho triggers an error on evaluation", {
  df <- cate_demo(n = 20, p = 2, seed = 4)
  df$rho[7] <- NA_real_
  x0 <- colMeans(df[, 1:2, drop = FALSE])

  sc <- score_cate_rho(x0 = x0, x_cols = 1:2, s = 1, k = 1, rho_col = "rho")
  ctx <- score_setup(sc, df)
  ctx$ptr_ctx <- NULL

  # Try to evaluate the problematic index directly
  expect_error(score_eval(sc, 7L, df, ctx), "non-finite", ignore.case = TRUE)
})
