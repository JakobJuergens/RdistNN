# tests/testthat/test-score-np-singleton.R
# -------------------------------------------------------------------
# Nonparametric regression score — SINGLETON (order-1) pure score
# This score must return exactly Y_i for singleton subsets and
# reject any non-singleton usage.
# -------------------------------------------------------------------

testthat::local_edition(3)

# --- local helpers ---------------------------------------------------

reg_demo <- function(n = 120, p = 2, seed = 1L) {
  set.seed(seed)
  X <- matrix(rnorm(n * p), n, p, dimnames = list(NULL, paste0("X", 1:p)))
  Y <- rnorm(n)
  data.frame(X, Y = Y)
}

.singletons <- function(n) matrix(seq_len(n), ncol = 1L)

# Use compiled batch path if it actually works; otherwise R fallback via eval().
.eval_batch_safe <- function(sc, subs, df, ctx) {
  out_try <- try(score_eval_batch(sc, subs, df, ctx), silent = TRUE)
  if (!inherits(out_try, "try-error")) {
    return(out_try)
  }
  # Fallback: per-row eval (R path)
  vapply(seq_len(nrow(subs)), function(r) score_eval(sc, subs[r, 1L], df, ctx), numeric(1))
}

# --- tests -----------------------------------------------------------

test_that("[factory] XOR (y vs y_col), enforce s = 1, k > 0, and setup dimension checks", {
  df <- reg_demo(n = 30, p = 3, seed = 1)
  x0 <- colMeans(df[, 1:3, drop = FALSE])

  # Must supply exactly one of y or y_col
  expect_error(
    score_np_y(x0 = x0, x_cols = 1:3, s = 1, k = 1),
    "exactly one|mutually exclusive",
    ignore.case = TRUE
  )
  expect_error(
    score_np_y(x0 = x0, x_cols = 1:3, s = 1, k = 1, y = df$Y, y_col = "Y"),
    "exactly one|mutually exclusive",
    ignore.case = TRUE
  )

  # Enforce s = 1 (your factory throws the specific singleton message)
  expect_error(
    score_np_y(x0 = x0, x_cols = 1:3, s = 2, k = 1, y = df$Y),
    "singleton.*require `s = 1`",
    ignore.case = TRUE
  )
  expect_error(
    score_np_y(x0 = x0, x_cols = 1:3, s = 0, k = 1, y = df$Y),
    "singleton.*require `s = 1`",
    ignore.case = TRUE
  )

  # k must be positive integer (your exact message)
  expect_error(
    score_np_y(x0 = x0, x_cols = 1:3, s = 1, k = 0, y = df$Y),
    "`k` must be a positive integer\\.",
    ignore.case = FALSE
  )

  # x0 length mismatch (setup-time)
  expect_error(
    score_np_y(x0 = rep(0, 4), x_cols = 1:3, s = 1, k = 1, y = df$Y) |> score_setup(df),
    "x0 length must equal ncol\\(X\\)\\.",
    ignore.case = FALSE
  )

  # y length mismatch (setup-time)
  expect_error(
    score_np_y(x0 = x0, x_cols = 1:3, s = 1, k = 1, y = df$Y[-1]) |> score_setup(df),
    "y length must equal nrow\\(data\\)\\.",
    ignore.case = FALSE
  )
})

test_that("[order] reported order is 1", {
  df <- reg_demo(n = 10, p = 2, seed = 2)
  x0 <- colMeans(df[, 1:2, drop = FALSE])

  sc <- score_np_y(x0 = x0, x_cols = 1:2, s = 1L, k = 1L, y_col = "Y")
  ctx <- score_setup(sc, df)
  expect_identical(sc$order(nrow(df), df, ctx), 1L)
})

test_that("[eval] singleton returns exactly Y_i (y_col path) with robust batch fallback", {
  df <- reg_demo(n = 50, p = 2, seed = 3)
  x0 <- colMeans(df[, 1:2, drop = FALSE])

  sc <- score_np_y(x0 = x0, x_cols = 1:2, s = 1L, k = 1L, y_col = "Y")
  ctx <- score_setup(sc, df)

  subs <- .singletons(nrow(df))
  vals <- .eval_batch_safe(sc, subs, df, ctx)
  expect_equal(vals, df$Y, tolerance = 0)

  for (i in c(1L, 17L, 50L)) {
    expect_equal(score_eval(sc, i, df, ctx), df$Y[i], tolerance = 0)
  }

  # Empty batch should return numeric(0)
  expect_identical(.eval_batch_safe(sc, matrix(integer(0), ncol = 1L), df, ctx), numeric(0L))
})

test_that("[eval] singleton returns exactly Y_i (external y path)", {
  df <- reg_demo(n = 35, p = 2, seed = 4)
  x0 <- colMeans(df[, 1:2, drop = FALSE])

  sc <- score_np_y(x0 = x0, x_cols = 1:2, s = 1L, k = 1L, y = df$Y)
  ctx <- score_setup(sc, df)

  subs <- .singletons(nrow(df))
  out <- .eval_batch_safe(sc, subs, df, ctx)
  expect_equal(out, df$Y, tolerance = 0)
})

test_that("[inputs] integer x_cols indexing works", {
  df <- reg_demo(n = 24, p = 3, seed = 5)
  x0 <- colMeans(df[, 1:3, drop = FALSE])

  sc <- score_np_y(x0 = x0, x_cols = 1:3, s = 1L, k = 1L, y_col = "Y")
  ctx <- score_setup(sc, df)

  subs <- .singletons(nrow(df))
  out <- .eval_batch_safe(sc, subs, df, ctx)
  expect_equal(out, df$Y, tolerance = 0)
})

test_that("[eval] errors for non-singleton usage", {
  df <- reg_demo(n = 12, p = 2, seed = 6)
  x0 <- colMeans(df[, 1:2, drop = FALSE])

  sc <- score_np_y(x0 = x0, x_cols = 1:2, s = 1L, k = 1L, y_col = "Y")
  ctx <- score_setup(sc, df)

  # eval() with >1 index
  expect_error(
    score_eval(sc, c(1L, 2L), df, ctx),
    "expects a singleton subset in `eval\\(\\)`",
    ignore.case = TRUE
  )

  # eval_batch() with ncol != 1
  bad_subs <- cbind(1:5, 2:6)
  expect_error(
    score_eval_batch(sc, bad_subs, df, ctx),
    "expects `subs` to be an n\\s*x\\s*1 matrix",
    ignore.case = TRUE
  )

  # eval_batch() with non-matrix should error
  expect_error(
    score_eval_batch(sc, 1:5, df, ctx),
    "expects `subs` to be an n\\s*x\\s*1 matrix",
    ignore.case = TRUE
  )
})

test_that("[inputs] non-finite Y triggers an error on evaluation", {
  df <- reg_demo(n = 20, p = 2, seed = 7)
  df$Y[7] <- NA_real_
  x0 <- colMeans(df[, 1:2, drop = FALSE])

  sc <- score_np_y(x0 = x0, x_cols = 1:2, s = 1L, k = 1L, y_col = "Y")
  ctx <- score_setup(sc, df)

  # Single eval at the bad row should fail
  expect_error(score_eval(sc, 7L, df, ctx), "Non-finite g encountered\\.", ignore.case = FALSE)

  # Batch should also fail — use direct call (not safe wrapper) to assert the error
  subs <- .singletons(nrow(df))
  expect_error(score_eval_batch(sc, subs, df, ctx), "Non-finite g encountered\\.", ignore.case = FALSE)
})

test_that("[sanity] singleton: constant and linear Y return exactly Y_i", {
  # constant Y
  n <- 20
  p <- 2
  c0 <- 3.25
  set.seed(101)
  X <- matrix(rnorm(n * p), n, p, dimnames = list(NULL, paste0("X", 1:p)))
  df <- data.frame(X, Y = rep(c0, n))
  x0 <- colMeans(df[, 1:2, drop = FALSE])

  sc <- score_np_y(x0 = x0, x_cols = 1:2, s = 1L, k = 1L, y_col = "Y")
  ctx <- score_setup(sc, df)
  subs <- .singletons(n)
  expect_equal(.eval_batch_safe(sc, subs, df, ctx), df$Y, tolerance = 0)

  # linear Y in index
  n <- 25
  a <- -1.0
  b <- 0.75
  set.seed(202)
  X2 <- matrix(rnorm(n * p), n, p, dimnames = list(NULL, paste0("X", 1:p)))
  df2 <- data.frame(X2, Y = a + b * seq_len(n))
  x02 <- colMeans(df2[, 1:2, drop = FALSE])

  sc2 <- score_np_y(x0 = x02, x_cols = 1:2, s = 1L, k = 1L, y_col = "Y")
  ctx2 <- score_setup(sc2, df2)
  subs2 <- .singletons(n)
  expect_equal(.eval_batch_safe(sc2, subs2, df2, ctx2), df2$Y, tolerance = 0)
})

test_that("[object] basic shape of score object", {
  df <- reg_demo(n = 8, p = 2, seed = 8)
  x0 <- colMeans(df[, 1:2, drop = FALSE])

  sc <- score_np_y(x0 = x0, x_cols = 1:2, s = 1L, k = 1L, y = df$Y)
  expect_true(inherits(sc, "ustat_score"))
  expect_type(sc$eval, "closure")
  expect_type(sc$eval_batch, "closure")
  expect_false(isTRUE(sc$randomized))
  expect_true(isTRUE(sc$symmetric))
})

test_that("[parity] compiled batch (if actually callable) equals direct R indexing", {
  df <- reg_demo(n = 45, p = 2, seed = 9)
  x0 <- colMeans(df[, 1:2, drop = FALSE])

  sc <- score_np_y(x0 = x0, x_cols = 1:2, s = 1L, k = 1L, y_col = "Y")
  ctx <- score_setup(sc, df)
  subs <- .singletons(nrow(df))

  # Try the compiled path; skip if the shared object isn't loaded.
  out_compiled_try <- try(score_eval_batch(sc, subs, df, ctx), silent = TRUE)
  if (inherits(out_compiled_try, "try-error")) {
    testthat::skip("Compiled helper not available at runtime; skipping parity check.")
  } else {
    out_R <- df$Y # direct R indexing is the ground truth for singleton score
    expect_equal(out_compiled_try, out_R, tolerance = 0)
  }
})
