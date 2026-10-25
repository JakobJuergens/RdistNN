# tests/testthat/test-score-np-singleton.R
# -------------------------------------------------------------------
# Nonparametric regression score — SINGLETON (order-1) pure score
# This score must return exactly Y_i for singleton subsets and
# reject any non-singleton usage.
# -------------------------------------------------------------------

testthat::local_edition(3)

# --- local helpers ---------------------------------------------------

reg_demo <- function(n = 120, p = 2, seed = 1L) {
  .test_rng_scope({
    X <- matrix(rnorm(n * p), n, p, dimnames = list(NULL, paste0("X", 1:p)))
    Y <- rnorm(n)
    data.frame(X, Y = Y)
  }, seed = seed)
}

.singletons <- function(n) matrix(seq_len(n), ncol = 1L)

# --- tests -----------------------------------------------------------

test_that("[factory] XOR (y vs y_col), enforce s = 1, k > 0, and setup dimension checks", {
  df <- reg_demo(n = 30, p = 3, seed = 1)
  x0 <- colMeans(df[, 1:3, drop = FALSE])

  # Must supply exactly one of y or y_col
  expect_error(
    score_np_y(x0 = x0, x_cols = 1:3, s = 1),
    "exactly one|mutually exclusive",
    ignore.case = TRUE
  )
  expect_error(
    score_np_y(x0 = x0, x_cols = 1:3, s = 1, y = df$Y, y_col = "Y"),
    "exactly one|mutually exclusive",
    ignore.case = TRUE
  )

  # Enforce s = 1 (your factory throws the specific singleton message)
  expect_error(
    score_np_y(x0 = x0, x_cols = 1:3, s = 2, y = df$Y),
    "singleton.*require `s = 1`",
    ignore.case = TRUE
  )
  expect_error(
    score_np_y(x0 = x0, x_cols = 1:3, s = 0, y = df$Y),
    "singleton.*require `s = 1`",
    ignore.case = TRUE
  )



  # x0 length mismatch (setup-time)
  expect_error(
    score_np_y(x0 = rep(0, 4), x_cols = 1:3, s = 1, y = df$Y) |> score_setup(df),
    "x0 length must equal ncol\\(X\\)\\.",
    ignore.case = FALSE
  )

  # y length mismatch (setup-time)
  expect_error(
    score_np_y(x0 = x0, x_cols = 1:3, s = 1, y = df$Y[-1]) |> score_setup(df),
    "y length must equal nrow\\(data\\)\\.",
    ignore.case = FALSE
  )
})

test_that("[order] reported order is 1", {
  df <- reg_demo(n = 10, p = 2, seed = 2)
  x0 <- colMeans(df[, 1:2, drop = FALSE])

  sc <- score_np_y(x0 = x0, x_cols = 1:2, s = 1L, y_col = "Y")
  ctx <- score_setup(sc, df)
  expect_identical(sc$order(nrow(df), df, ctx), 1L)
})

test_that("[eval] singleton returns exactly Y_i (y_col path)", {
  df <- reg_demo(n = 50, p = 2, seed = 3)
  x0 <- colMeans(df[, 1:2, drop = FALSE])

  sc <- score_np_y(x0 = x0, x_cols = 1:2, s = 1L, y_col = "Y")
  ctx <- score_setup(sc, df)

  subs <- .singletons(nrow(df))
  vals <- score_eval_batch(sc, subs, df, ctx)
  expect_equal(vals, df$Y, tolerance = 1e-12)

  for (i in c(1L, 17L, 50L)) {
    expect_equal(score_eval(sc, i, df, ctx), df$Y[i], tolerance = 1e-12)
  }

  # Empty batch should return numeric(0)
  expect_identical(score_eval_batch(sc, matrix(integer(0), ncol = 1L), df, ctx), numeric(0L))
})

test_that("[eval] singleton returns exactly Y_i (external y path)", {
  df <- reg_demo(n = 35, p = 2, seed = 4)
  x0 <- colMeans(df[, 1:2, drop = FALSE])

  sc <- score_np_y(x0 = x0, x_cols = 1:2, s = 1L, y = df$Y)
  ctx <- score_setup(sc, df)

  subs <- .singletons(nrow(df))
  out <- score_eval_batch(sc, subs, df, ctx)
  expect_equal(out, df$Y, tolerance = 1e-12)
})

test_that("[inputs] integer x_cols indexing works", {
  df <- reg_demo(n = 24, p = 3, seed = 5)
  x0 <- colMeans(df[, 1:3, drop = FALSE])

  sc <- score_np_y(x0 = x0, x_cols = 1:3, s = 1L, y_col = "Y")
  ctx <- score_setup(sc, df)

  subs <- .singletons(nrow(df))
  out <- score_eval_batch(sc, subs, df, ctx)
  expect_equal(out, df$Y, tolerance = 1e-12)
})

test_that("[eval] errors for non-singleton usage", {
  df <- reg_demo(n = 12, p = 2, seed = 6)
  x0 <- colMeans(df[, 1:2, drop = FALSE])

  sc <- score_np_y(x0 = x0, x_cols = 1:2, s = 1L, y_col = "Y")
  ctx <- score_setup(sc, df)

  # eval() with >1 index
  expect_error(
    score_eval(sc, c(1L, 2L), df, ctx),
    "idxs.*1 elements",
    ignore.case = TRUE
  )

  # eval_batch() with ncol != 1
  bad_subs <- cbind(1:5, 2:6)
  expect_error(
    score_eval_batch(sc, bad_subs, df, ctx),
    "subs.*1 columns",
    ignore.case = TRUE
  )

  # eval_batch() treats a vector as one subset, not as a batch of singletons.
  expect_error(
    score_eval_batch(sc, 1:5, df, ctx),
    "idxs.*1 elements",
    ignore.case = TRUE
  )
})

test_that("[inputs] NA in Y triggers NA-specific error on evaluation", {
  df <- reg_demo(n = 20, p = 2, seed = 7)
  df$Y[7] <- NA_real_
  x0 <- colMeans(df[, 1:2, drop = FALSE])

  sc <- score_np_y(x0 = x0, x_cols = 1:2, s = 1L, y_col = "Y")
  ctx <- score_setup(sc, df)

  # Single eval at the bad row should fail
  expect_error(score_eval(sc, 7L, df, ctx), "g contains `NA` values\\.", ignore.case = FALSE)

  # Batch should also fail — use direct call (not safe wrapper) to assert the error
  subs <- .singletons(nrow(df))
  expect_error(score_eval_batch(sc, subs, df, ctx), "g contains `NA` values\\.", ignore.case = FALSE)
})

test_that("[inputs] non-finite Y triggers non-finite-specific error on evaluation", {
  df <- reg_demo(n = 20, p = 2, seed = 7)
  df$Y[7] <- Inf
  x0 <- colMeans(df[, 1:2, drop = FALSE])

  sc <- score_np_y(x0 = x0, x_cols = 1:2, s = 1L, y_col = "Y")
  ctx <- score_setup(sc, df)

  expect_error(score_eval(sc, 7L, df, ctx), "g contains non-finite values", ignore.case = FALSE)
  subs <- .singletons(nrow(df))
  expect_error(score_eval_batch(sc, subs, df, ctx), "g contains non-finite values", ignore.case = FALSE)
})

test_that("[sanity] singleton: constant and linear Y return exactly Y_i", {
  .test_rng_scope({
    # constant Y
    n <- 20
    p <- 2
    c0 <- 3.25
    set.seed(101)
    X <- matrix(rnorm(n * p), n, p, dimnames = list(NULL, paste0("X", 1:p)))
    df <- data.frame(X, Y = rep(c0, n))
    x0 <- colMeans(df[, 1:2, drop = FALSE])

    sc <- score_np_y(x0 = x0, x_cols = 1:2, s = 1L, y_col = "Y")
    ctx <- score_setup(sc, df)
    subs <- .singletons(n)
    expect_equal(score_eval_batch(sc, subs, df, ctx), df$Y, tolerance = 1e-12)

    # linear Y in index
    n <- 25
    a <- -1.0
    b <- 0.75
    set.seed(202)
    X2 <- matrix(rnorm(n * p), n, p, dimnames = list(NULL, paste0("X", 1:p)))
    df2 <- data.frame(X2, Y = a + b * seq_len(n))
    x02 <- colMeans(df2[, 1:2, drop = FALSE])

    sc2 <- score_np_y(x0 = x02, x_cols = 1:2, s = 1L, y_col = "Y")
    ctx2 <- score_setup(sc2, df2)
    subs2 <- .singletons(n)
    expect_equal(score_eval_batch(sc2, subs2, df2, ctx2), df2$Y, tolerance = 1e-12)
  })
})

test_that("[object] basic shape of score object", {
  df <- reg_demo(n = 8, p = 2, seed = 8)
  x0 <- colMeans(df[, 1:2, drop = FALSE])

  sc <- score_np_y(x0 = x0, x_cols = 1:2, s = 1L, y = df$Y)
  expect_true(inherits(sc, "rdistnn_ustat_score"))
  expect_type(sc$eval, "closure")
  expect_type(sc$eval_batch, "closure")
  expect_false(isTRUE(sc$randomized))
  expect_true(isTRUE(sc$symmetric))
})

test_that("[parity] compiled batch equals direct R indexing", {
  df <- reg_demo(n = 45, p = 2, seed = 9)
  x0 <- colMeans(df[, 1:2, drop = FALSE])

  sc <- score_np_y(x0 = x0, x_cols = 1:2, s = 1L, y_col = "Y")
  ctx <- score_setup(sc, df)
  subs <- .singletons(nrow(df))

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
  out_R <- df$Y # direct R indexing is the ground truth for singleton score
  expect_equal(out_compiled, out_R, tolerance = 1e-12)
})
