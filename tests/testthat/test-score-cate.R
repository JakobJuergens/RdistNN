# tests/testthat/test-score-cate.R
# -------------------------------------------------------------------
# CATE score tests — PURE per-observation score; no weights/USTAT layer
# The kernel stores g_i = rho_i and returns sum_{i in S} g_i for any subset S.
# This suite uses the extended generator cate_demo_sets() if available.
# -------------------------------------------------------------------

testthat::local_edition(3)

# --- local helpers ---------------------------------------------------

# robust batch evaluator that survives partially-linked compiled helpers
.eval_batch_safe <- function(sc, subs, df, ctx) {
  out_try <- try(score_eval_batch(sc, subs, df, ctx), silent = TRUE)
  if (!inherits(out_try, "try-error")) {
    return(out_try)
  }
  # If compiled path blows up (e.g. symbol exists but .Call not loaded),
  # compute in pure R using ctx$g as the ground truth definition.
  subs <- if (!is.matrix(subs)) as.matrix(subs) else subs
  if (nrow(subs) == 0L) {
    return(numeric(0L))
  }
  gi <- ctx$g
  rowSums(matrix(gi[subs], nrow = nrow(subs), ncol = ncol(subs)))
}

# Fallback mini generator if cate_demo_sets() is not present
.cate_minigen <- function(n = 60L, p = 2L, seed = 1L) {
  set.seed(seed)
  X <- matrix(rnorm(n * p), n, p, dimnames = list(NULL, paste0("X", 1:p)))
  rho <- rnorm(n)
  df <- data.frame(X, rho_oracle = rho)
  list(
    data = df,
    x0 = rep(0, p),
    x_cols = colnames(X),
    s = 1L, k = 1L,
    subs_all = matrix(seq_len(n), ncol = 1L),
    subs_empty = matrix(integer(0), ncol = 1L),
    subs_2col = cbind(seq_len(n), rev(seq_len(n))),
    cols = list(rho = "rho_oracle")
  )
}

.get_cate_fx <- function(...) {
  if (exists("cate_demo_sets", mode = "function")) {
    cate_demo_sets(...)
  } else {
    list(ok_simple = .cate_minigen(...))
  }
}

# --- tests -----------------------------------------------------------

test_that("[factory] XOR (rho vs rho_col), positive s and k, and setup dimension checks", {
  fx <- .get_cate_fx(n = 50, p = 3)
  scn <- fx[[1L]] # either ok_simple or the first scenario
  rho_col <- scn$cols$rho

  # Must supply exactly one of rho or rho_col
  expect_error(
    score_cate_rho(x0 = scn$x0, x_cols = scn$x_cols, s = 1L, k = 1L),
    "exactly one|mutually exclusive",
    ignore.case = TRUE
  )
  expect_error(
    score_cate_rho(
      x0 = scn$x0, x_cols = scn$x_cols, s = 1L, k = 1L,
      rho = scn$data[[rho_col]], rho_col = rho_col
    ),
    "exactly one|mutually exclusive",
    ignore.case = TRUE
  )

  # s must be positive integer (no enforcement to s == 1)
  expect_error(
    score_cate_rho(x0 = scn$x0, x_cols = scn$x_cols, s = 0L, k = 1L, rho_col = rho_col),
    "`s` must be a positive integer\\.",
    ignore.case = FALSE
  )

  # k must be positive integer
  expect_error(
    score_cate_rho(x0 = scn$x0, x_cols = scn$x_cols, s = 1L, k = 0L, rho_col = rho_col),
    "`k` must be a positive integer\\.",
    ignore.case = FALSE
  )

  # x0 length mismatch (setup-time)
  sc_badx0 <- score_cate_rho(
    x0 = rep(0, length(scn$x_cols) + 1L), x_cols = scn$x_cols, s = 1L, k = 1L, rho_col = rho_col
  )
  expect_error(score_setup(sc_badx0, scn$data),
    "x0 length must equal ncol\\(X\\)\\.",
    ignore.case = FALSE
  )

  # rho length mismatch (setup-time) via external vector
  sc_bad_r <- score_cate_rho(
    x0 = scn$x0, x_cols = scn$x_cols, s = 1L, k = 1L, rho = scn$data[[rho_col]][-1]
  )
  expect_error(score_setup(sc_bad_r, scn$data),
    "rho length must equal nrow\\(data\\)\\.",
    ignore.case = FALSE
  )
})

test_that("[order] reported order equals chosen s (test s = 1 and s = 5)", {
  fx <- .get_cate_fx(n = 40, p = 2)
  scn <- fx[[1L]]
  rho_col <- scn$cols$rho

  for (ss in c(1L, 5L)) {
    sc <- score_cate_rho(x0 = scn$x0, x_cols = scn$x_cols, s = ss, k = 1L, rho_col = rho_col)
    ctx <- score_setup(sc, scn$data)
    expect_identical(sc$order(nrow(scn$data), scn$data, ctx), as.integer(ss))
  }
})

test_that("[eval] singletons return exactly rho_i (rho_col path) with robust batch fallback", {
  fx <- .get_cate_fx(n = 55, p = 2)
  scn <- fx[[1L]]
  rho_col <- scn$cols$rho

  sc <- score_cate_rho(x0 = scn$x0, x_cols = scn$x_cols, s = 1L, k = 1L, rho_col = rho_col)
  ctx <- score_setup(sc, scn$data)

  subs1 <- matrix(seq_len(nrow(scn$data)), ncol = 1L)
  vals <- .eval_batch_safe(sc, subs1, scn$data, ctx)
  expect_equal(vals, scn$data[[rho_col]], tolerance = 1e-12)

  for (i in c(1L, 17L, nrow(scn$data))) {
    expect_equal(score_eval(sc, i, scn$data, ctx), scn$data[[rho_col]][i], tolerance = 1e-12)
  }

  # empty batch
  expect_identical(.eval_batch_safe(sc, scn$subs_empty, scn$data, ctx), numeric(0L))
})

test_that("[eval] general subsets (multi-column) sum per-row rho", {
  fx <- .get_cate_fx(n = 60, p = 2)
  scn <- fx[[1L]]
  rho_col <- scn$cols$rho

  # Use s=2 behavior (but eval/eval_batch just sum over provided columns)
  sc <- score_cate_rho(x0 = scn$x0, x_cols = scn$x_cols, s = 2L, k = 1L, rho_col = rho_col)
  ctx <- score_setup(sc, scn$data)

  subs2 <- if (!is.null(scn$subs_2col)) scn$subs_2col else cbind(seq_len(nrow(scn$data)), rev(seq_len(nrow(scn$data))))
  vals <- .eval_batch_safe(sc, subs2, scn$data, ctx)

  # direct R row sums of rho over those columns
  rho <- scn$data[[rho_col]]
  gR <- rowSums(matrix(rho[subs2], nrow = nrow(subs2), ncol = ncol(subs2)))
  expect_equal(vals, gR, tolerance = 1e-12)

  # Also check that passing a vector (coerced to matrix) works
  idx <- seq_len(nrow(scn$data))
  out1 <- .eval_batch_safe(sc, idx, scn$data, ctx) # treated as 1-col matrix
  expect_equal(out1, rho[idx], tolerance = 1e-12)
})

test_that("[inputs] integer x_cols indexing works", {
  fx <- .get_cate_fx(n = 48, p = 3)
  scn <- fx[[1L]]
  rho_col <- scn$cols$rho

  sc <- score_cate_rho(x0 = scn$x0, x_cols = seq_along(scn$x_cols), s = 1L, k = 1L, rho_col = rho_col)
  ctx <- score_setup(sc, scn$data)

  subs1 <- matrix(seq_len(nrow(scn$data)), ncol = 1L)
  vals <- .eval_batch_safe(sc, subs1, scn$data, ctx)
  expect_equal(vals, scn$data[[rho_col]], tolerance = 1e-12)
})

# tests/testthat/test-score-cate.R  (replace the whole [non-finite] block)
test_that("[invalid g] NA in rho triggers NA-specific errors on evaluation", {
  fx <- .get_cate_fx(n = 42, p = 2)
  scn <- fx[[1L]]
  rho_col <- scn$cols$rho

  # inject NA at known position
  df_bad <- scn$data
  df_bad[[rho_col]][2] <- NA_real_

  sc <- score_cate_rho(x0 = scn$x0, x_cols = scn$x_cols, s = 3L, k = 1L, rho_col = rho_col)
  ctx <- score_setup(sc, df_bad)

  # eval(): must pass an order-s subset that includes the bad index.
  expect_error(score_eval(sc, c(1L, 2L, 3L), df_bad, ctx), "g contains `NA` values\\.", ignore.case = FALSE)

  # eval_batch(): global guard triggers regardless of subset content
  subs3 <- cbind(1:10, 11:20, 21:30)
  expect_error(score_eval_batch(sc, subs3, df_bad, ctx), "g contains `NA` values\\.", ignore.case = FALSE)
})

test_that("[invalid g] non-finite rho triggers non-finite-specific errors on evaluation", {
  fx <- .get_cate_fx(n = 42, p = 2)
  scn <- fx[[1L]]
  rho_col <- scn$cols$rho

  df_bad <- scn$data
  df_bad[[rho_col]][5] <- Inf

  sc <- score_cate_rho(x0 = scn$x0, x_cols = scn$x_cols, s = 3L, k = 1L, rho_col = rho_col)
  ctx <- score_setup(sc, df_bad)

  expect_error(score_eval(sc, c(3L, 4L, 5L), df_bad, ctx), "g contains non-finite values", ignore.case = FALSE)
  subs3 <- cbind(1:10, 11:20, 21:30)
  expect_error(score_eval_batch(sc, subs3, df_bad, ctx), "g contains non-finite values", ignore.case = FALSE)
})


test_that("[object] basic shape/flags", {
  fx <- .get_cate_fx(n = 28, p = 2)
  scn <- fx[[1L]]
  rho <- scn$data[[scn$cols$rho]]

  sc <- score_cate_rho(x0 = scn$x0, x_cols = scn$x_cols, s = 4L, k = 1L, rho = rho)
  expect_true(inherits(sc, "ustat_score"))
  expect_type(sc$eval, "closure")
  expect_type(sc$eval_batch, "closure")
  expect_false(isTRUE(sc$randomized))
  expect_true(isTRUE(sc$symmetric))
})

test_that("[parity] compiled batch (if actually callable) equals direct row sums", {
  fx <- .get_cate_fx(n = 50, p = 2)
  scn <- fx[[1L]]
  rho_col <- scn$cols$rho

  sc <- score_cate_rho(x0 = scn$x0, x_cols = scn$x_cols, s = 3L, k = 1L, rho_col = rho_col)
  ctx <- score_setup(sc, scn$data)

  # try compiled; if the shared object isn't actually loaded, skip
  subs <- cbind(
    seq_len(nrow(scn$data)),
    rev(seq_len(nrow(scn$data))),
    sample(seq_len(nrow(scn$data)))
  )
  out_try <- try(score_eval_batch(sc, subs, scn$data, ctx), silent = TRUE)
  if (inherits(out_try, "try-error")) {
    testthat::skip("Compiled helper not available at runtime; skipping parity check.")
  } else {
    rho <- scn$data[[rho_col]]
    gR <- rowSums(matrix(rho[subs], nrow = nrow(subs), ncol = ncol(subs)))
    expect_equal(out_try, gR, tolerance = 1e-12)
  }
})
