# tests/testthat/test-score-cate.R
# -------------------------------------------------------------------
# CATE score tests — PURE per-observation score; no weights/USTAT layer
# The kernel stores g_i = rho_i and returns sum_{i in S} g_i for any subset S.
# This suite uses the extended generator cate_demo_sets() if available.
# -------------------------------------------------------------------

testthat::local_edition(3)

# --- local helpers ---------------------------------------------------

# Fallback mini generator if cate_demo_sets() is not present
.cate_minigen <- function(n = 60L, p = 2L, seed = 1L) {
  .test_rng_scope({
    X <- matrix(rnorm(n * p), n, p, dimnames = list(NULL, paste0("X", 1:p)))
    rho <- rnorm(n)
    df <- data.frame(X, aipw_oracle = rho)
    list(
      data = df,
      x0 = rep(0, p),
      x_cols = colnames(X),
      s = 1L,
      subs_all = matrix(seq_len(n), ncol = 1L),
      subs_empty = matrix(integer(0), ncol = 1L),
      subs_2col = cbind(seq_len(n), rev(seq_len(n))),
      cols = list(rho = "aipw_oracle")
    )
  }, seed = seed)
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
    score_cate_rho(x0 = scn$x0, x_cols = scn$x_cols, s = 1L),
    "exactly one|mutually exclusive",
    ignore.case = TRUE
  )
  expect_error(
    score_cate_rho(
      x0 = scn$x0, x_cols = scn$x_cols, s = 1L,
      rho = scn$data[[rho_col]], rho_col = rho_col
    ),
    "exactly one|mutually exclusive",
    ignore.case = TRUE
  )

  # s must be positive integer (no enforcement to s == 1)
  expect_error(
    score_cate_rho(x0 = scn$x0, x_cols = scn$x_cols, s = 0L, rho_col = rho_col),
    "`s` must be a positive integer\\.",
    ignore.case = FALSE
  )

  # `k` was a compatibility-only argument and has been removed
  expect_error(
    score_cate_rho(x0 = scn$x0, x_cols = scn$x_cols, s = 1L, k = 1L, rho_col = rho_col),
    "unused argument"
  )

  # x0 length mismatch (setup-time)
  sc_badx0 <- score_cate_rho(
    x0 = rep(0, length(scn$x_cols) + 1L), x_cols = scn$x_cols, s = 1L, rho_col = rho_col
  )
  expect_error(score_setup(sc_badx0, scn$data),
    "x0 length must equal ncol\\(X\\)\\.",
    ignore.case = FALSE
  )

  # rho length mismatch (setup-time) via external vector
  sc_bad_r <- score_cate_rho(
    x0 = scn$x0, x_cols = scn$x_cols, s = 1L, rho = scn$data[[rho_col]][-1]
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
    sc <- score_cate_rho(x0 = scn$x0, x_cols = scn$x_cols, s = ss, rho_col = rho_col)
    ctx <- score_setup(sc, scn$data)
    expect_identical(sc$order(nrow(scn$data), scn$data, ctx), as.integer(ss))
  }
})

test_that("[eval] singletons return exactly rho_i (rho_col path)", {
  fx <- .get_cate_fx(n = 55, p = 2)
  scn <- fx[[1L]]
  rho_col <- scn$cols$rho

  sc <- score_cate_rho(x0 = scn$x0, x_cols = scn$x_cols, s = 1L, rho_col = rho_col)
  ctx <- score_setup(sc, scn$data)

  subs1 <- matrix(seq_len(nrow(scn$data)), ncol = 1L)
  vals <- score_eval_batch(sc, subs1, scn$data, ctx)
  expect_equal(vals, scn$data[[rho_col]], tolerance = 1e-12)

  for (i in c(1L, 17L, nrow(scn$data))) {
    expect_equal(score_eval(sc, i, scn$data, ctx), scn$data[[rho_col]][i], tolerance = 1e-12)
  }

  # empty batch
  expect_identical(score_eval_batch(sc, scn$subs_empty, scn$data, ctx), numeric(0L))
})

test_that("[eval] general subsets (multi-column) sum per-row rho", {
  fx <- .get_cate_fx(n = 60, p = 2)
  scn <- fx[[1L]]
  rho_col <- scn$cols$rho

  # Each subset has the declared order s = 2.
  sc <- score_cate_rho(x0 = scn$x0, x_cols = scn$x_cols, s = 2L, rho_col = rho_col)
  ctx <- score_setup(sc, scn$data)

  subs2 <- if (!is.null(scn$subs_2col)) scn$subs_2col else cbind(seq_len(nrow(scn$data)), rev(seq_len(nrow(scn$data))))
  vals <- score_eval_batch(sc, subs2, scn$data, ctx)

  # direct R row sums of rho over those columns
  rho <- scn$data[[rho_col]]
  gR <- rowSums(matrix(rho[subs2], nrow = nrow(subs2), ncol = ncol(subs2)))
  expect_equal(vals, gR, tolerance = 1e-12)

  # A length-s vector represents one subset, as in the scalar evaluator.
  idx <- c(1L, 17L)
  out1 <- score_eval_batch(sc, idx, scn$data, ctx)
  expect_equal(out1, sum(rho[idx]), tolerance = 1e-12)
  shape_error <- expect_error(
    score_eval_batch(sc, seq_len(nrow(scn$data)), scn$data, ctx),
    "must have 2 elements", class = "rdistnn_input_error"
  )
  expect_identical(shape_error$code, "argument_invalid")
  expect_identical(shape_error$argument, "idxs")
})

test_that("[inputs] integer x_cols indexing works", {
  fx <- .get_cate_fx(n = 48, p = 3)
  scn <- fx[[1L]]
  rho_col <- scn$cols$rho

  sc <- score_cate_rho(x0 = scn$x0, x_cols = seq_along(scn$x_cols), s = 1L, rho_col = rho_col)
  ctx <- score_setup(sc, scn$data)

  subs1 <- matrix(seq_len(nrow(scn$data)), ncol = 1L)
  vals <- score_eval_batch(sc, subs1, scn$data, ctx)
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

  sc <- score_cate_rho(x0 = scn$x0, x_cols = scn$x_cols, s = 3L, rho_col = rho_col)
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

  sc <- score_cate_rho(x0 = scn$x0, x_cols = scn$x_cols, s = 3L, rho_col = rho_col)
  ctx <- score_setup(sc, df_bad)

  expect_error(score_eval(sc, c(3L, 4L, 5L), df_bad, ctx), "g contains non-finite values", ignore.case = FALSE)
  subs3 <- cbind(1:10, 11:20, 21:30)
  expect_error(score_eval_batch(sc, subs3, df_bad, ctx), "g contains non-finite values", ignore.case = FALSE)
})


test_that("[object] basic shape/flags", {
  fx <- .get_cate_fx(n = 28, p = 2)
  scn <- fx[[1L]]
  rho <- scn$data[[scn$cols$rho]]

  sc <- score_cate_rho(x0 = scn$x0, x_cols = scn$x_cols, s = 4L, rho = rho)
  expect_true(inherits(sc, "rdistnn_ustat_score"))
  expect_type(sc$eval, "closure")
  expect_type(sc$eval_batch, "closure")
  expect_false(isTRUE(sc$randomized))
  expect_true(isTRUE(sc$symmetric))
})

test_that("[parity] compiled batch equals direct row sums", {
  .test_rng_scope({
    fx <- .get_cate_fx(n = 50, p = 2)
    scn <- fx[[1L]]
    rho_col <- scn$cols$rho

    sc <- score_cate_rho(x0 = scn$x0, x_cols = scn$x_cols, s = 3L, rho_col = rho_col)
    ctx <- score_setup(sc, scn$data)

    subs <- cbind(
      seq_len(nrow(scn$data)),
      rev(seq_len(nrow(scn$data))),
      sample(seq_len(nrow(scn$data)))
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
    out_compiled <- score_eval_batch(sc, subs, scn$data, ctx)
    expect_identical(calls, 1L)
    rho <- scn$data[[rho_col]]
    gR <- rowSums(matrix(rho[subs], nrow = nrow(subs), ncol = ncol(subs)))
    expect_equal(out_compiled, gR, tolerance = 1e-12)
  }, seed = 123L)
})
