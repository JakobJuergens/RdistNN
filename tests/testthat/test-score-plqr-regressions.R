# Pin large decimal fixtures to binary64 bits: R builds can parse decimals
# with different intermediate precision. Exact score assertions stay exact.
test_that("PLQR subset sums retain signed cancellation in scalar and batch engines", {
  old <- options(RdistNN.use_compiled = TRUE)
  on.exit(options(old), add = TRUE)
  dat <- data.frame(X = 1:3, D = c(2^101, 2, -2^101), mD = 0, g = 0, Y = 1)
  permutations <- rbind(c(1L, 2L, 3L), c(1L, 3L, 2L), c(2L, 1L, 3L),
    c(2L, 3L, 1L), c(3L, 1L, 2L), c(3L, 2L, 1L))
  score <- score_plqr_component(0, "X", "D", "g", "mD", "Y",
    theta = 0, tau = .5, s = 3L)
  for (native in c(TRUE, FALSE)) {
    options(RdistNN.use_compiled = native)
    if (native) {
      expect_false(is.null(.rdistnn_compiled_fun("rustat_plqr_ctx_create")))
      expect_false(is.null(.rdistnn_compiled_fun("rustat_plqr_g_all")))
      expect_false(is.null(.rdistnn_compiled_fun("rustat_subset_row_sums_g")))
    }
    for (row_order in seq_len(nrow(permutations))) {
      local_data <- dat[permutations[row_order, ], , drop = FALSE]
      ctx <- score_setup(score, local_data)
      expect_false("ptr_ctx" %in% names(ctx))
      restored <- unserialize(serialize(ctx, NULL))
      expect_identical(restored, ctx)
      expect_identical(score_eval(score, 1:3, local_data, restored), 1)
      # The represented contributions are exactly 2^100, 1 and -2^100.
      expect_identical(score_eval(score, 1:3, local_data, ctx), 1)
      expect_identical(score_eval_batch(score, permutations, local_data, ctx), rep(1, 6))
    }
  }
})

test_that("PLQR score factories reject missing and nonfinite theta before setup", {
  old <- options(RdistNN.use_compiled = TRUE)
  on.exit(options(old), add = TRUE)
  for (native in c(TRUE, FALSE)) {
    options(RdistNN.use_compiled = native)
    for (theta in c(NA_real_, NaN, Inf, -Inf)) {
      expect_error(score_plqr_component(0, "X", "D", "g", "mD", "Y",
        theta = theta, s = 1L), "`theta`.*(NA|non-finite)")
    }
  }
})

test_that("PLQR score component validation identifies invalid j without coercion warnings", {
  for (j in c(NA_real_, NaN, Inf, -Inf, 2^40)) {
    expect_no_warning(expect_error(score_plqr_component(0, "X", "D", "g", "mD", "Y",
      theta = 0, j = j, s = 1L), "`j` must be an integer in 1..length"))
  }
})

test_that("exact PLQR threshold cancellation retains inclusive equality", {
  old <- options(RdistNN.use_compiled = TRUE)
  on.exit(options(old), add = TRUE)
  dat <- data.frame(X = 0, D1 = 0x1.1ccf385ebc8a0p+1023, D2 = 0x1.1ccf385ebc8a0p+1023,
    mD1 = 0, mD2 = 0, g = 0, Y = 0)
  score <- score_plqr_component(0, "X", c("D1", "D2"), "g",
    c("mD1", "mD2"), "Y", theta = c(2, -2), s = 1L)
  for (native in c(TRUE, FALSE)) {
    options(RdistNN.use_compiled = native)
    if (native) expect_false(is.null(.rdistnn_compiled_fun("rustat_plqr_g_all")))
    expect_identical(score_setup(score, dat)$score, -0x1.1ccf385ebc8a0p+1022)
  }
  options(RdistNN.use_compiled = TRUE)
  create <- .rdistnn_compiled_fun("rustat_plqr_ctx_create")
  ctx <- create(dat$Y, as.matrix(dat[c("D1", "D2")]), dat$g,
    as.matrix(dat[c("mD1", "mD2")]), c(2, -2), .5, 1L)
  expect_identical(.rdistnn_compiled_fun("rustat_plqr_g_all")(ctx), -0x1.1ccf385ebc8a0p+1022)
})
test_that("PLQR score factories capture selectors and targets before setup", {
  dat <- data.frame(X = 1:3, D = 1:3, g = 0, mD = 0, Y = 1:3)
  x0 <- 0
  x <- "X"
  d <- "D"
  g <- "g"
  md <- "mD"
  y <- "Y"
  score <- score_plqr_component(x0, x, d, g, md, y, theta = 0, s = 1L)
  x0 <- c(0, 0)
  x <- d <- g <- md <- y <- "missing"
  ctx <- score_setup(score, dat)
  expect_identical(ctx$score, c(.5, 1, 1.5))
})

test_that("PLQR validation aggregation preserves unexpected condition identities", {
  failure <- structure(list(message = "resource sentinel", call = NULL),
    class = c("resourceError", "error", "condition"))
  local_mocked_bindings(.validate_scalar_prob = function(...) stop(failure), .package = "RdistNN")
  caught <- tryCatch(score_plqr_component(0, "X", "D", "g", "mD", "Y", theta = 0, s = 1L),
    error = identity)
  expect_identical(caught, failure)
})
