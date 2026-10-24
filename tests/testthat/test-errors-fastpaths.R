# Required native boundaries execute the loaded package DLL. Installing or
# running a binary package does not require local build tools.

test_that("rustat_subset_row_sums_g handles valid inputs and rejects out-of-bounds indices", {
  g <- as.numeric(1:5)
  subs_ok <- matrix(c(1L, 2L, 5L, 4L), nrow = 2L, byrow = TRUE)
  expect_equal(rustat_subset_row_sums_g(subs_ok, g), c(3, 9))

  subs_zero <- matrix(c(0L, 2L), nrow = 1L)
  expect_error(rustat_subset_row_sums_g(subs_zero, g), "Index out of bounds")

  subs_hi <- matrix(c(6L, 1L), nrow = 1L)
  expect_error(rustat_subset_row_sums_g(subs_hi, g), "Index out of bounds")
})

test_that("rustat_plqr_* validates context creation and numeric scores", {
  Y <- c(0, 1, 2)
  D <- matrix(c(1, 0, 2), ncol = 1)
  mD <- matrix(c(0.2, 0.1, 0.4), ncol = 1)
  qx <- c(0, 0, 0)
  beta <- 0.5
  tau <- 0.5

  # shape/range guards at context construction
  expect_error(
    rustat_plqr_ctx_create(Y, D[1:2, , drop = FALSE], qx, mD, beta, tau, 1L),
    "must share the same number of rows"
  )
  expect_error(
    rustat_plqr_ctx_create(Y, D, qx, cbind(mD, mD), beta, tau, 1L),
    "same number of columns"
  )
  expect_error(
    rustat_plqr_ctx_create(Y, D, qx, mD, c(beta, beta), tau, 1L),
    "beta length must equal ncol\\(D\\)"
  )
  expect_error(
    rustat_plqr_ctx_create(Y, D, qx, mD, beta, tau, 2L),
    "j is out of range"
  )
  expect_error(
    rustat_plqr_ctx_create(Y, D, qx, mD, beta, 1.0, 1L),
    "tau must be in \\(0,1\\)"
  )

  ctx <- rustat_plqr_ctx_create(Y, D, qx, mD, beta, tau, 1L)
  g_all <- as.numeric(rustat_plqr_g_all(ctx))
  expect_identical(g_all, c(-.4, -.05, .8))
  changed <- rustat_plqr_ctx_create(Y, D, qx, mD, 1.25, tau, 1L)
  expect_identical(as.numeric(rustat_plqr_g_all(changed)), c(-.4, -.05, -.8))
})

# Invalid native contexts run separately: an abort must fail this expectation,
# without terminating the test runner or masking the required R error.
test_that("PLQR native contexts reject null, restored and foreign pointers safely", {
  for (kind in c("null", "restored", "foreign")) for (n in c(3L, 1100L)) {
    script <- withr::local_tempfile(fileext = ".R")
    output <- withr::local_tempfile(fileext = ".log")
    writeLines(c(
      paste0(".libPaths(", paste(deparse(.libPaths()), collapse = ""), ")"),
      'library(RdistNN)',
      'options(RdistNN.num_threads = 2L, RdistNN.use_compiled = TRUE)',
      sprintf('n <- %dL', n),
      'ctx <- RdistNN:::rustat_plqr_ctx_create(rep(1, n), matrix(1, n), rep(0, n), matrix(0, n), 0, .5, 1L)',
      'stopifnot(identical(RdistNN:::rustat_plqr_g_all(ctx), rep(.5, n)))',
      switch(kind,
        null = 'ctx <- new("externalptr")',
        restored = 'ctx <- unserialize(serialize(ctx, NULL))',
        foreign = 'ctx <- getNativeSymbolInfo("_RdistNN_rustat_dist_l2", PACKAGE = "RdistNN")$address'),
      'message <- tryCatch({RdistNN:::rustat_plqr_g_all(ctx); "no error"}, error = conditionMessage)',
      'stopifnot(grepl("Invalid PLQR context", message, fixed = TRUE))',
      'cat("EXPECTED_CONTEXT_ERROR_AND_NORMAL_COMPLETION\\n")'
    ), script)
    status <- suppressWarnings(system2(file.path(R.home("bin"), "Rscript"),
      c("--vanilla", shQuote(script)), stdout = output, stderr = output, timeout = 30))
    detail <- paste(readLines(output, warn = FALSE), collapse = "\n")
    expect_identical(status, 0L, info = paste(kind, n, detail))
    expect_match(detail, "EXPECTED_CONTEXT_ERROR_AND_NORMAL_COMPLETION", fixed = TRUE)
  }
})
