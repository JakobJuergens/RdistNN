# tests/testthat/test-fastpaths-errors.R
# ------------------------------------------------------------------------------

.skip_if_stale_compiled <- function() {
  if (requireNamespace("pkgbuild", quietly = TRUE) &&
      !pkgbuild::has_build_tools(debug = FALSE)) {
    skip("Build tools are unavailable; compiled fastpath edge tests are skipped.")
  }

  dll <- getLoadedDLLs()[["RdistNN"]]
  if (is.null(dll)) {
    skip("RdistNN DLL is not loaded in this session.")
  }

  dll_path <- dll[["path"]]
  src_dir <- file.path(getwd(), "src")
  if (!file.exists(dll_path) || !dir.exists(src_dir)) {
    return(invisible(NULL))
  }

  src_files <- list.files(src_dir, pattern = "\\.(cpp|hpp)$", full.names = TRUE)
  if (length(src_files) == 0L) {
    return(invisible(NULL))
  }

  dll_mtime <- file.info(dll_path)$mtime
  src_mtime <- max(file.info(src_files)$mtime, na.rm = TRUE)
  if (!is.na(src_mtime) && dll_mtime < src_mtime) {
    skip("Compiled DLL appears stale relative to src/ changes; skipping fastpath failure-mode tests.")
  }
}

test_that("rustat_subset_row_sums_g handles valid inputs and rejects out-of-bounds indices", {
  .skip_if_stale_compiled()
  skip_if_not(exists("rustat_subset_row_sums_g", mode = "function"))

  g <- as.numeric(1:5)
  subs_ok <- matrix(c(1L, 2L, 5L, 4L), nrow = 2L, byrow = TRUE)
  expect_equal(rustat_subset_row_sums_g(subs_ok, g), c(3, 9))

  subs_zero <- matrix(c(0L, 2L), nrow = 1L)
  expect_error(rustat_subset_row_sums_g(subs_zero, g), "Index out of bounds")

  subs_hi <- matrix(c(6L, 1L), nrow = 1L)
  expect_error(rustat_subset_row_sums_g(subs_hi, g), "Index out of bounds")
})

test_that("rustat_cate_* validates context and indexed access bounds", {
  .skip_if_stale_compiled()
  skip_if_not(exists("rustat_cate_ctx_create", mode = "function"))
  skip_if_not(exists("rustat_cate_g_all", mode = "function"))
  skip_if_not(exists("rustat_cate_g_indexed", mode = "function"))

  expect_error(rustat_cate_ctx_create(numeric(0)), "rho must be non-empty")

  rho <- c(0.5, -1, 2)
  ctx <- rustat_cate_ctx_create(rho)
  expect_equal(as.numeric(rustat_cate_g_all(ctx)), rho)
  expect_equal(as.numeric(rustat_cate_g_indexed(c(3L, 1L), ctx)), c(2, 0.5))

  expect_error(rustat_cate_g_indexed(0L, ctx), "Index out of bounds")
  expect_error(rustat_cate_g_indexed(4L, ctx), "Index out of bounds")
})

test_that("rustat_plqr_* validates context creation and indexed access bounds", {
  .skip_if_stale_compiled()
  skip_if_not(exists("rustat_plqr_ctx_create", mode = "function"))
  skip_if_not(exists("rustat_plqr_g_all", mode = "function"))
  skip_if_not(exists("rustat_plqr_g_indexed", mode = "function"))
  skip_if_not(exists("rustat_plqr_ctx_set_beta", mode = "function"))

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
  g_idx <- as.numeric(rustat_plqr_g_indexed(c(3L, 1L), ctx))
  expect_equal(g_idx, g_all[c(3, 1)])

  expect_error(rustat_plqr_g_indexed(0L, ctx), "Index out of bounds")
  expect_error(rustat_plqr_g_indexed(4L, ctx), "Index out of bounds")

  # beta update should affect g-values
  rustat_plqr_ctx_set_beta(ctx, 1.25)
  g_new <- as.numeric(rustat_plqr_g_all(ctx))
  expect_false(isTRUE(all.equal(g_all, g_new)))
})
