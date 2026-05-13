.iv_repo_root <- function() {
  candidates <- c(
    getwd(),
    normalizePath(file.path(getwd(), ".."), mustWork = FALSE),
    normalizePath(file.path(getwd(), "..", ".."), mustWork = FALSE)
  )
  for (cand in candidates) {
    if (file.exists(file.path(cand, "DESCRIPTION")) &&
        dir.exists(file.path(cand, "inst", "paper-replication", "iv"))) {
      return(cand)
    }
  }
  stop("Could not locate package root for IV paper-replication tests.", call. = FALSE)
}

.iv_repl_path <- function(...) {
  file.path(.iv_repo_root(), "inst", "paper-replication", "iv", ...)
}

.set_env_for_test <- function(values) {
  old <- Sys.getenv(names(values), unset = NA_character_)
  do.call(Sys.setenv, as.list(values))
  function() {
    for (nm in names(old)) {
      if (is.na(old[[nm]])) {
        Sys.unsetenv(nm)
      } else {
        do.call(Sys.setenv, stats::setNames(list(old[[nm]]), nm))
      }
    }
  }
}

test_that("IV paper replication grid has the intended main cells", {
  env <- new.env(parent = globalenv())
  source(.iv_repl_path("grid.R"), local = env)

  grid <- env$iv_paper_grid(smoke = FALSE)
  cells <- env$iv_paper_cells(grid)

  expect_equal(nrow(cells), 36L)
  expect_equal(env$iv_expected_cells(grid), 36L)
  expect_false("dgp4_iv_weak_local" %in% grid$designs)
  expect_equal(
    grid$methods,
    c("dnn_iv_fullsample", "dnn_iv_crossfit", "local_2sls", "2sls_interact")
  )

  justid <- cells[cells$design == "dgp1_iv_justid_linear", , drop = FALSE]
  overid <- cells[cells$design != "dgp1_iv_justid_linear", , drop = FALSE]
  expect_equal(unique(justid$q), 1L)
  expect_equal(nrow(justid), 4L)
  expect_equal(sort(unique(overid$q)), c(2L, 3L))
  expect_equal(nrow(overid), 32L)

  diagnostic <- env$iv_diagnostic_grid(smoke = FALSE)
  expect_equal(diagnostic$designs, "dgp4_iv_weak_local")
})

test_that("IV paper replication smoke grid is small and script-safe", {
  env <- new.env(parent = globalenv())
  source(.iv_repl_path("grid.R"), local = env)

  grid <- env$iv_paper_grid(smoke = TRUE)
  cells <- env$iv_paper_cells(grid)

  expect_equal(nrow(cells), 2L)
  expect_equal(grid$n_reps, 2L)
  expect_equal(grid$methods, "iv_oracle")
  expect_equal(unique(cells$q[cells$design == "dgp1_iv_justid_linear"]), 1L)
})

test_that("IV single-machine replication script writes smoke outputs", {
  out_dir <- tempfile("rdistnn-iv-run-")
  dir.create(out_dir, recursive = TRUE)
  restore <- .set_env_for_test(c(
    RDISTNN_IV_SMOKE = "1",
    SIM_OUTDIR = out_dir
  ))
  on.exit(restore(), add = TRUE)

  env <- new.env(parent = globalenv())
  invisible(utils::capture.output(
    suppressMessages(source(.iv_repl_path("run.R"), local = env)),
    type = "output"
  ))

  raw_file <- file.path(out_dir, "sim_iv_raw.rds")
  tab_file <- file.path(out_dir, "sim_iv_table.csv")
  expect_true(file.exists(raw_file))
  expect_true(file.exists(tab_file))

  raw <- readRDS(raw_file)
  tab <- utils::read.csv(tab_file)
  expect_equal(nrow(raw), 4L)
  expect_equal(nrow(tab), 2L)
  expect_true(all(c("method", "design", "n", "kx", "q", "n_reps") %in% names(tab)))
})

test_that("IV HPC task and merge scripts smoke-run one cell", {
  out_dir <- tempfile("rdistnn-iv-hpc-")
  dir.create(out_dir, recursive = TRUE)
  restore <- .set_env_for_test(c(
    RDISTNN_IV_SMOKE = "1",
    SIM_OUTDIR = out_dir,
    SLURM_ARRAY_TASK_ID = "1",
    SLURM_CPUS_PER_TASK = "1"
  ))
  on.exit(restore(), add = TRUE)

  task_env <- new.env(parent = globalenv())
  suppressMessages(source(.iv_repl_path("hpc", "task.R"), local = task_env))

  cell_file <- file.path(out_dir, "cells", "cell_01.rds")
  expect_true(file.exists(cell_file))
  expect_equal(nrow(readRDS(cell_file)), 2L)

  merge_env <- new.env(parent = globalenv())
  expect_warning(
    invisible(utils::capture.output(
      suppressMessages(source(.iv_repl_path("hpc", "merge.R"), local = merge_env)),
      type = "output"
    )),
    "expected cell files"
  )

  expect_true(file.exists(file.path(out_dir, "sim_iv_raw.rds")))
  expect_true(file.exists(file.path(out_dir, "sim_iv_table.csv")))
})
