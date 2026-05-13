# inst/paper-replication/iv/hpc/task.R
#
# SLURM array task: runs all n_reps for one (design, n, kx, q) cell and saves
# the raw per-replicate results as an .rds file.
#
# Cell index is read from SLURM_ARRAY_TASK_ID or from the first command-line
# argument:
#
#   Rscript task.R 3
#
# Output is written to $SIM_OUTDIR/cells/cell_<idx>.rds.

library(RdistNN)

.iv_script_dir <- function() {
  file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (length(file_arg) > 0L) {
    return(dirname(normalizePath(sub("^--file=", "", file_arg[1L]), mustWork = FALSE)))
  }

  frames <- sys.frames()
  files <- vapply(frames, function(frame) {
    file <- frame$ofile
    if (is.null(file)) NA_character_ else file
  }, character(1L))
  files <- files[!is.na(files)]
  if (length(files) > 0L) {
    return(dirname(normalizePath(files[length(files)], mustWork = FALSE)))
  }

  getwd()
}

.iv_grid_file <- function() {
  local_file <- file.path(.iv_script_dir(), "..", "grid.R")
  if (file.exists(local_file)) {
    return(normalizePath(local_file, mustWork = FALSE))
  }

  installed_file <- system.file("paper-replication/iv/grid.R", package = "RdistNN")
  if (nzchar(installed_file)) {
    return(installed_file)
  }

  stop("Could not locate paper-replication/iv/grid.R.", call. = FALSE)
}

source(.iv_grid_file())

cell_idx <- suppressWarnings(as.integer(
  Sys.getenv("SLURM_ARRAY_TASK_ID", unset = commandArgs(trailingOnly = TRUE)[1L])
))
if (is.na(cell_idx) || cell_idx < 1L) {
  stop(
    "Could not determine cell index. Set SLURM_ARRAY_TASK_ID ",
    "or pass it as the first command-line argument.",
    call. = FALSE
  )
}

GRID <- iv_paper_grid()
cells <- iv_paper_cells(GRID)
n_cells <- nrow(cells)

if (cell_idx > n_cells) {
  stop(
    sprintf("cell_idx=%d exceeds the number of cells (%d).", cell_idx, n_cells),
    call. = FALSE
  )
}

cell <- cells[cell_idx, , drop = FALSE]
d <- cell$design
n <- cell$n
kx <- cell$kx
q <- cell$q
s <- GRID$s_fn(n, kx)

ncores <- suppressWarnings(as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", unset = "1")))
if (is.na(ncores) || ncores < 1L) {
  ncores <- 1L
}

if (ncores > 1L &&
    requireNamespace("future", quietly = TRUE) &&
    requireNamespace("future.apply", quietly = TRUE)) {
  future::plan(future::multisession, workers = ncores)
  on.exit(future::plan(future::sequential), add = TRUE)
  message(sprintf("[task %02d] multisession workers=%d", cell_idx, ncores))
} else {
  message(sprintf("[task %02d] sequential", cell_idx))
}

message(sprintf(
  "[task %02d] design=%-32s n=%4d kx=%d q=%d s=%d",
  cell_idx, d, n, kx, q, s
))

raw <- sim_iv_run(
  designs = d,
  n_vals = n,
  kx_vals = kx,
  q_vals = q,
  n_reps = GRID$n_reps,
  master_seed = GRID$master_seed,
  k = GRID$k,
  s_fn = GRID$s_fn,
  methods = GRID$methods,
  alpha = GRID$alpha,
  sigma = GRID$sigma,
  clip = GRID$clip,
  n_folds = GRID$n_folds,
  verbose = FALSE
)

out_dir <- file.path(Sys.getenv("SIM_OUTDIR", unset = "results"), "cells")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
out_file <- file.path(out_dir, sprintf("cell_%02d.rds", cell_idx))
saveRDS(raw, out_file)
message(sprintf("[task %02d] saved -> %s", cell_idx, out_file))
