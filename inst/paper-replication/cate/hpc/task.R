# inst/paper-replication/cate/hpc/task.R
#
# SLURM array task: runs all n_reps for one (design, n, kx) cell and saves
# the raw per-replicate results as an .rds file.
#
# Cell index is read from SLURM_ARRAY_TASK_ID (set by SLURM) or from the
# first command-line argument (for interactive or non-SLURM use):
#
#   Rscript task.R 3          # cell index 3, interactive / non-SLURM
#
# Within-task rep parallelism uses all cores allocated by SLURM
# (--cpus-per-task).  Requires future + future.apply; falls back to
# sequential if they are not installed.
#
# Output is written to $SIM_OUTDIR/cells/cell_<idx>.rds.
# Set the SIM_OUTDIR environment variable in submit.sh (default: "results").

library(RdistNN)

# ── Cell index ────────────────────────────────────────────────────────────────

cell_idx <- suppressWarnings(as.integer(
  Sys.getenv("SLURM_ARRAY_TASK_ID",
             unset = commandArgs(trailingOnly = TRUE)[1L])
))
if (is.na(cell_idx) || cell_idx < 1L)
  stop("Could not determine cell index.  ",
       "Set SLURM_ARRAY_TASK_ID or pass it as the first command-line argument.")

# ── Paper grid ────────────────────────────────────────────────────────────────
# Must match merge.R.  Edit both files if you change the grid.

GRID <- list(
  designs     = c("dgp1_linear", "dgp2_nonlinear", "dgp3_sparse_overlap",
                  "dgp4_clustered", "dgp5_correlated"),
  n_vals      = c(500L, 1000L, 2000L),
  kx_vals     = 2L,
  n_reps      = 500L,
  master_seed = 42L,
  k           = 5L,
  s_fn        = function(n, kx) floor(n^((kx + 1L) / (kx + 2L))),
  # Must match run.R and merge.R.  dnn_oracle_mahal is especially informative
  # for dgp5_correlated (shape adaptivity).
  methods     = c("dnn_fullsample", "dnn_oracle_mahal", "dnn_crossfit",
                  "local_constant", "ols_interact"),
  alpha       = 0.05,
  sigma       = 0.5,
  clip        = 1e-3
)

cells <- expand.grid(
  design = GRID$designs,
  n      = GRID$n_vals,
  kx     = GRID$kx_vals,
  stringsAsFactors = FALSE
)
n_cells <- nrow(cells)

if (cell_idx > n_cells)
  stop(sprintf("cell_idx=%d exceeds the number of cells (%d).",
               cell_idx, n_cells))

d  <- cells$design[cell_idx]
n  <- cells$n[cell_idx]
kx <- cells$kx[cell_idx]
s  <- GRID$s_fn(n, kx)

# ── Within-task parallelism ───────────────────────────────────────────────────
# Uses cores allocated by SLURM (--cpus-per-task).  On non-SLURM systems
# defaults to 1 (sequential), which is also the safe fallback.

ncores <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", unset = "1"))

if (ncores > 1L &&
    requireNamespace("future",       quietly = TRUE) &&
    requireNamespace("future.apply", quietly = TRUE)) {
  future::plan(future::multisession, workers = ncores)
  on.exit(future::plan(future::sequential), add = TRUE)
  message(sprintf("[task %02d]  multisession  workers=%d", cell_idx, ncores))
} else {
  message(sprintf("[task %02d]  sequential", cell_idx))
}

# ── Run ───────────────────────────────────────────────────────────────────────

message(sprintf("[task %02d]  design=%-22s  n=%4d  kx=%d  s=%d",
                cell_idx, d, n, kx, s))

raw <- sim_cate_run(
  designs     = d,
  n_vals      = n,
  kx_vals     = kx,
  n_reps      = GRID$n_reps,
  master_seed = GRID$master_seed,
  k           = GRID$k,
  s_fn        = GRID$s_fn,
  methods     = GRID$methods,
  alpha       = GRID$alpha,
  sigma       = GRID$sigma,
  clip        = GRID$clip,
  verbose     = FALSE
)

# ── Save ──────────────────────────────────────────────────────────────────────

out_dir <- file.path(Sys.getenv("SIM_OUTDIR", unset = "results"), "cells")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
out_file <- file.path(out_dir, sprintf("cell_%02d.rds", cell_idx))
saveRDS(raw, out_file)
message(sprintf("[task %02d]  saved -> %s", cell_idx, out_file))
