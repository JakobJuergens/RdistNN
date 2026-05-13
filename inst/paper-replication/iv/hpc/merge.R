# inst/paper-replication/iv/hpc/merge.R
#
# Collect per-cell .rds files written by the SLURM array tasks, merge them
# into a single raw data.frame, and write the summary table.
#
# Usage:
#   Rscript merge.R
#   SIM_OUTDIR=/custom/path Rscript merge.R

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

GRID <- iv_paper_grid()
expected_cells <- iv_expected_cells(GRID)

out_dir <- Sys.getenv("SIM_OUTDIR", unset = "results")
cells_dir <- file.path(out_dir, "cells")

files <- sort(list.files(
  cells_dir,
  pattern = "^cell_[0-9]+\\.rds$",
  full.names = TRUE
))

if (length(files) == 0L) {
  stop(
    "No cell_*.rds files found in: ", cells_dir,
    "\nCheck that the SLURM array tasks completed successfully.",
    call. = FALSE
  )
}

if (length(files) < expected_cells) {
  warning(
    sprintf("%d of %d expected cell files found. ", length(files), expected_cells),
    "Results will be partial.",
    call. = FALSE
  )
}

message(sprintf("Merging %d cell file(s) from %s ...", length(files), cells_dir))
raw <- do.call(rbind, lapply(files, readRDS))

tab <- summarise_iv_results(raw)

if (any(tab$n_failed > 0L)) {
  warning(
    sprintf("%d replicate(s) failed across all cells. ", sum(tab$n_failed)),
    "See n_failed column in the summary table.",
    call. = FALSE
  )
}

cat("\n=== IV-GMM simulation results ===\n")
print(tab, digits = 3)

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
saveRDS(raw, file.path(out_dir, "sim_iv_raw.rds"))
write.csv(tab, file.path(out_dir, "sim_iv_table.csv"), row.names = FALSE)
message("\nSaved: ", file.path(out_dir, "sim_iv_raw.rds"))
message("Saved: ", file.path(out_dir, "sim_iv_table.csv"))
