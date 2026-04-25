# inst/paper-replication/cate/hpc/merge.R
#
# Collect per-cell .rds files written by the SLURM array tasks, merge them
# into a single raw data.frame, and write the summary table.
#
# Run after all array tasks have completed:
#
#   Rscript merge.R
#   SIM_OUTDIR=/custom/path Rscript merge.R
#
# Or interactively:
#
#   source(system.file("paper-replication/cate/hpc/merge.R", package = "RdistNN"))

library(RdistNN)

out_dir <- Sys.getenv("SIM_OUTDIR", unset = "results")
cells_dir <- file.path(out_dir, "cells")

# Collect cell files.

files <- sort(list.files(
  cells_dir,
  pattern = "^cell_[0-9]+\\.rds$",
  full.names = TRUE
))

if (length(files) == 0L) {
  stop(
    "No cell_*.rds files found in: ", cells_dir,
    "\nCheck that the SLURM array tasks completed successfully."
  )
}

# Warn about any missing cells in the default main production grid.
expected_cells <- 16L  # 4 designs x 4 n_vals x 1 kx_val; update if grid changes
if (length(files) < expected_cells) {
  warning(
    sprintf("%d of %d expected cell files found.  ", length(files), expected_cells),
    "Results will be partial.",
    call. = FALSE
  )
}

message(sprintf("Merging %d cell file(s) from %s ...", length(files), cells_dir))
raw <- do.call(rbind, lapply(files, readRDS))

# Summary table.

tab <- summarise_sim_results(raw)

if (any(tab$n_failed > 0L)) {
  warning(
    sprintf("%d replicate(s) failed across all cells.  ", sum(tab$n_failed)),
    "See n_failed column in the summary table.",
    call. = FALSE
  )
}

cat("\n=== Simulation results ===\n")
print(tab, digits = 3)

# Save.

saveRDS(raw, file.path(out_dir, "sim_cate_raw.rds"))
write.csv(tab, file.path(out_dir, "sim_cate_table.csv"), row.names = FALSE)
message("\nSaved: ", file.path(out_dir, "sim_cate_raw.rds"))
message("Saved: ", file.path(out_dir, "sim_cate_table.csv"))
