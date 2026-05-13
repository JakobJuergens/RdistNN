# inst/paper-replication/iv/run.R
#
# Full paper replication script for the local IV-GMM simulation study.
#
# Assumes RdistNN is installed. Run from any working directory:
#
#   Rscript system.file("paper-replication/iv/run.R", package = "RdistNN")
#
# or interactively:
#
#   source(system.file("paper-replication/iv/run.R", package = "RdistNN"))
#
# Results are written to ./results/sim_iv_raw.rds and
# ./results/sim_iv_table.csv in the current working directory by default.
# Set SIM_OUTDIR to write results elsewhere.

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
  local_file <- file.path(.iv_script_dir(), "grid.R")
  if (file.exists(local_file)) {
    return(local_file)
  }

  installed_file <- system.file("paper-replication/iv/grid.R", package = "RdistNN")
  if (nzchar(installed_file)) {
    return(installed_file)
  }

  stop("Could not locate paper-replication/iv/grid.R.", call. = FALSE)
}

source(.iv_grid_file())

PAPER_GRID <- iv_paper_grid()
PAPER_CELLS <- iv_paper_cells(PAPER_GRID)

cat("Available methods:", paste(sim_iv_methods(), collapse = ", "), "\n\n")
cat("=== Running IV-GMM simulation ===\n")
cat("Cells:", nrow(PAPER_CELLS), "\n")
cat("Replications per cell:", PAPER_GRID$n_reps, "\n")
cat("Methods:", paste(PAPER_GRID$methods, collapse = ", "), "\n\n")

# Parallelise across reps by setting, for example:
#   future::plan(future::multisession)
raw <- sim_iv_run(
  designs = PAPER_GRID$designs,
  n_vals = PAPER_GRID$n_vals,
  kx_vals = PAPER_GRID$kx_vals,
  q_vals = PAPER_GRID$q_vals,
  n_reps = PAPER_GRID$n_reps,
  master_seed = PAPER_GRID$master_seed,
  k = PAPER_GRID$k,
  s_fn = PAPER_GRID$s_fn,
  methods = PAPER_GRID$methods,
  alpha = PAPER_GRID$alpha,
  sigma = PAPER_GRID$sigma,
  clip = PAPER_GRID$clip,
  n_folds = PAPER_GRID$n_folds,
  verbose = TRUE
)

table_iv <- summarise_iv_results(raw)

if (any(table_iv$n_failed > 0L)) {
  warning(
    sprintf("%d replicate(s) failed across all cells. ", sum(table_iv$n_failed)),
    "See n_failed column in the summary table.",
    call. = FALSE
  )
}

cat("\n=== IV-GMM simulation results ===\n")
print(table_iv, digits = 3)

out_dir <- Sys.getenv("SIM_OUTDIR", unset = "")
if (!nzchar(out_dir)) {
  out_dir <- file.path(getwd(), "results")
}
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

saveRDS(raw, file.path(out_dir, "sim_iv_raw.rds"))
write.csv(table_iv, file.path(out_dir, "sim_iv_table.csv"), row.names = FALSE)

cat("\nResults saved to", out_dir, "\n")
