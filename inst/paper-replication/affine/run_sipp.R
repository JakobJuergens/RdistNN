#!/usr/bin/env Rscript
script <- sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE))
if (length(script) != 1L) stop("Run with Rscript.", call. = FALSE)
bundle <- dirname(normalizePath(script))
source(file.path(bundle, "helpers.R"))
opts <- replication_args(commandArgs(TRUE), c("mode", "output-dir", "help"))
if (opts$help) {
  cat(
    "Usage: run_sipp.R --mode=smoke|production --output-dir=NEW_DIR\n",
    "Smoke uses 300 households, 50 trees and two folds; exports are labeled validation.\n",
    "Production fits raw-dollar outcomes for all 9,915 households, 2,000 trees and five folds.\n"
  )
} else {
  replication_sipp(bundle, opts$output_dir, opts$mode)
  cat("Completed", opts$mode, "SIPP workflow:", opts$output_dir, "\n")
}
