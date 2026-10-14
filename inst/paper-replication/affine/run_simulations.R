#!/usr/bin/env Rscript
script <- sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE))
if (length(script) != 1L) stop("Run with Rscript.", call. = FALSE)
bundle <- dirname(normalizePath(script))
source(file.path(bundle, "helpers.R"))
opts <- replication_args(commandArgs(TRUE))
if (opts$help) {
  cat(
    "Usage: run_simulations.R --mode=smoke|production --output-dir=NEW_DIR [--profiles=dnn,kernel_comparison,M3,settings_comparison] [--resume]\n",
    "Smoke uses two preflight draw families at n=128; it is validation only.\n",
    "Production uses all four accepted profiles and 10,000 draws per cell, with Linux process measurement.\n"
  )
} else {
  result <- replication_simulations(bundle, opts$output_dir, opts$mode, opts$profiles, opts$resume)
  cat(result$status, ":", opts$output_dir, "\n")
}
