#!/usr/bin/env Rscript
script <- sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE))
if (length(script) != 1L) stop("Run with Rscript.", call. = FALSE)
bundle <- dirname(normalizePath(script))
source(file.path(bundle, "helpers.R"))
opts <- replication_args(commandArgs(TRUE), c("output-dir", "help"))
if (opts$help) {
  cat(
    "Usage: render_evidence.R --output-dir=NEW_DIR\n",
    "Verify the saved evidence and regenerate 21 current manuscript tables and six PDF figures. No data generation or estimation.\n"
  )
} else {
  replication_render(bundle, opts$output_dir)
  cat("Verified current affine manuscript evidence:", opts$output_dir, "\n")
}
