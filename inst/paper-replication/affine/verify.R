#!/usr/bin/env Rscript
script <- sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE))
if (length(script) != 1L) stop("Run with Rscript.", call. = FALSE)
bundle <- dirname(normalizePath(script))
source(file.path(bundle, "helpers.R"))
replication_verify(bundle)
cat("Verified pinned affine source, production configurations and evidence archive.\n")
