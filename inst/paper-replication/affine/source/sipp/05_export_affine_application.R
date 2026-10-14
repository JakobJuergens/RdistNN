# Export a declared affine SIPP run without fitting. The historical exporter
# 04 remains separate. Input and a new output directory must be explicit.
.script <- sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE)[1L])
source(file.path(dirname(normalizePath(.script)), "R", "estimation_helpers.R"))
args <- commandArgs(trailingOnly = TRUE)
options <- list(fit_cache = NULL, diagnostics_cache = NULL, output_dir = NULL,
                profile_ids = NULL, allow_validation = FALSE)
for (arg in args) {
  if (arg == "--allow-validation") { options$allow_validation <- TRUE; next }
  key <- sub("=.*$", "", arg); value <- sub("^[^=]*=", "", arg)
  field <- switch(key, "--fit-cache" = "fit_cache", "--diagnostics-cache" = "diagnostics_cache",
                  "--output-dir" = "output_dir", "--profiles" = "profile_ids",
                  stop("Unknown affine export argument: ", arg, call. = FALSE))
  options[[field]] <- if (field == "profile_ids") strsplit(value, ",", fixed = TRUE)[[1L]] else value
}
if (is.null(options$fit_cache) || is.null(options$output_dir))
  stop("Provide --fit-cache and a new --output-dir.", call. = FALSE)
fit <- readRDS(options$fit_cache)
diagnostics <- if (is.null(options$diagnostics_cache)) make_sipp_401k_profile_diagnostics(fit) else
  readRDS(options$diagnostics_cache)
export_sipp_401k_affine(fit, diagnostics, options$output_dir, options$profile_ids, options$allow_validation)
cat("Exported affine SIPP evidence:", options$output_dir, "\n")
