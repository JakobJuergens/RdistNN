# Fit the declared SIPP eligibility CATE analysis. The default affine family
# compares literal all-tree, single-fit/OOB and held-out DNN AIPW predictions.
# Use --method-family=legacy for the historical OOB/five-fold/forest family.
# Every saved run has its own directory; historical flat caches are untouched.
# See README.md for economic specifications and complete command examples.
.script_path <- sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE)[1L])
.script_dir <- dirname(normalizePath(.script_path))
source(file.path(.script_dir, "01_prepare_pension.R"))
source(file.path(.script_dir, "R", "estimation_helpers.R"))
require_sipp_401k_fit_packages(c("hdm", "grf", "RdistNN", "digest", "jsonlite"))
fit_options <- parse_sipp_401k_fit_args(commandArgs(trailingOnly = TRUE))
config <- fit_options$config
sipp_401k <- apply_sipp_401k_fit_options(prepare_sipp_401k(fit_options$spec), fit_options)
rdistnn_source <- detect_sipp_401k_rdistnn_source(script_dir = .script_dir)
provenance <- capture_sipp_401k_provenance(.script_dir)
run_info <- list(
  r_version = R.version.string, hdm_version = sipp_401k$source$hdm_version,
  grf_version = as.character(utils::packageVersion("grf")),
  rdistnn_version = rdistnn_source$installed_version, rdistnn_source = rdistnn_source,
  executable_provenance = provenance,
  timestamp = format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"),
  smoke = fit_options$smoke, sample_n = fit_options$sample_n,
  outcome_scales = names(sipp_401k$outcome_cols),
  command_args = commandArgs(trailingOnly = TRUE))
run <- NULL
if (fit_options$write_cache) {
  cache_root <- if (is.null(fit_options$cache_dir)) file.path(.script_dir, "cache") else fit_options$cache_dir
  run <- start_sipp_401k_run(cache_root, fit_options$run_id,
    sipp_401k$specification, config, provenance)
  cat("Run directory:", run$directory, "\n")
}
tryCatch({
  set.seed(config$seed)
  sipp_401k_fits <- fit_sipp_401k_profiles(sipp_401k, config, run_info)
  cat("\nSIPP eligibility profile estimates; family:", config$method_family, "\n")
  print(sipp_401k_fits$estimates, row.names = FALSE)
  if (!is.null(run)) {
    path <- complete_sipp_401k_run(run, sipp_401k_fits)
    cat("Cached fit object:", path, "\n")
  } else cat("Cache writing skipped by --no-cache.\n")
}, error = function(error) {
  if (!is.null(run)) {
    manifest <- readRDS(file.path(run$directory, "manifest.rds"))
    if (manifest$status == "running") fail_sipp_401k_run(run, conditionMessage(error))
  }
  stop(conditionMessage(error), call. = FALSE)
})
