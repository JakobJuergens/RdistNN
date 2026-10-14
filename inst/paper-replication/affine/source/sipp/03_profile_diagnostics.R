# Summarize local support and estimator-agreement diagnostics for the SIPP 401(k)
# profile CATE fit cache.
#
# This script consumes the local cache written by `02_fit_cate_profiles.R` and
# writes a diagnostic RDS object under `cache/` by default. It does not export
# manuscript tables or figures.
#
# Usage:
#   Rscript Codes/Applications/SIPP_401k/03_profile_diagnostics.R
#
# Useful overrides:
#   --fit-cache=<path>                 Input fit cache
#   --output=<path>                    Output diagnostics cache
#   --cache-dir=<path>                 Directory containing/writing default caches
#   --no-cache                         Print only, do not write diagnostics cache
#   --min-ess=<number>                 Minimum local ESS screening threshold
#   --min-arm-mass=<number>            Minimum treated/untreated local mass
#   --max-top5-weight-share=<number>   Maximum top-five weight concentration
#   --propensity-lower=<number>        Lower propensity-screen boundary
#   --propensity-upper=<number>        Upper propensity-screen boundary
#   --max-se-ratio=<number>            DNN/causal-forest SE ratio screen
#   --max-tau-gap-se-units=<number>    DNN/causal-forest point gap screen

.script_path <- if (interactive()) {
  attr(sys.call(-1L), "srcfile")$filename
} else {
  args <- commandArgs(trailingOnly = FALSE)
  sub("--file=", "", grep("--file=", args, value = TRUE)[1L])
}
.script_dir <- tryCatch(
  dirname(normalizePath(.script_path)),
  error = function(e) getwd()
)

source(file.path(.script_dir, "R", "estimation_helpers.R"))

diag_options <- parse_sipp_401k_diagnostic_args(
  commandArgs(trailingOnly = TRUE),
  script_dir = .script_dir
)
if (diag_options$write_cache && file.exists(diag_options$output))
  stop("Diagnostic output already exists. Choose a new --output or use --no-cache.",
       call. = FALSE)

if (!file.exists(diag_options$fit_cache)) {
  stop(
    "Missing SIPP 401(k) fit cache: ",
    diag_options$fit_cache,
    "\nRun 02_fit_cate_profiles.R first or pass --fit-cache=<path>.",
    call. = FALSE
  )
}

sipp_401k_fits <- readRDS(diag_options$fit_cache)
sipp_401k_diagnostics <- make_sipp_401k_profile_diagnostics(
  sipp_401k_fits,
  config = diag_options$config
)

cat("\nSIPP 401(k) profile diagnostics\n")
cat("--------------------------------\n")
cat("Fit cache:", diag_options$fit_cache, "\n")
cat("Support thresholds:\n")
print(
  data.frame(
    threshold = names(diag_options$config),
    value = unlist(diag_options$config, use.names = FALSE),
    row.names = NULL
  ),
  row.names = FALSE
)

cat("\nMechanical support checks:\n")
profile_screen <- sipp_401k_diagnostics$profile_screen
print(
  profile_screen[
    ,
    c(
      "outcome_scale", "profile_id", "n_methods", "flagged_methods",
      "passes_mechanical_checks"
    ),
    drop = FALSE
  ],
  row.names = FALSE
)

cat("\nHistorical comparator descriptive flags:\n")
se_rows <- sipp_401k_diagnostics$competitor_comparison
if (nrow(se_rows) > 0L) {
print(
  se_rows[
    ,
    c(
      "outcome_scale", "profile_id", "dnn_method", "competitor_method", "se_ratio",
      "tau_gap_se_units", "se_disagreement", "large_tau_gap"
    ),
    drop = FALSE
  ],
  row.names = FALSE
)
} else cat("No historical comparator in this method family.\n")

if (!is.null(sipp_401k_diagnostics$frozen)) {
  cat("\nFrozen-score information diagnostics:\n")
  print(sipp_401k_diagnostics$frozen$information[, c(
    "outcome_scale", "profile_id", "method", "jk_top1_sq_share",
    "ipw_treated_kish", "ipw_control_kish", "local_clipped_mass")], row.names = FALSE)
}

if (isTRUE(diag_options$write_cache)) {
  dir.create(dirname(diag_options$output), recursive = TRUE, showWarnings = FALSE)
  saveRDS(sipp_401k_diagnostics, diag_options$output)
  cat("\nCached diagnostic object:", diag_options$output, "\n")
} else {
  cat("\nDiagnostic cache writing skipped by --no-cache.\n")
}
