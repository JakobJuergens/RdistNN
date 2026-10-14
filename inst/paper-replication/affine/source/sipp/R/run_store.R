# Run identity, actual executable provenance and immutable terminal records.
sipp_401k_fit_fingerprint <- function(fit) {
  digest::digest(fit[c("schema_version", "estimates", "analysis_data", "dnn_weights",
    "signals", "nuisances", "tuning_grids", "profile_diagnostics", "profiles",
    "specification", "method_metadata", "config", "source", "run_info")], algo = "sha256")
}

sipp_401k_file_hash <- function(path) {
  if (!file.exists(path)) stop("Cannot hash missing file: ", path, call. = FALSE)
  digest::digest(file = path, algo = "sha256", serialize = FALSE)
}

capture_sipp_401k_provenance <- function(app_dir) {
  packages <- lapply(c("RdistNN", "grf", "hdm"), function(package) {
    requireNamespace(package, quietly = TRUE)
    path <- find.package(package)
    files <- list.files(path, recursive = TRUE, full.names = TRUE)
    files <- files[grepl("/(R|libs|data)/|/(DESCRIPTION|NAMESPACE)$", gsub("\\\\", "/", files))]
    names(files) <- substring(gsub("\\\\", "/", files), nchar(gsub("\\\\", "/", path)) + 2L)
    list(version = as.character(utils::packageVersion(package)),
         installed_path = normalizePath(path, winslash = "/"),
         files = lapply(files, sipp_401k_file_hash))
  })
  names(packages) <- c("RdistNN", "grf", "hdm")
  files <- list.files(app_dir, pattern = "\\.R$", full.names = TRUE)
  files <- c(files, list.files(file.path(app_dir, "R"), pattern = "\\.R$", full.names = TRUE))
  names(files) <- substring(files, nchar(app_dir) + 2L)
  list(captured_utc = format(Sys.time(), tz = "UTC", usetz = TRUE),
       R = R.version.string, platform = R.version$platform,
       packages = packages, application_files = lapply(files, sipp_401k_file_hash),
       loaded_libraries = lapply(getLoadedDLLs(), function(dll) dll[["path"]]))
}

write_sipp_401k_manifest <- function(directory, manifest) {
  saveRDS(manifest, file.path(directory, "manifest.rds"))
  jsonlite::write_json(manifest, file.path(directory, "manifest.json"),
                       pretty = TRUE, auto_unbox = TRUE, null = "null", na = "null", digits = NA)
}

start_sipp_401k_run <- function(cache_root, run_id = NULL, specification, config, provenance) {
  if (is.null(run_id)) run_id <- paste0(format(Sys.time(), "%Y%m%dT%H%M%S", tz = "UTC"),
    "-", substr(digest::digest(list(specification, config), algo = "sha256"), 1L, 10L))
  if (length(run_id) != 1L || is.na(run_id) || nchar(run_id) > 120L ||
      !grepl("^[A-Za-z0-9][A-Za-z0-9_.-]*$", run_id))
    stop("Run identifier must be a safe single directory name.", call. = FALSE)
  parent <- file.path(cache_root, "runs")
  dir.create(parent, recursive = TRUE, showWarnings = FALSE)
  directory <- file.path(parent, run_id)
  if (!dir.create(directory, showWarnings = FALSE))
    stop("Run already exists or cannot be created: ", directory, call. = FALSE)
  manifest <- list(schema_version = 2L, run_id = run_id, status = "running",
    started_utc = format(Sys.time(), tz = "UTC", usetz = TRUE),
    specification = specification, config = config, provenance = provenance)
  write_sipp_401k_manifest(directory, manifest)
  list(run_id = run_id, directory = normalizePath(directory, winslash = "/"))
}

complete_sipp_401k_run <- function(run, fit) {
  lock <- file.path(run$directory, ".completion")
  if (!dir.create(lock, showWarnings = FALSE))
    stop("Run completion is already in progress.", call. = FALSE)
  on.exit(unlink(lock, recursive = TRUE), add = TRUE)
  manifest <- readRDS(file.path(run$directory, "manifest.rds"))
  if (manifest$status != "running") stop("Only a running run can complete.", call. = FALSE)
  intended <- file.path(run$directory, c("models.rds", "02_fit_cate_profiles.rds"))
  if (any(file.exists(intended)))
    stop("Run output already exists; use a new run.", call. = FALSE)
  fit$run_info$run_id <- run$run_id
  if (!is.null(fit$models)) {
    path <- file.path(run$directory, "models.rds")
    saveRDS(fit$models, path)
    fit$model_files <- list(path = "models.rds", sha256 = sipp_401k_file_hash(path))
    fit$models <- NULL
  }
  path <- file.path(run$directory, "02_fit_cate_profiles.rds")
  saveRDS(fit, path)
  manifest$status <- "complete"
  manifest$completed_utc <- format(Sys.time(), tz = "UTC", usetz = TRUE)
  manifest$selected_row_sha256 <- digest::digest(fit$analysis_data$row_id, algo = "sha256")
  manifest$analysis_data_sha256 <- digest::digest(fit$analysis_data, algo = "sha256")
  manifest$outputs <- list(fit = list(path = basename(path), sha256 = sipp_401k_file_hash(path)),
                          models = fit$model_files)
  write_sipp_401k_manifest(run$directory, manifest)
  path
}

fail_sipp_401k_run <- function(run, reason) {
  manifest <- readRDS(file.path(run$directory, "manifest.rds"))
  if (manifest$status != "running") stop("Cannot change a terminal run.", call. = FALSE)
  manifest$status <- "failed"; manifest$error <- reason
  manifest$finished_utc <- format(Sys.time(), tz = "UTC", usetz = TRUE)
  write_sipp_401k_manifest(run$directory, manifest)
  invisible(run)
}
