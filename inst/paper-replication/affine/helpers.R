# Portable launchers around the paper's unchanged scientific source.
replication_require <- function(packages) {
  missing <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing)) stop("Install replication dependencies: ", paste(missing, collapse = ", "), call. = FALSE)
  invisible(TRUE)
}

replication_manifest <- function(bundle) {
  replication_require(c("digest", "jsonlite"))
  jsonlite::read_json(file.path(bundle, "manifest.json"), simplifyVector = TRUE)
}

replication_verify <- function(bundle) {
  manifest <- replication_manifest(bundle)
  for (name in names(manifest$files)) {
    path <- file.path(bundle, name)
    if (!file.exists(path) || !identical(digest::digest(file = path, algo = "sha256"), manifest$files[[name]])) {
      stop("Replication hash mismatch: ", name, call. = FALSE)
    }
  }
  invisible(TRUE)
}

replication_args <- function(args, allowed = c("mode", "output-dir", "profiles", "resume", "help")) {
  out <- list(mode = "smoke", output_dir = NULL, profiles = NULL, resume = FALSE, help = FALSE)
  seen <- character()
  for (arg in args) {
    key <- sub("=.*$", "", sub("^--", "", arg))
    if (!startsWith(arg, "--") || !key %in% allowed) stop("Unknown replication argument: ", arg, call. = FALSE)
    if (key %in% seen) stop("Repeated replication argument: --", key, call. = FALSE)
    seen <- c(seen, key)
    if (key %in% c("resume", "help")) {
      if (arg != paste0("--", key)) stop("Use --", key, " without a value.", call. = FALSE)
      out[[key]] <- TRUE
    } else {
      if (!grepl("=.+$", arg)) stop("Supply a value for --", key, call. = FALSE)
      out[[gsub("-", "_", key)]] <- sub("^[^=]*=", "", arg)
    }
  }
  if (!out$mode %in% c("smoke", "production")) stop("`mode` must be smoke or production.", call. = FALSE)
  if (!is.null(out$profiles)) out$profiles <- strsplit(out$profiles, ",", fixed = TRUE)[[1L]]
  out
}

replication_new_output <- function(path) {
  if (is.null(path) || length(path) != 1L || is.na(path) || !nzchar(path)) {
    stop("Supply --output-dir with a new writable directory.", call. = FALSE)
  }
  path <- normalizePath(path, winslash = "/", mustWork = FALSE)
  if (file.exists(path) || dir.exists(path)) {
    stop("Output directory already exists; choose a new directory.", call. = FALSE)
  }
  if (!dir.create(path, recursive = TRUE)) stop("Cannot create output directory: ", path, call. = FALSE)
  normalizePath(path, winslash = "/")
}

replication_source <- function(bundle) {
  env <- new.env(parent = globalenv())
  directory <- file.path(bundle, "source", "sim")
  for (name in c(
    "affine_config.R", "affine_models.R", "affine_series.R", "affine_local.R",
    "affine_kernel.R", "affine_diagnostics.R", "affine_summary.R", "affine_runner.R"
  )) {
    sys.source(file.path(directory, name), envir = env)
  }
  env
}

replication_configs <- function(bundle, mode = c("smoke", "production")) {
  mode <- match.arg(mode)
  profiles <- c("dnn", "kernel_comparison", "M3", "settings_comparison")
  configs <- setNames(lapply(profiles, function(profile) {
    readRDS(file.path(bundle, "configs", profile, "config.rds"))
  }), profiles)
  if (mode == "production") {
    return(configs)
  }
  code <- replication_source(bundle)
  lapply(configs, function(cfg) {
    cfg$stage <- NULL
    cfg$cells <- data.frame(
      n = 128L, order_mult = 1,
      design_id = if (cfg$study == "settings_comparison") "C0H" else "C0"
    )
    cfg$replications <- 2L
    cfg$control_reps <- 0L
    cfg$intensive_reps <- 0L
    cfg$workers <- 1L
    cfg$chunk_size <- 1L
    cfg$cpu_budget <- 300
    cfg$campaign_id <- paste0("VALIDATION_", cfg$campaign_id)
    code$affine_config("preflight", cfg)
  })
}

replication_with_env <- function(values, expr) {
  old <- Sys.getenv(names(values), unset = NA_character_)
  on.exit(for (name in names(old)) {
    if (is.na(old[[name]])) Sys.unsetenv(name) else do.call(Sys.setenv, setNames(list(old[[name]]), name))
  }, add = TRUE)
  do.call(Sys.setenv, as.list(values))
  force(expr)
}

replication_run_r <- function(script, args = character(), log = "") {
  status <- system2(file.path(R.home("bin"), "Rscript"), c(shQuote(script), shQuote(args)),
    stdout = if (nzchar(log)) log else "", stderr = if (nzchar(log)) log else ""
  )
  if (status != 0L) {
    stop("Replication script failed: ", basename(script),
      if (nzchar(log)) paste0("; see ", log),
      call. = FALSE
    )
  }
  invisible(TRUE)
}

replication_restore_config <- function(frozen, runtime) {
  if (!identical(readRDS(frozen), readRDS(runtime))) {
    stop("Runtime configuration settings differ from the frozen production configuration.", call. = FALSE)
  }
  # RDS embeds its writer's R version. Restore byte identity for the upstream
  # audit only after proving that the scientific configuration is unchanged.
  if (!file.copy(frozen, runtime, overwrite = TRUE)) {
    stop("Cannot restore frozen configuration bytes: ", runtime, call. = FALSE)
  }
  invisible(TRUE)
}

replication_stage_valid <- function(path, stage, audit_file = NULL) {
  tryCatch(
    {
      file <- if (stage == "audit") "audit.json" else "analysis_manifest.json"
      if (!file.exists(file.path(path, file))) {
        return(FALSE)
      }
      manifest <- jsonlite::read_json(file.path(path, file), simplifyVector = TRUE)
      expected <- if (stage == "audit") "completed_and_audited_production" else "analysis_of_audited_saved_production"
      if (!identical(manifest$status, expected)) {
        return(FALSE)
      }
      if (stage == "audit") {
        if (manifest$total_rows != 2516400 || manifest$common_independent_draw_families != 10000 ||
          !file.exists(file.path(path, "SHA256SUMS"))) {
          return(FALSE)
        }
        entries <- readLines(file.path(path, "SHA256SUMS"), warn = FALSE)
        if (!length(entries) || any(!grepl("^[0-9a-f]{64}  .+", entries))) {
          return(FALSE)
        }
        hashes <- substr(entries, 1L, 64L)
        names(hashes) <- sub("^[0-9a-f]{64}  ", "", entries)
      } else {
        if (isTRUE(manifest$estimation_performed)) {
          return(FALSE)
        }
        if (!is.null(audit_file) && !identical(
          manifest$source_audit_sha256,
          digest::digest(file = audit_file, algo = "sha256")
        )) {
          return(FALSE)
        }
        h <- manifest$data_hashes
        if (!is.data.frame(h) || !nrow(h) || !all(c("file", "sha256") %in% names(h))) {
          return(FALSE)
        }
        hashes <- setNames(h$sha256, h$file)
      }
      all(vapply(names(hashes), function(name) {
        file.exists(file.path(path, name)) &&
          identical(digest::digest(file = file.path(path, name), algo = "sha256"), hashes[[name]])
      }, logical(1)))
    },
    error = function(error) FALSE
  )
}

replication_stage <- function(destination, stage, run, audit_file = NULL) {
  if (replication_stage_valid(destination, stage, audit_file)) {
    return(invisible(TRUE))
  }
  if (file.exists(destination) || dir.exists(destination)) {
    preserved <- tempfile(paste0(basename(destination), ".interrupted-"), tmpdir = dirname(destination))
    if (!file.rename(destination, preserved)) stop("Cannot preserve interrupted stage: ", destination, call. = FALSE)
  }
  staging <- tempfile(paste0(basename(destination), ".staging-"), tmpdir = dirname(destination))
  run(staging)
  if (!replication_stage_valid(staging, stage, audit_file)) {
    stop("Stage completion receipt or output hashes are invalid: ", staging, call. = FALSE)
  }
  if (!file.rename(staging, destination)) stop("Cannot publish verified stage: ", destination, call. = FALSE)
  invisible(TRUE)
}

replication_simulations <- function(bundle, output, mode = c("smoke", "production"), profiles = NULL, resume = FALSE) {
  mode <- match.arg(mode)
  replication_require(c("RdistNN", "digest", "jsonlite", "future", "future.apply"))
  replication_verify(bundle)
  configs <- replication_configs(bundle, mode)
  if (is.null(profiles)) profiles <- names(configs)
  if (!length(profiles) || anyDuplicated(profiles) || any(!profiles %in% names(configs))) {
    stop("`profiles` must select dnn, kernel_comparison, M3 or settings_comparison.", call. = FALSE)
  }
  if (resume && mode != "production") stop("Resume is available for production campaigns only.", call. = FALSE)
  if (mode == "production" && (Sys.info()[["sysname"]] != "Linux" || !nzchar(Sys.which("python3")))) {
    stop(
      "Measured production launch requires Linux and python3; smoke runs and evidence rendering are portable.",
      call. = FALSE
    )
  }
  if (resume) {
    output <- normalizePath(output, winslash = "/", mustWork = TRUE)
    previous <- jsonlite::read_json(file.path(output, "replication.json"), simplifyVector = TRUE)
    if (previous$mode != mode || previous$bundle_manifest_sha256 !=
      digest::digest(file = file.path(bundle, "manifest.json"), algo = "sha256")) {
      stop("Resume bundle or mode mismatch.", call. = FALSE)
    }
  } else {
    output <- replication_new_output(output)
    dir.create(file.path(output, "configs"))
    file.copy(file.path(bundle, "configs", "proposal.json"), file.path(output, "configs", "proposal.json"))
    for (profile in names(configs)) {
      dir.create(file.path(output, "configs", profile))
      if (mode == "production") {
        file.copy(
          file.path(bundle, "configs", profile, "config.rds"),
          file.path(output, "configs", profile, "config.rds")
        )
      } else {
        saveRDS(configs[[profile]], file.path(output, "configs", profile, "config.rds"))
      }
    }
    jsonlite::write_json(list(
      status = "running", mode = mode,
      scope = if (mode == "smoke") {
        "validation_only; 128 observations; two preflight draw families"
      } else {
        "fresh replication of the accepted 28-cell scientific configuration"
      },
      bundle_manifest_sha256 = digest::digest(file = file.path(bundle, "manifest.json"), algo = "sha256"),
      paper_commit = replication_manifest(bundle)$paper_commit
    ), file.path(output, "replication.json"), auto_unbox = TRUE, pretty = TRUE)
  }
  source_dir <- file.path(bundle, "source", "sim")
  code <- replication_source(bundle)
  for (profile in profiles) {
    path <- file.path(output, "runs", profile)
    if (mode == "smoke") {
      replication_with_env(
        c(SIM_RDISTNN_SOURCE_DIR = "", TZ = "UTC"),
        code$affine_run(configs[[profile]], path, source_dir = source_dir)
      )
    } else {
      # The paper's sampled RSS/wall measurement tool requires Linux /proc.
      if (Sys.info()[["sysname"]] != "Linux" || !nzchar(Sys.which("python3"))) {
        stop(
          "Measured production launch requires Linux and python3; smoke runs and evidence rendering are portable.",
          call. = FALSE
        )
      }
      measure <- file.path(output, "runs", "measurements", paste0(profile, "_", format(Sys.time(), "%Y%m%dT%H%M%S")))
      spec <- jsonlite::read_json(file.path(bundle, "configs", "proposal.json"), simplifyVector = TRUE)
      settings <- spec$profiles[[profile]]
      values <- c(
        SIM_RDISTNN_SOURCE_DIR = "", SIM_AFFINE_CONFIG_RDS = file.path(output, "configs", profile, "config.rds"),
        SIM_AFFINE_OUTPUT_DIR = path,
        SIM_AFFINE_RESUME = if (resume && file.exists(file.path(path, "freeze.rds"))) "true" else "false",
        SIM_AFFINE_REPLICATIONS = "", SIM_AFFINE_RUN_NAME = "", SIM_AFFINE_WORKERS = "", SIM_AFFINE_CHUNK_SIZE = "",
        SIM_AFFINE_CPU_BUDGET = "", SIM_AFFINE_MAX_CHUNKS = "", OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1",
        MKL_NUM_THREADS = "1", VECLIB_MAXIMUM_THREADS = "1", NUMEXPR_NUM_THREADS = "1", TZ = "UTC"
      )
      args <- c(
        file.path(source_dir, "campaign_tools", "measure_process.py"), "--output", measure,
        "--wall-limit", as.character(settings$wall_limit_seconds), "--rss-limit-mib", as.character(spec$rss_limit_mib),
        "--", file.path(R.home("bin"), "Rscript"), file.path(source_dir, "run_affine_scores.R")
      )
      status <- replication_with_env(values, system2(Sys.which("python3"), shQuote(args)))
      if (status != 0L) {
        stop(
          "Production profile stopped: ", profile,
          "; inspect its measured process and sealed chunks before resuming.",
          call. = FALSE
        )
      }
    }
  }
  all_profiles_present <- all(vapply(names(configs), function(profile) {
    file.exists(file.path(output, "runs", profile, "manifest.json"))
  }, logical(1)))
  if (mode == "production" && all_profiles_present) {
    for (profile in names(configs)) {
      replication_restore_config(
        file.path(output, "configs", profile, "config.rds"), file.path(output, "runs", profile, "config.rds")
      )
    }
    audit <- file.path(output, "audited")
    replication_stage(audit, "audit", function(staging) {
      replication_run_r(file.path(source_dir, "campaign_tools", "audit_production.R"), c(output, staging))
    })
    analysis <- file.path(output, "analysis")
    replication_stage(analysis, "analysis", function(staging) {
      replication_run_r(file.path(source_dir, "campaign_tools", "analyze_production.R"), c(output, audit, staging))
    }, audit_file = file.path(audit, "audit.json"))
  }
  status <- if (mode == "smoke") {
    "validation_only"
  } else if (all_profiles_present) {
    "completed_and_audited_production"
  } else {
    "partial_production"
  }
  receipt <- jsonlite::read_json(file.path(output, "replication.json"), simplifyVector = TRUE)
  receipt$status <- status
  receipt$profiles_requested <- profiles
  jsonlite::write_json(receipt, file.path(output, "replication.json"), auto_unbox = TRUE, pretty = TRUE)
  invisible(list(status = status, mode = mode, output = output, profiles = profiles))
}

replication_extract <- function(bundle, output) {
  manifest <- replication_manifest(bundle)
  archive <- file.path(bundle, "evidence.zip")
  entries <- utils::unzip(archive, list = TRUE)$Name
  if (anyDuplicated(entries) || !setequal(entries, names(manifest$evidence_files)) ||
    any(grepl("(^/|^[A-Za-z]:|(^|/)\\.\\.(/|$))", entries))) {
    stop("Unexpected evidence archive paths.", call. = FALSE)
  }
  utils::unzip(archive, exdir = output)
  for (name in names(manifest$evidence_files)) {
    if (!identical(digest::digest(file = file.path(output, name), algo = "sha256"), manifest$evidence_files[[name]])) {
      stop("Evidence hash mismatch: ", name, call. = FALSE)
    }
  }
  invisible(TRUE)
}

replication_run_renderer <- function(script, args, paper_commit) {
  # Preserve the upstream renderer and supply its portable script path. Its Git
  # provenance refers to the pinned paper snapshot, rather than this output tree.
  env <- new.env(parent = globalenv())
  env$commandArgs <- function(trailingOnly = FALSE) {
    if (trailingOnly) args else c(paste0("--file=", script), args)
  }
  env$system2 <- function(command, ...) {
    if (identical(command, "git")) paper_commit else base::system2(command, ...)
  }
  env$normalizePath <- function(path, winslash = "/", mustWork = NA) {
    base::normalizePath(path, winslash = winslash, mustWork = mustWork)
  }
  env$sub <- function(pattern, replacement, x, ...) {
    # The upstream simulation renderer strips its repository prefix with sub().
    # Treat that prefix literally, including Windows drives and punctuation.
    if (identical(replacement, "") && exists("repo", envir = env, inherits = FALSE)) {
      prefix <- paste0(env$repo, "/")
      if (identical(pattern, paste0("^", prefix))) {
        keep <- startsWith(x, prefix)
        x[keep] <- substring(x[keep], nchar(prefix) + 1L)
        return(x)
      }
    }
    base::sub(pattern, replacement, x, ...)
  }
  source(script, local = env)
  invisible(TRUE)
}

replication_render <- function(bundle, output) {
  replication_require(c("digest", "jsonlite", "ggplot2", "patchwork"))
  replication_verify(bundle)
  output <- replication_new_output(output)
  manifest <- replication_manifest(bundle)
  for (stored in names(manifest$source_mapping)) {
    target <- file.path(output, manifest$source_mapping[[stored]])
    dir.create(dirname(target), recursive = TRUE, showWarnings = FALSE)
    if (!file.copy(file.path(bundle, stored), target)) {
      stop("Cannot copy the pinned source into the output directory: ", stored, call. = FALSE)
    }
  }
  replication_extract(bundle, output)
  scripts <- c(
    "Codes/Simulations/affine_scores/campaign_tools/render_manuscript_simulations.R",
    "Codes/Applications/SIPP_401k/07_render_affine_manuscript.R"
  )
  for (script in scripts) {
    replication_with_env(
      c(TZ = "UTC"), replication_run_renderer(file.path(output, script), character(), manifest$paper_commit)
    )
    replication_with_env(
      c(TZ = "UTC"), replication_run_renderer(file.path(output, script), "--check", manifest$paper_commit)
    )
  }
  for (name in names(manifest$reference_tables)) {
    actual <- readLines(file.path(output, name), warn = FALSE)
    expected <- readLines(file.path(output, "reference", name), warn = FALSE)
    if (!identical(actual, expected)) {
      stop("Manuscript table differs from the accepted reference: ", name, call. = FALSE)
    }
  }
  jsonlite::write_json(
    list(
      status = "verified", paper_commit = manifest$paper_commit,
      bundle_manifest_sha256 = digest::digest(file = file.path(bundle, "manifest.json"), algo = "sha256"),
      tables_matched = length(manifest$reference_tables), simulation_performed = FALSE, estimation_performed = FALSE,
      provenance_adapter = paste0(
        "unchanged renderers; portable script/relative paths supplied by launcher; ",
        "Git source identity taken from the pinned paper manifest"
      )
    ),
    file.path(output, "replication.json"),
    auto_unbox = TRUE, pretty = TRUE
  )
  invisible(TRUE)
}

replication_sipp <- function(bundle, output, mode = c("smoke", "production")) {
  mode <- match.arg(mode)
  replication_require(c("hdm", "grf", "RdistNN", "digest", "jsonlite"))
  replication_verify(bundle)
  output <- replication_new_output(output)
  app <- file.path(bundle, "source", "sipp")
  fit_args <- c(
    "--run-id=affine-current", "--outcomes=raw", "--method-family=affine",
    paste0("--cache-dir=", output)
  )
  if (mode == "smoke") fit_args <- c(fit_args, "--smoke", "--num-threads=1")
  fit_log <- file.path(output, "fit.log")
  replication_run_r(file.path(app, "02_fit_cate_profiles.R"), fit_args, fit_log)
  run <- file.path(output, "runs", "affine-current")
  fit_cache <- file.path(run, "02_fit_cate_profiles.rds")
  fit <- readRDS(fit_cache)
  if (mode == "production" && nrow(fit$analysis_data) != 9915L) {
    stop(
      "SIPP analysis sample differs from the accepted 9,915-household sample; inspect the data and package versions.",
      call. = FALSE
    )
  }
  replication_run_r(
    file.path(app, "03_profile_diagnostics.R"),
    paste0("--fit-cache=", fit_cache), file.path(output, "diagnostics.log")
  )
  export_args <- c(paste0("--fit-cache=", fit_cache), paste0("--output-dir=", file.path(run, "affine_export")))
  if (mode == "smoke") export_args <- c(export_args, "--allow-validation")
  replication_run_r(file.path(app, "05_export_affine_application.R"), export_args, file.path(output, "export.log"))
  orders <- if (mode == "production") {
    c(99L, 149L, 198L, 248L, 298L, 397L, 496L, 745L)
  } else {
    floor(c(.1, .15, .2, .25, .3, .4, .5, .75) * nrow(fit$analysis_data)^.75)
  }
  localization_args <- c(
    paste0("--fit-cache=", fit_cache), paste0("--output-dir=", file.path(run, "localization")),
    paste0("--orders=", paste(orders, collapse = ","))
  )
  if (mode == "smoke") localization_args <- c(localization_args, "--allow-validation")
  replication_run_r(
    file.path(app, "06_localization_sensitivity.R"), localization_args, file.path(output, "localization.log")
  )
  jsonlite::write_json(
    list(
      status = "completed", mode = mode, sample_n = nrow(fit$analysis_data),
      purpose = if (mode == "smoke") {
        "validation_only"
      } else {
        "fresh current affine SIPP replication; not automatic manuscript promotion"
      },
      paper_commit = replication_manifest(bundle)$paper_commit,
      fit_manifest_sha256 = digest::digest(file = file.path(run, "manifest.json"), algo = "sha256")
    ),
    file.path(output, "replication.json"),
    auto_unbox = TRUE, pretty = TRUE
  )
  invisible(list(mode = mode, output = output, fit_cache = fit_cache))
}
