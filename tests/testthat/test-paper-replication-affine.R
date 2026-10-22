.affine_bundle <- function() {
  path <- system.file("paper-replication", "affine", package = "RdistNN")
  expect_true(nzchar(path), info = "Current affine replication bundle must be installed")
  path
}

.affine_helpers <- function() {
  path <- .affine_bundle()
  if (!nzchar(path)) {
    return(NULL)
  }
  env <- new.env(parent = globalenv())
  sys.source(file.path(path, "helpers.R"), envir = env)
  list(path = path, env = env)
}

test_that("renderer adapters preserve caller bindings when nested scripts reject arguments", {
  bundle <- .affine_helpers()
  if (is.null(bundle)) return(invisible(NULL))
  before_names <- ls(globalenv(), all.names = TRUE)
  before <- as.list(globalenv(), all.names = TRUE)
  on.exit({
    added <- setdiff(ls(globalenv(), all.names = TRUE), before_names)
    if (length(added)) rm(list = added, envir = globalenv())
    list2env(before, envir = globalenv())
  }, add = TRUE)
  sentinels <- c("plot_sipp_401k_affine_summaries", "plot_sipp_401k_localization")
  for (name in sentinels) assign(name, function(...) "caller helper", envir = globalenv())
  expected <- mget(sentinels, envir = globalenv(), inherits = FALSE)
  expected_names <- ls(globalenv(), all.names = TRUE)
  script <- file.path(bundle$path, "source", "sipp", "07_render_affine_manuscript.R")
  expect_error(bundle$env$replication_run_renderer(script, "--unknown-audit-probe", "unused"),
    "Unknown argument")
  expect_identical(mget(sentinels, envir = globalenv(), inherits = FALSE), expected)
  expect_setequal(ls(globalenv(), all.names = TRUE), expected_names)
})

test_that("renderer adapters keep nested dependencies in their requested environments", {
  bundle <- .affine_helpers()
  if (is.null(bundle)) return(invisible(NULL))
  directory <- tempfile("renderer scope (test) ")
  dir.create(directory)
  on.exit(unlink(directory, recursive = TRUE), add = TRUE)
  script <- file.path(directory, "script.R")
  first <- file.path(directory, "first.R")
  second <- file.path(directory, "second.R")
  output <- file.path(directory, "result.rds")
  literal <- function(path) encodeString(path, quote = '"')
  writeLines(".rdistnn_renderer_scope_helper <- function() 'private dependency'", second)
  writeLines(paste0("source(", literal(second), ")"), first)
  writeLines(c(
    paste0("source(", literal(first), ")"),
    "nested_local <- identical(environment(.rdistnn_renderer_scope_helper), environment())",
    "explicit <- new.env(parent = environment())",
    paste0("source(", literal(second), ", local = explicit)"),
    "explicit_local <- identical(environment(explicit$.rdistnn_renderer_scope_helper), explicit)",
    paste0("source(", literal(second), ", local = TRUE)"),
    "true_local <- identical(environment(.rdistnn_renderer_scope_helper), environment())",
    paste0("saveRDS(list(value = .rdistnn_renderer_scope_helper(), nested = nested_local, ",
      "explicit = explicit_local, local_true = true_local), ", literal(output), ")")
  ), script)
  before_names <- ls(globalenv(), all.names = TRUE)
  before <- as.list(globalenv(), all.names = TRUE)
  on.exit({
    added <- setdiff(ls(globalenv(), all.names = TRUE), before_names)
    if (length(added)) rm(list = added, envir = globalenv())
    list2env(before, envir = globalenv())
  }, add = TRUE)
  expect_true(bundle$env$replication_run_renderer(script, character(), "unused"))
  expect_identical(readRDS(output), list(value = "private dependency", nested = TRUE,
    explicit = TRUE, local_true = TRUE))
  expect_setequal(ls(globalenv(), all.names = TRUE), before_names)
})

test_that("affine replication verifies scientific source and rejects corruption", {
  skip_if_not_installed("digest")
  skip_if_not_installed("jsonlite")
  bundle <- .affine_helpers()
  if (is.null(bundle)) {
    return(invisible(NULL))
  }
  expect_true(bundle$env$replication_verify(bundle$path))
  copy <- tempfile("affine-corruption-")
  dir.create(copy)
  on.exit(unlink(copy, recursive = TRUE), add = TRUE)
  file.copy(list.files(bundle$path, full.names = TRUE), copy, recursive = TRUE)
  file <- file.path(copy, "source", "sim", "affine_models.R")
  cat("\n# changed scientific source\n", file = file, append = TRUE)
  expect_error(bundle$env$replication_verify(copy), "hash mismatch")
})

test_that("production configuration covers the accepted affine campaign", {
  skip_if_not_installed("digest")
  skip_if_not_installed("jsonlite")
  bundle <- .affine_helpers()
  if (is.null(bundle)) {
    return(invisible(NULL))
  }
  configs <- bundle$env$replication_configs(bundle$path, "production")
  expect_identical(names(configs), c("dnn", "kernel_comparison", "M3", "settings_comparison"))
  expect_equal(vapply(configs, function(x) nrow(x$cells), integer(1)), c(dnn = 14L, kernel_comparison = 9L, M3 = 1L, settings_comparison = 4L))
  expect_true(all(vapply(configs, function(x) all(x$cells$replications == 10000L), logical(1))))
  expect_equal(configs$M3$cells$n, 531441L)
  expect_equal(configs$M3$cells$M, 3L)
  rows <- sum(vapply(configs, function(x) sum(x$cells$replications) * 6 * length(x$localization_methods) + sum(x$cells$control_replications) * 4 * length(x$localization_methods), numeric(1)))
  expect_equal(rows, 2516400)
})

test_that("affine simulation smoke run executes both applications and all prediction paths", {
  for (pkg in c("digest", "jsonlite", "future", "future.apply")) skip_if_not_installed(pkg)
  bundle <- .affine_helpers()
  if (is.null(bundle)) {
    return(invisible(NULL))
  }
  output <- tempfile("affine-smoke-")
  on.exit(unlink(output, recursive = TRUE), add = TRUE)
  result <- bundle$env$replication_simulations(bundle$path, output, mode = "smoke", profiles = "dnn")
  raw <- readRDS(file.path(output, "runs", "dnn", "raw.rds"))
  expect_equal(nrow(raw), 12L)
  expect_setequal(raw$application, c("cate", "slope"))
  expect_setequal(raw$training_method, c("full", "pooled_five_fold", "oracle"))
  expect_true(all(raw$status == "ok"))
  expect_true(all(is.finite(raw$estimate)))
  expect_true(all(raw$point_target == 1 / 8))
  expect_identical(result$mode, "smoke")
  expect_error(bundle$env$replication_simulations(bundle$path, output, mode = "smoke", profiles = "dnn"), "already exists")
})

test_that("replication arguments reject unsupported or ambiguous requests early", {
  bundle <- .affine_helpers()
  if (is.null(bundle)) {
    return(invisible(NULL))
  }
  expect_error(bundle$env$replication_args(c("--mode=old")), "mode")
  expect_error(bundle$env$replication_args(c("--unknown=yes")), "Unknown")
  expect_error(bundle$env$replication_args(c("--mode=smoke", "--mode=production")), "Repeated")
})

test_that("accepted evidence regenerates all current manuscript tables", {
  for (pkg in c("digest", "jsonlite", "ggplot2", "patchwork")) skip_if_not_installed(pkg)
  bundle <- .affine_helpers()
  if (is.null(bundle)) {
    return(invisible(NULL))
  }
  output <- tempfile("affine [evidence] (test) ")
  on.exit(unlink(output, recursive = TRUE), add = TRUE)
  expect_true(bundle$env$replication_render(bundle$path, output))
  tables <- list.files(file.path(output, "Tables"), pattern = "\\.tex$", recursive = TRUE)
  expect_length(tables, 21L)
  for (file in tables) expect_identical(readLines(file.path(output, "Tables", file)), readLines(file.path(output, "reference", "Tables", file)))
  figures <- list.files(file.path(output, "Artifacts"), pattern = "\\.pdf$", recursive = TRUE)
  expect_length(figures, 6L)
  expect_error(bundle$env$replication_render(bundle$path, output), "already exists")
})

test_that("SIPP smoke replication retains three distinct prediction mechanisms", {
  for (pkg in c("digest", "jsonlite", "hdm", "grf")) skip_if_not_installed(pkg)
  bundle <- .affine_helpers()
  if (is.null(bundle)) {
    return(invisible(NULL))
  }
  output <- tempfile("sipp smoke ")
  on.exit(unlink(output, recursive = TRUE), add = TRUE)
  bundle$env$replication_sipp(bundle$path, output, mode = "smoke")
  fit <- readRDS(file.path(output, "runs", "paper-current", "02_fit_cate_profiles.rds"))
  expect_equal(nrow(fit$analysis_data), 300L)
  expect_equal(nrow(fit$estimates), 21L)
  expect_setequal(fit$estimates$method, c("dnn_fullsample", "dnn_oob", "dnn_crossfit"))
  expect_true(all(is.finite(fit$estimates$tau)))
  expect_true(isTRUE(fit$run_info$smoke))
  expect_gt(max(abs(fit$nuisances$fullsample$e_hat - fit$nuisances$oob$e_hat)), 0)
  localization <- read.csv(file.path(output, "runs", "paper-current", "localization", "localization_all.csv"))
  expect_equal(nrow(localization), 168L)
  expect_true(file.exists(file.path(output, "runs", "paper-current", "affine_export", "information.csv")))
})

test_that("audit adapter preserves settings while restoring frozen serialization bytes", {
  for (pkg in c("digest", "jsonlite")) skip_if_not_installed(pkg)
  bundle <- .affine_helpers()
  if (is.null(bundle)) {
    return(invisible(NULL))
  }
  root <- tempfile("affine-config-")
  dir.create(root)
  on.exit(unlink(root, recursive = TRUE), add = TRUE)
  frozen <- file.path(bundle$path, "configs", "dnn", "config.rds")
  path <- file.path(root, "config.rds")
  cfg <- readRDS(frozen)
  saveRDS(cfg, path, version = 2L)
  expect_identical(readRDS(path), cfg)
  expect_false(identical(readBin(path, "raw", n = file.info(path)$size), readBin(frozen, "raw", n = file.info(frozen)$size)))
  bundle$env$replication_restore_config(frozen, path)
  .expect_identical_fast(readBin(path, "raw", n = file.info(path)$size), readBin(frozen, "raw", n = file.info(frozen)$size))
  cfg$master_seed <- 1L
  saveRDS(cfg, path)
  expect_error(bundle$env$replication_restore_config(frozen, path), "settings differ")
})

test_that("an interrupted analysis directory cannot masquerade as completed evidence", {
  for (pkg in c("digest", "jsonlite")) skip_if_not_installed(pkg)
  bundle <- .affine_helpers()
  if (is.null(bundle)) {
    return(invisible(NULL))
  }
  root <- tempfile("affine-analysis-")
  dir.create(root)
  on.exit(unlink(root, recursive = TRUE), add = TRUE)
  destination <- file.path(root, "analysis")
  dir.create(destination)
  writeLines("interrupted", file.path(destination, "partial.csv"))
  called <- FALSE
  render <- function(path) {
    called <<- TRUE
    dir.create(path)
    writeLines("complete", file.path(path, "metric.csv"))
    jsonlite::write_json(
      list(status = "analysis_of_audited_saved_production", data_hashes = data.frame(
        file = "metric.csv", sha256 = digest::digest(file = file.path(path, "metric.csv"), algo = "sha256")
      )),
      file.path(path, "analysis_manifest.json"),
      auto_unbox = TRUE
    )
  }
  bundle$env$replication_stage(destination, "analysis", render)
  expect_true(called)
  expect_true(file.exists(file.path(destination, "analysis_manifest.json")))
  expect_true(any(grepl("analysis.interrupted", list.files(root), fixed = TRUE)))
  called <- FALSE
  bundle$env$replication_stage(destination, "analysis", render)
  expect_false(called)
  unlink(file.path(destination, "metric.csv"))
  bundle$env$replication_stage(destination, "analysis", render)
  expect_true(called)
})

test_that("saved analysis must reference the audit it was computed from", {
  for (pkg in c("digest", "jsonlite")) skip_if_not_installed(pkg)
  bundle <- .affine_helpers()
  if (is.null(bundle)) {
    return(invisible(NULL))
  }
  root <- tempfile("affine-audit-binding-")
  dir.create(root)
  on.exit(unlink(root, recursive = TRUE), add = TRUE)
  audit <- file.path(root, "audit.json")
  writeLines("current audit", audit)
  writeLines("complete metrics", file.path(root, "metric.csv"))
  jsonlite::write_json(
    list(
      status = "analysis_of_audited_saved_production", source_audit_sha256 = "stale audit",
      data_hashes = data.frame(file = "metric.csv", sha256 = digest::digest(file = file.path(root, "metric.csv"), algo = "sha256"))
    ),
    file.path(root, "analysis_manifest.json"),
    auto_unbox = TRUE
  )
  expect_false(bundle$env$replication_stage_valid(root, "analysis", audit))
})
