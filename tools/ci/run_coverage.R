#!/usr/bin/env Rscript

# covr cleans and recompiles its input. Always instrument a retained source copy.
coverage_stop <- function(...) stop(..., call. = FALSE)

coverage_options <- function(args) {
  options <- list(filter = NULL, output_dir = NULL, workspace = NULL,
                  install_path = NULL, cobertura = NULL, html = FALSE,
                  require_native = FALSE)
  i <- 1L
  while (i <= length(args)) {
    arg <- args[[i]]
    if (arg %in% c("--help", "-h")) {
      cat(paste0("Usage: Rscript tools/ci/run_coverage.R [--filter REGEX] ",
                 "[--output-dir DIR] [--workspace DIR] [--install-path DIR] ",
                 "[--cobertura FILE] [--html] [--require-native]\n",
                 "Without --filter, run the full installed test suite. ",
                 "Destinations must not already exist.\n"))
      return(NULL)
    }
    key <- sub("^--", "", arg)
    key <- gsub("-", "_", key, fixed = TRUE)
    if (!startsWith(arg, "--") || !key %in% names(options)) {
      coverage_stop("Unknown coverage argument: ", arg, ". Use --help.")
    }
    if (key %in% c("html", "require_native")) {
      options[[key]] <- TRUE
    } else {
      i <- i + 1L
      if (i > length(args) || !nzchar(args[[i]]) || startsWith(args[[i]], "--")) {
        coverage_stop(arg, " requires a nonempty value.")
      }
      options[[key]] <- args[[i]]
    }
    i <- i + 1L
  }
  options
}

# Resolve existing ancestors too, so symlinks cannot disguise a source path.
coverage_path <- function(path) {
  path <- path.expand(path)
  if (!grepl("^(/|[A-Za-z]:[/\\\\])", path)) path <- file.path(getwd(), path)
  if (file.exists(path) || dir.exists(path)) {
    return(normalizePath(path, winslash = "/", mustWork = TRUE))
  }
  parent <- dirname(path)
  if (identical(parent, path)) coverage_stop("Cannot resolve output path: ", path)
  leaf <- basename(path)
  if (leaf == ".") return(coverage_path(parent))
  if (leaf == "..") return(dirname(coverage_path(parent)))
  file.path(coverage_path(parent), leaf)
}

coverage_inside <- function(path, parent) {
  identical(path, parent) || startsWith(path, paste0(parent, "/"))
}

coverage_destinations <- function(options, root) {
  if (is.null(options$output_dir)) {
    options$output_dir <- file.path(root, "archives", "coverage",
                                   paste0(format(Sys.time(), "%Y%m%dT%H%M%S"), "-", Sys.getpid()))
  }
  options$output_dir <- coverage_path(options$output_dir)
  if (is.null(options$workspace)) options$workspace <- file.path(options$output_dir, "source")
  if (is.null(options$install_path)) options$install_path <- file.path(options$output_dir, "library")
  if (is.null(options$cobertura)) options$cobertura <- file.path(options$output_dir, "cobertura.xml")
  for (key in c("workspace", "install_path", "cobertura")) {
    options[[key]] <- coverage_path(options[[key]])
  }
  for (key in c("output_dir", "workspace", "install_path", "cobertura")) {
    path <- options[[key]]
    link <- Sys.readlink(path)
    if (file.exists(path) || dir.exists(path) || (!is.na(link) && nzchar(link))) {
      coverage_stop("Coverage destination already exists (", key, "): ", path,
                    ". Choose a new destination; existing files are never replaced.")
    }
    allowed_local <- coverage_inside(path, file.path(root, "archives", "coverage"))
    # CI historically publishes this single report at the checkout root.
    ci_report <- key == "cobertura" && identical(path, file.path(root, "cobertura.xml"))
    if ((coverage_inside(path, root) && !allowed_local && !ci_report) ||
        coverage_inside(root, path)) {
      coverage_stop("Unsafe coverage destination (", key, "): ", path,
                    ". Use a new directory under archives/coverage or outside the checkout.")
    }
    ancestor <- dirname(path)
    while (!file.exists(ancestor) && !dir.exists(ancestor)) ancestor <- dirname(ancestor)
    if (!dir.exists(ancestor) || file.access(ancestor, 2L) != 0L) {
      coverage_stop("Coverage destination has no writable parent directory (", key, "): ", path)
    }
  }
  if (coverage_inside(options$workspace, options$install_path) ||
      coverage_inside(options$install_path, options$workspace)) {
    coverage_stop("Coverage workspace and install path must be separate directories.")
  }
  for (path in c(options$output_dir, options$cobertura)) {
    if (coverage_inside(path, options$workspace) || coverage_inside(path, options$install_path)) {
      coverage_stop("Coverage reports must be outside the source workspace and installed library.")
    }
  }
  for (key in c("output_dir", "workspace", "install_path")) {
    if (coverage_inside(options[[key]], options$cobertura)) {
      coverage_stop("Cobertura must be a file, not the ", key,
                    " directory or its ancestor; choose separate paths.")
    }
  }
  reserved <- file.path(options$output_dir, c("coverage.rds", "coverage.csv", "summary.json",
                                             "provenance.json", "installed-execution.rds", "coverage.html"))
  if (any(vapply(reserved, function(path) {
      coverage_inside(path, options$cobertura) || coverage_inside(options$cobertura, path) ||
        coverage_inside(path, options$workspace) || coverage_inside(options$workspace, path) ||
        coverage_inside(path, options$install_path) || coverage_inside(options$install_path, path)
    }, logical(1)))) {
    coverage_stop("Coverage destinations collide with reserved report files; choose separate paths.")
  }
  options
}

coverage_selected_tests <- function(root, filter) {
  files <- list.files(file.path(root, "tests", "testthat"),
                      pattern = "^test.*\\.[Rr]$", full.names = FALSE)
  if (!is.null(filter)) {
    names <- sub("\\.[Rr]$", "", sub("^test[-_]", "", files))
    selected <- tryCatch(grepl(filter, names),
                         error = function(e) coverage_stop("Invalid --filter regex: ", conditionMessage(e)))
    files <- files[selected]
  }
  if (!length(files)) coverage_stop("Coverage selected no test files; check --filter and tests/testthat.")
  files
}

coverage_dependencies <- function(root, html) {
  required <- c(covr = "3.6.5", testthat = "3.2.0", xml2 = "0.0.0", jsonlite = "0.0.0", digest = "0.0.0")
  desc <- read.dcf(file.path(root, "DESCRIPTION"))
  for (field in intersect(c("Depends", "Imports", "LinkingTo"), colnames(desc))) {
    entries <- trimws(strsplit(desc[1L, field], ",", fixed = TRUE)[[1L]])
    for (entry in entries) {
      package <- sub("\\s*\\(.*", "", entry)
      if (identical(package, "R")) next
      if (!package %in% names(required)) required[[package]] <- "0.0.0"
      if (grepl("\\(", entry)) {
        condition <- sub(".*\\(\\s*(.*?)\\s*\\).*", "\\1", entry, perl = TRUE)
        parts <- strsplit(condition, "\\s+")[[1L]]
        if (length(parts) != 2L || !requireNamespace(package, quietly = TRUE) ||
            !do.call(parts[[1L]], list(utils::packageVersion(package), package_version(parts[[2L]])))) {
          coverage_stop("Coverage requires ", package, " (", condition, "); install a compatible dependency first.")
        }
      }
    }
  }
  if (html) required <- c(required, DT = "0.0.0", htmltools = "0.0.0")
  for (tool in c("make", "gcov")) {
    if (!nzchar(Sys.which(tool))) {
      coverage_stop("Coverage requires '", tool, "' on PATH; install the native coverage toolchain first.")
    }
  }
  for (package in names(required)) {
    if (!requireNamespace(package, quietly = TRUE)) {
      coverage_stop("Coverage", if (html && package %in% c("DT", "htmltools")) " --html" else "",
                    " requires package '", package, "'; install it before running coverage.")
    }
  }
  for (package in names(required)) {
    if (utils::packageVersion(package) < package_version(required[[package]])) {
      coverage_stop("Coverage requires ", package, " >= ", required[[package]],
                    "; installed version is ", as.character(utils::packageVersion(package)), ".")
    }
  }
  setNames(lapply(names(required), function(x) as.character(utils::packageVersion(x))), names(required))
}

coverage_source_files <- function(root) {
  ignore <- c("^\\.git($|/)", "^archives($|/)", "^build($|/)",
              "(^|/)\\.ruff_cache($|/)", "(^|/)__pycache__($|/)",
              "\\.(o|so|dll|dylib|gcda|gcno|gcov|pyc)$", "^cobertura\\.xml$")
  if (file.exists(file.path(root, ".Rbuildignore"))) {
    ignore <- c(ignore, readLines(file.path(root, ".Rbuildignore"), warn = FALSE))
  }
  # R CMD build excludes matched directories and their entire contents. Pruning
  # here also avoids traversing old coverage/native workspaces under archives.
  walk <- function(relative = "") {
    entries <- list.files(file.path(root, relative), all.files = TRUE, no.. = TRUE)
    paths <- if (nzchar(relative)) file.path(relative, entries) else entries
    for (pattern in ignore[nzchar(ignore)]) {
      paths <- paths[!grepl(pattern, paths, ignore.case = TRUE, perl = TRUE)]
    }
    directories <- dir.exists(file.path(root, paths))
    nested <- unlist(lapply(paths[directories], walk), use.names = FALSE)
    c(paths[!directories], nested)
  }
  sort(walk())
}

coverage_snapshot <- function(root, workspace, files) {
  dir.create(workspace, recursive = TRUE, showWarnings = FALSE)
  for (file in files) {
    destination <- file.path(workspace, file)
    dir.create(dirname(destination), recursive = TRUE, showWarnings = FALSE)
    if (!file.copy(file.path(root, file), destination, copy.mode = TRUE)) {
      coverage_stop("Cannot snapshot source file: ", file)
    }
  }
  setNames(vapply(files, function(file) digest::digest(file = file.path(workspace, file), algo = "sha256"), ""), files)
}

coverage_command <- function(command, args) {
  if (!nzchar(command)) return(list(command = command, args = args, status = 127L, output = "Tool not found."))
  output <- suppressWarnings(system2(command, args, stdout = TRUE, stderr = TRUE))
  list(command = command, args = args, status = if (is.null(attr(output, "status"))) 0L else attr(output, "status"),
       output = unname(output))
}

coverage_main <- function(args = commandArgs(trailingOnly = TRUE)) {
  options <- coverage_options(args)
  if (is.null(options)) return(invisible(NULL))
  script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (length(script_arg) != 1L) coverage_stop("Run coverage with Rscript tools/ci/run_coverage.R.")
  root <- normalizePath(file.path(dirname(sub("^--file=", "", script_arg)), "..", ".."),
                        winslash = "/", mustWork = TRUE)
  if (!file.exists(file.path(root, "DESCRIPTION")) ||
      read.dcf(file.path(root, "DESCRIPTION"))[1L, "Package"] != "RdistNN") {
    coverage_stop("Coverage runner must belong to an RdistNN package checkout.")
  }
  selected <- coverage_selected_tests(root, options$filter)
  options <- coverage_destinations(options, root)
  files <- coverage_source_files(root)
  if (!all(c("DESCRIPTION", "NAMESPACE", "tests/testthat.R") %in% files)) {
    coverage_stop("Coverage source snapshot is missing DESCRIPTION, NAMESPACE, or tests/testthat.R.")
  }
  if (!all(file.path("tests", "testthat", selected) %in% files)) {
    coverage_stop("Selected coverage tests are excluded from the source snapshot; check .Rbuildignore.")
  }
  versions <- coverage_dependencies(root, options$html)
  r_binary <- file.path(R.home("bin"), "R")
  compiler <- coverage_command(r_binary, c("CMD", "config", "CXX17"))
  compiler_name <- strsplit(trimws(paste(compiler$output, collapse = " ")), "\\s+")[[1L]][[1L]]
  if (compiler$status != 0L || !nzchar(Sys.which(compiler_name))) {
    coverage_stop("Coverage cannot resolve R's CXX17 compiler: ", paste(compiler$output, collapse = " "))
  }
  # Preflight is complete. No prior operation has created or cleaned outputs.
  dir.create(options$output_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(dirname(options$install_path), recursive = TRUE, showWarnings = FALSE)
  dir.create(dirname(options$cobertura), recursive = TRUE, showWarnings = FALSE)
  source_hashes <- coverage_snapshot(root, options$workspace, files)
  provenance_file <- file.path(options$output_dir, "provenance.json")
  provenance <- list(status = "running", started_utc = format(Sys.time(), tz = "UTC", usetz = TRUE),
                     mode = if (is.null(options$filter)) "full" else "focused",
                     filter = options$filter, selected_tests = selected, source_root = root,
                     workspace = options$workspace, install_path = options$install_path,
                     reports = options$output_dir, R = unclass(R.version), tools = versions,
                     compiler = compiler, compiler_version = coverage_command(compiler_name, "--version"),
                     gcov = coverage_command(Sys.which("gcov"), "--version"),
                     make = coverage_command(Sys.which("make"), "--version"),
                     lib_paths = .libPaths(), makeflags = "-j2",
                     runtime_threads = as.list(Sys.getenv(c("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS",
                                                            "MKL_NUM_THREADS"))),
                     covr_flags = getOption("covr.flags"), source_sha256 = as.list(source_hashes),
                     git_head = coverage_command(Sys.which("git"), c("-C", shQuote(root), "rev-parse", "HEAD")),
                     git_status = coverage_command(Sys.which("git"), c("-C", shQuote(root), "status", "--porcelain")))
  write_provenance <- function() {
    jsonlite::write_json(provenance, provenance_file, pretty = TRUE, auto_unbox = TRUE, null = "null")
  }
  write_provenance()
  completed <- FALSE
  on.exit({
    if (!completed) {
      provenance$status <- "failed"
      provenance$finished_utc <- format(Sys.time(), tz = "UTC", usetz = TRUE)
      write_provenance()
    }
  }, add = TRUE)
  old_makeflags <- Sys.getenv("MAKEFLAGS", unset = NA_character_)
  Sys.setenv(MAKEFLAGS = "-j2")
  on.exit(if (is.na(old_makeflags)) Sys.unsetenv("MAKEFLAGS") else Sys.setenv(MAKEFLAGS = old_makeflags), add = TRUE)
  code <- character()
  execution_file <- file.path(options$output_dir, "installed-execution.rds")
  # This code runs in covr's instrumented installed library, never pkgload.
  if (!is.null(options$filter)) {
    code <- c(code, sprintf(paste0("results <- testthat::test_package('RdistNN', filter = %s, ",
                                  "reporter = 'summary', stop_on_failure = TRUE)"), deparse(options$filter)),
              "if (!length(results)) stop('Coverage executed no tests.')")
  }
  code <- c(code, sprintf(paste0("saveRDS(list(package_path = find.package('RdistNN'), ",
                                "library_paths = .libPaths(), ",
                                "session_info = capture.output(sessionInfo()), ",
                                "dlls = lapply(getLoadedDLLs(), function(x) x[['path']])), %s)"),
                          deparse(execution_file)))
  cov <- covr::package_coverage(path = options$workspace,
                               type = if (is.null(options$filter)) "tests" else "none",
                               code = code, quiet = FALSE, clean = FALSE,
                               install_path = options$install_path)
  execution <- readRDS(execution_file)
  expected_package <- normalizePath(file.path(options$install_path, "RdistNN"), winslash = "/", mustWork = TRUE)
  if (!identical(normalizePath(execution$package_path, winslash = "/"), expected_package)) {
    coverage_stop("Coverage tests loaded a package outside the instrumented installed library.")
  }
  table <- as.data.frame(cov)
  if (!nrow(table) || !any(grepl("^R/", table$filename) & table$value > 0)) {
    coverage_stop("Coverage recorded no positive R source hits; inspect the retained workspace and test logs.")
  }
  native <- grepl("^src/.*\\.(cpp|cc|cxx)$", table$filename) & !grepl("(^|/)RcppExports\\.cpp$", table$filename)
  native_hit <- native & table$value > 0
  summary <- list(mode = provenance$mode, filter = options$filter, selected_tests = selected,
                  percent = covr::percent_coverage(cov),
                  R_hit_expressions = sum(grepl("^R/", table$filename) & table$value > 0),
                  native_hit_expressions = sum(native_hit), native_hit_files = unique(table$filename[native_hit]))
  dlls <- unlist(execution$dlls, use.names = TRUE)
  dlls <- dlls[file.exists(dlls) & startsWith(dlls, paste0(expected_package, "/"))]
  provenance$execution <- execution
  provenance$installed_dll_sha256 <- as.list(setNames(vapply(dlls, function(x) {
    digest::digest(file = x, algo = "sha256")
  }, ""), dlls))
  provenance$summary <- summary
  if ((is.null(options$filter) || options$require_native) && !any(native_hit)) {
    coverage_stop("Coverage recorded no positive owned C++ source hits; native execution is required for this run.")
  }
  saveRDS(cov, file.path(options$output_dir, "coverage.rds"))
  utils::write.csv(table, file.path(options$output_dir, "coverage.csv"), row.names = FALSE)
  jsonlite::write_json(summary, file.path(options$output_dir, "summary.json"),
                       pretty = TRUE, auto_unbox = TRUE, null = "null")
  covr::to_cobertura(cov, filename = options$cobertura)
  if (options$html) covr::report(cov, file = file.path(options$output_dir, "coverage.html"), browse = FALSE)
  provenance$status <- "passed"
  provenance$finished_utc <- format(Sys.time(), tz = "UTC", usetz = TRUE)
  write_provenance()
  completed <- TRUE
  print(cov)
  cat("Coverage reports: ", options$output_dir, "\n", sep = "")
  invisible(cov)
}

if (sys.nframe() == 0L) coverage_main()
