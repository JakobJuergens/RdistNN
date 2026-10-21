#!/usr/bin/env Rscript
# Run installed tests with only Rcpp, testthat and their strong dependencies.
main <- function(args = commandArgs(TRUE)) {
  if (length(args) != 1L || !file.exists(args[1L])) {
    stop("Usage: Rscript tools/ci/check_minimal_dependencies.R PACKAGE.tar.gz", call. = FALSE)
  }
  archive <- normalizePath(args[1L], mustWork = TRUE)
  temporary <- tempfile("minimal-dependencies-")
  dir.create(temporary)
  on.exit(unlink(temporary, recursive = TRUE, force = TRUE), add = TRUE)
  library <- file.path(temporary, "library")
  dir.create(library)
  installed <- utils::installed.packages()
  roots <- c("Rcpp", "testthat")
  dependencies <- tools::package_dependencies(roots, db = installed,
    which = c("Depends", "Imports", "LinkingTo"), recursive = TRUE)
  packages <- setdiff(unique(c(roots, unlist(dependencies, use.names = FALSE))), "R")
  for (package in packages) {
    path <- find.package(package)
    if (dirname(path) == .Library) next
    target <- file.path(library, package)
    if (!isTRUE(suppressWarnings(file.symlink(path, target))) &&
        !isTRUE(file.copy(path, library, recursive = TRUE))) {
      stop("Cannot isolate dependency: ", package, call. = FALSE)
    }
  }
  run <- function(command, arguments, env = character()) {
    status <- system2(command, arguments, env = env)
    if (status != 0L) stop("Minimal-dependency check failed (exit ", status, ").", call. = FALSE)
  }
  run(file.path(R.home("bin"), "R"),
    c("CMD", "INSTALL", "--install-tests", "--no-multiarch",
      paste0("--library=", shQuote(library)), shQuote(archive)),
    env = paste0(c("R_LIBS=", "R_LIBS_USER=", "R_LIBS_SITE="), shQuote(library)))
  runner <- file.path(temporary, "tests.R")
  writeLines(c(
    ".libPaths(c(commandArgs(TRUE)[1L], .Library), include.site = FALSE)",
    "optional <- c('ggplot2', 'patchwork', 'grf', 'hdm', 'future.apply', 'tibble')",
    "if (any(vapply(optional, requireNamespace, logical(1), quietly = TRUE)))",
    "  stop('Optional feature dependencies must be absent from the isolated library.')",
    "testthat::test_package('RdistNN', reporter = 'summary', stop_on_failure = TRUE)"
  ), runner)
  run(file.path(R.home("bin"), "Rscript"),
    c("--vanilla", shQuote(runner), shQuote(library)))
  cat("Minimal-dependency installed tests passed.\n")
}

main()
