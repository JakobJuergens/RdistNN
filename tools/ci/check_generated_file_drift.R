#!/usr/bin/env Rscript

if (!requireNamespace("roxygen2", quietly = TRUE)) {
  stop("Package 'roxygen2' is required for this check.", call. = FALSE)
}
if (!requireNamespace("Rcpp", quietly = TRUE)) {
  stop("Package 'Rcpp' is required for this check.", call. = FALSE)
}

copy_if_exists <- function(path, dest_root) {
  if (!file.exists(path)) {
    return(invisible(FALSE))
  }
  if (dir.exists(path)) {
    ok <- file.copy(path, dest_root, recursive = TRUE, overwrite = TRUE)
  } else {
    dest <- file.path(dest_root, path)
    dir.create(dirname(dest), recursive = TRUE, showWarnings = FALSE)
    ok <- file.copy(path, dest, overwrite = TRUE)
  }
  if (!isTRUE(ok)) {
    stop(sprintf("Failed to copy '%s' into temporary check directory.", path), call. = FALSE)
  }
  invisible(TRUE)
}

read_raw_file <- function(path) {
  if (!file.exists(path)) {
    return(NULL)
  }
  sz <- file.info(path)$size
  con <- file(path, open = "rb")
  on.exit(close(con), add = TRUE)
  if (is.na(sz) || sz <= 0L) {
    raw(0L)
  } else {
    readBin(con, what = "raw", n = sz)
  }
}

repo_root <- normalizePath(".", winslash = "/", mustWork = TRUE)
tmp_root <- tempfile("generated-drift-")
dir.create(tmp_root, recursive = TRUE, showWarnings = FALSE)
on.exit(unlink(tmp_root, recursive = TRUE, force = TRUE), add = TRUE)

# Minimal package skeleton needed to regenerate namespace + Rcpp exports.
for (p in c("DESCRIPTION", "NAMESPACE", "R", "src", "man")) {
  copy_if_exists(p, tmp_root)
}

tmp_pkg <- tmp_root
has_src <- dir.exists(file.path(tmp_pkg, "src")) &&
  length(list.files(file.path(tmp_pkg, "src"), recursive = TRUE, all.files = TRUE)) > 0L

roxygen2::roxygenise(
  package.dir = tmp_pkg,
  roclets = c("namespace"),
  load = "source"
)

if (has_src) {
  Rcpp::compileAttributes(pkgdir = tmp_pkg, verbose = FALSE)
}

targets <- c("NAMESPACE", "R/RcppExports.R", "src/RcppExports.cpp")
drift <- character(0)
missing_after <- character(0)

for (target in targets) {
  orig <- read_raw_file(file.path(repo_root, target))
  regen <- read_raw_file(file.path(tmp_pkg, target))

  if (is.null(regen)) {
    missing_after <- c(missing_after, target)
    next
  }
  if (is.null(orig) || !identical(orig, regen)) {
    drift <- c(drift, target)
  }
}

if (length(missing_after)) {
  cat("Generated-file check could not find expected regenerated targets:\n")
  for (p in unique(missing_after)) cat(sprintf("- %s\n", p))
  quit(status = 1L)
}

if (length(drift)) {
  cat("Generated-file drift detected. Regenerate and commit these files:\n")
  for (p in unique(drift)) cat(sprintf("- %s\n", p))
  cat(
    "\nSuggested local commands:\n",
    "  R -q -e \"roxygen2::roxygenise(roclets='namespace')\"\n",
    "  R -q -e \"Rcpp::compileAttributes()\"\n",
    sep = ""
  )
  quit(status = 1L)
}

cat("No generated-file drift detected for NAMESPACE and RcppExports files.\n")
