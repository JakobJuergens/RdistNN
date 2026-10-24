args <- commandArgs(trailingOnly = FALSE)
script_arg <- grep("^--file=", args, value = TRUE)
if (length(script_arg) != 1L) {
  stop("Run this check with Rscript.", call. = FALSE)
}

script_path <- sub("^--file=", "", script_arg)
repo_root <- normalizePath(
  file.path(dirname(script_path), "..", ".."),
  winslash = "/",
  mustWork = TRUE
)
setwd(repo_root)

errors <- character()
add_error <- function(message) {
  errors <<- c(errors, message)
}

# DESCRIPTION owns the current version; history and installation pins do not.
description <- if (file.exists("DESCRIPTION")) readLines("DESCRIPTION", warn = FALSE) else character()
declarations <- trimws(sub("^Version:", "", grep("^Version:", description, value = TRUE)))
if (length(declarations) != 1L ||
    !grepl("^[0-9]+[.][0-9]+[.][0-9]+([.][0-9]+)?$", declarations)) {
  add_error("DESCRIPTION: declare exactly one three- or four-component package Version")
} else {
  version <- declarations[[1L]]
  patterns <- c(
    "NEWS.md" = "^# RdistNN ([0-9]+([.][0-9]+){2,3})( \\(([0-9]{4}-[0-9]{2}-[0-9]{2})\\))?[[:space:]]*$",
    "README.Rmd" = "^Version \\*\\*([0-9]+([.][0-9]+){2,3})\\*\\*",
    "README.md" = "^Version \\*\\*([0-9]+([.][0-9]+){2,3})\\*\\*"
  )
  for (path in names(patterns)) {
    if (!file.exists(path)) {
      add_error(sprintf("%s: missing current package version owner", path))
      next
    }
    lines <- readLines(path, warn = FALSE)
    if (identical(path, "NEWS.md")) {
      lines <- head(grep("^# ", lines, value = TRUE), 1L)
    }
    matches <- regmatches(lines, regexec(patterns[[path]], lines))
    matches <- matches[lengths(matches) > 0L]
    if (length(matches) != 1L) {
      add_error(sprintf("%s: declare one current package version (%s)", path, version))
    } else if (matches[[1L]][[2L]] != version) {
      add_error(sprintf("%s: current package version %s differs from DESCRIPTION %s",
                        path, matches[[1L]][[2L]], version))
    }
    if (identical(path, "NEWS.md") && length(matches) == 1L) {
      date <- matches[[1L]][[5L]]
      if (nzchar(date)) {
        parsed_date <- as.Date(date, format = "%Y-%m-%d")
        if (is.na(parsed_date) || format(parsed_date, "%Y-%m-%d") != date) {
          add_error("NEWS.md: date suffix must be a valid ISO date (YYYY-MM-DD)")
        }
      }
    }
  }
}

# Direct namespace definitions have one owner; function-local closures may use
# the same names because they do not replace package-level bindings.
function_owners <- list()
for (path in list.files("R", "[.]R$", full.names = TRUE)) {
  for (expr in parse(path, keep.source = FALSE)) {
    if (is.call(expr) && is.symbol(expr[[1L]]) &&
        as.character(expr[[1L]]) %in% c("<-", "=")) {
      # Inspect the RHS only for direct assignments, not bodies or local scopes.
      if (is.symbol(expr[[2L]]) && is.call(expr[[3L]]) &&
          identical(expr[[3L]][[1L]], as.name("function"))) {
        name <- as.character(expr[[2L]])
        function_owners[[name]] <- c(function_owners[[name]], path)
      }
    }
  }
}
for (name in names(function_owners)) {
  owners <- function_owners[[name]]
  if (length(owners) > 1L) {
    add_error(sprintf("duplicate function definition %s: %s", name,
                      paste(owners, collapse = ", ")))
  }
}

tracked <- system2("git", c("ls-files", "--cached", "--others", "--exclude-standard"), stdout = TRUE, stderr = TRUE)
git_status <- attr(tracked, "status")
if (!is.null(git_status) && git_status != 0L) {
  stop("Could not enumerate repository inputs with git ls-files.", call. = FALSE)
}

forbidden_tracked <- c(
  ".claude/settings.local.json",
  "Jackknife-GenUStat_JJuergens.pdf",
  "Local_Z_Est_JJuergens.pdf"
)
for (path in forbidden_tracked) {
  if (path %in% tracked) {
    add_error(sprintf("local-only file is tracked: %s", path))
  }
}

for (pattern in c(
  "(^|/)__pycache__/",
  "[.]py[co]$",
  "^build/",
  "^[^/]+[.]log$",
  "^GH_logs/",
  "^archives/",
  "[.]Rcheck/",
  "[.](o|so|dll|dylib)$",
  "^README[.]html$"
)) {
  matches <- grep(pattern, tracked, value = TRUE)
  if (length(matches)) {
    add_error(sprintf(
      "generated repository inputs match %s: %s",
      pattern,
      paste(matches, collapse = ", ")
    ))
  }
}

ignore_lines <- readLines(".gitignore", warn = FALSE)
required_ignores <- c(
  ".claude/",
  "/Jackknife-GenUStat_JJuergens.pdf",
  "/Local_Z_Est_JJuergens.pdf",
  "GH_logs/",
  "archives/",
  "README.html",
  "/build/"
)
for (entry in required_ignores) {
  if (!entry %in% ignore_lines) {
    add_error(sprintf("missing .gitignore entry: %s", entry))
  }
}

workflow_files <- list.files(
  ".github/workflows",
  pattern = "[.](yml|yaml)$",
  full.names = TRUE
)
workflow_text <- paste(
  unlist(lapply(workflow_files, readLines, warn = FALSE), use.names = FALSE),
  collapse = "\n"
)
for (pattern in c(
  "commit_logs:",
  "EndBug/add-and-commit",
  "git[[:space:]]+commit",
  "git[[:space:]]+push"
)) {
  if (grepl(pattern, workflow_text)) {
    add_error(sprintf("obsolete repository-writing workflow pattern: %s", pattern))
  }
}

for (workflow in workflow_files) {
  parsed <- yaml::read_yaml(workflow)
  for (job_name in names(parsed$jobs)) {
    steps <- parsed$jobs[[job_name]]$steps
    if (is.null(steps)) {
      next
    }
    for (i in seq_along(steps)) {
      uses <- steps[[i]]$uses
      if (is.null(uses) || !grepl("^actions/upload-artifact@", uses)) {
        next
      }
      retention <- steps[[i]]$with[["retention-days"]]
      if (is.null(retention)) {
        add_error(sprintf(
          "%s job %s upload step %d has no retention-days",
          workflow,
          job_name,
          i
        ))
      }
    }
  }
}

if (length(errors)) {
  cat("Repository hygiene check failed:\n")
  cat(sprintf("- %s\n", errors))
  quit(status = 1L)
}

cat(sprintf(
  "Repository hygiene check passed: %d repository inputs, %d workflows.\n",
  length(tracked),
  length(workflow_files)
))
