#!/usr/bin/env Rscript

if (!requireNamespace("yaml", quietly = TRUE)) {
  stop(
    "Package 'yaml' is required for this check. Install it first.",
    call. = FALSE
  )
}

flatten_contents <- function(x) {
  if (is.null(x)) {
    return(character(0))
  }
  if (is.character(x)) {
    return(x)
  }
  if (is.list(x)) {
    return(unlist(lapply(x, flatten_contents), use.names = FALSE))
  }
  character(0)
}

match_selector <- function(entry, topics) {
  parsed <- regexec("^([a-z_]+)\\((\"[^\"]*\"|'[^']*')\\)$", entry, perl = TRUE)
  bits <- regmatches(entry, parsed)[[1L]]
  if (!length(bits)) {
    return(list(recognized = FALSE, matches = character(0)))
  }

  fn <- bits[2L]
  arg <- bits[3L]
  arg <- substring(arg, 2L, nchar(arg) - 1L)

  matched <- switch(
    fn,
    starts_with = topics[startsWith(topics, arg)],
    ends_with = topics[endsWith(topics, arg)],
    contains = topics[grepl(arg, topics, fixed = TRUE)],
    matches = {
      ok <- tryCatch(grepl(arg, topics, perl = TRUE), error = function(e) NULL)
      if (is.null(ok)) character(0) else topics[ok]
    },
    character(0)
  )

  list(recognized = fn %in% c("starts_with", "ends_with", "contains", "matches"), matches = matched)
}

is_internal_topic <- function(rd_path) {
  lines <- readLines(rd_path, warn = FALSE, encoding = "UTF-8")
  any(grepl("\\\\keyword\\{internal\\}", lines, perl = TRUE))
}

cfg_path <- "_pkgdown.yml"
if (!file.exists(cfg_path)) {
  stop("_pkgdown.yml not found.", call. = FALSE)
}

cfg <- yaml::read_yaml(cfg_path)
if (is.null(cfg$reference) || !length(cfg$reference)) {
  stop("_pkgdown.yml has no `reference` index blocks.", call. = FALSE)
}

rd_files <- list.files("man", pattern = "\\.Rd$", full.names = TRUE)
if (!length(rd_files)) {
  stop("No .Rd files found under man/.", call. = FALSE)
}

topics_all <- sub("\\.Rd$", "", basename(rd_files))
internal <- topics_all[vapply(rd_files, is_internal_topic, logical(1L))]
topics_public <- sort(setdiff(topics_all, internal))

entries <- unique(unlist(lapply(cfg$reference, function(block) {
  flatten_contents(block$contents)
}), use.names = FALSE))
entries <- trimws(entries)
entries <- entries[nzchar(entries)]

matched_topics <- character(0)
unknown_entries <- character(0)

for (entry in entries) {
  selector <- match_selector(entry, topics_all)
  if (selector$recognized) {
    matched_topics <- c(matched_topics, selector$matches)
    next
  }
  if (entry %in% topics_all) {
    matched_topics <- c(matched_topics, entry)
    next
  }
  unknown_entries <- c(unknown_entries, entry)
}

matched_topics <- unique(matched_topics)
missing_topics <- sort(setdiff(topics_public, matched_topics))
unknown_entries <- unique(unknown_entries)

if (length(unknown_entries)) {
  cat("Unknown entries in _pkgdown.yml reference contents:\n")
  for (entry in sort(unknown_entries)) {
    cat(sprintf("- %s\n", entry))
  }
  cat("\n")
}

if (length(missing_topics)) {
  cat("Documentation topics missing from _pkgdown.yml reference index:\n")
  for (topic in missing_topics) {
    cat(sprintf("- %s\n", topic))
  }
  cat(
    "\nAdd these topics (or a selector that matches them) to `reference:` blocks",
    "in _pkgdown.yml.\n",
    sep = " "
  )
  quit(status = 1L)
}

if (length(unknown_entries)) {
  quit(status = 1L)
}

cat(sprintf(
  "Pkgdown reference coverage is complete: %d public topics are indexed.\n",
  length(topics_public)
))
