#!/usr/bin/env Rscript

parse_desc_packages <- function(desc_path = "DESCRIPTION") {
  if (!file.exists(desc_path)) {
    stop(sprintf("DESCRIPTION not found at '%s'.", desc_path), call. = FALSE)
  }

  dcf <- read.dcf(desc_path)
  pkg_name <- trimws(dcf[1L, "Package"])
  fields <- c("Depends", "Imports", "Suggests", "Enhances")

  declared <- pkg_name
  for (field in fields) {
    if (!field %in% colnames(dcf)) {
      next
    }
    raw <- dcf[1L, field]
    if (is.na(raw) || !nzchar(raw)) {
      next
    }
    vals <- trimws(strsplit(raw, ",", fixed = TRUE)[[1L]])
    vals <- sub("\\s*\\(.*\\)\\s*$", "", vals, perl = TRUE)
    vals <- vals[nzchar(vals)]
    declared <- c(declared, vals)
  }

  unique(declared)
}

extract_pkg_refs <- function(line) {
  pattern <- "\\b([A-Za-z][A-Za-z0-9.]*)::{1,2}"
  hits <- gregexpr(pattern, line, perl = TRUE)[[1L]]
  if (length(hits) == 1L && hits[1L] == -1L) {
    return(character(0))
  }
  tokens <- regmatches(line, gregexpr(pattern, line, perl = TRUE))[[1L]]
  unique(sub(":+$", "", tokens, perl = TRUE))
}

scan_test_refs <- function(tests_dir = "tests") {
  if (!dir.exists(tests_dir)) {
    return(data.frame(
      pkg = character(0),
      file = character(0),
      line = integer(0),
      stringsAsFactors = FALSE
    ))
  }

  files <- list.files(
    tests_dir,
    pattern = "\\.[Rr]$",
    recursive = TRUE,
    full.names = TRUE
  )

  rows <- vector("list", length = 0L)
  idx <- 1L

  for (path in files) {
    lines <- readLines(path, warn = FALSE, encoding = "UTF-8")
    for (i in seq_along(lines)) {
      pkgs <- extract_pkg_refs(lines[[i]])
      if (!length(pkgs)) {
        next
      }
      for (pkg in pkgs) {
        rows[[idx]] <- data.frame(
          pkg = pkg,
          file = path,
          line = i,
          stringsAsFactors = FALSE
        )
        idx <- idx + 1L
      }
    }
  }

  if (!length(rows)) {
    return(data.frame(
      pkg = character(0),
      file = character(0),
      line = integer(0),
      stringsAsFactors = FALSE
    ))
  }

  out <- do.call(rbind, rows)
  out[order(out$pkg, out$file, out$line), , drop = FALSE]
}

base_recommended <- c(
  "base", "compiler", "datasets", "graphics", "grDevices", "grid", "methods",
  "parallel", "splines", "stats", "stats4", "tcltk", "tools", "utils",
  "boot", "class", "cluster", "codetools", "foreign", "KernSmooth", "lattice",
  "MASS", "Matrix", "mgcv", "nlme", "nnet", "rpart", "spatial", "survival"
)

declared <- parse_desc_packages("DESCRIPTION")
refs <- scan_test_refs("tests")

if (!nrow(refs)) {
  cat("No pkg::/pkg::: references found under tests/. Nothing to check.\n")
  quit(status = 0L)
}

unknown <- sort(setdiff(unique(refs$pkg), c(declared, base_recommended)))

if (length(unknown)) {
  cat("Undeclared package references were found in tests/.\n")
  cat("Add these packages to DESCRIPTION (usually Suggests for tests):\n\n")
  for (pkg in unknown) {
    cat(sprintf("- %s\n", pkg))
    rows <- refs[refs$pkg == pkg, , drop = FALSE]
    locs <- unique(sprintf("  - %s:%d", rows$file, rows$line))
    cat(paste(locs, collapse = "\n"))
    cat("\n")
  }
  quit(status = 1L)
}

cat("All pkg::/pkg::: references in tests/ are declared in DESCRIPTION.\n")
