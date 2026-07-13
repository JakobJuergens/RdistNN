test_that("CITATION metadata is valid and release-ready", {
  root <- .rdistnn_source_root(required = c("DESCRIPTION", "inst/CITATION"))
  description <- as.list(read.dcf(file.path(root, "DESCRIPTION"))[1, ])
  citation_path <- file.path(root, "inst", "CITATION")
  read_citation <- getFromNamespace("readCitationFile", "utils")

  expect_no_error(entries <- read_citation(citation_path, meta = description))
  expect_s3_class(entries, "citation")
  expect_s3_class(entries, "bibentry")
  expect_length(entries, 3L)
  expect_true(all(vapply(seq_along(entries), function(i) {
    inherits(entries[i], "bibentry")
  }, logical(1))))

  fields <- lapply(seq_along(entries), function(i) unclass(entries[i])[[1]])
  types <- vapply(fields, function(entry) {
    attr(entry, "bibtype", exact = TRUE)
  }, character(1))
  normalize_dcf_text <- function(x) trimws(gsub("[[:space:]]+", " ", x))

  manual <- fields[[which(types == "Manual")]]
  expect_identical(
    normalize_dcf_text(manual$title),
    normalize_dcf_text(description$Title)
  )
  expect_identical(manual$year, "2026")
  expect_identical(manual$url, "https://jakobjuergens.github.io/RdistNN/")

  titles <- vapply(fields, `[[`, character(1), "title")
  conditional_z <- fields[[match(
    "Robust Inference for Conditional Z-Estimation via Distributional Nearest Neighbors",
    titles
  )]]
  expect_identical(conditional_z$year, "2026")
  expect_identical(
    conditional_z$url,
    "https://jakobjuergens.com/research/#conditional-z-estimation-dnn"
  )

  jackknife <- fields[[match(
    "Jackknife Variance Estimation for Hájek-Dominated Generalized U-Statistics",
    titles
  )]]
  expect_identical(jackknife$year, "2025")
  expect_identical(jackknife$doi, "10.48550/arXiv.2509.12356")
  expect_identical(jackknife$url, "https://arxiv.org/abs/2509.12356")

  source_text <- paste(readLines(citation_path, warn = FALSE), collapse = "\n")
  stale_patterns <- c(
    "replace with actual",
    "placeholder",
    "github.com/JakobJuergens/Cond_Z_Estimation",
    "RdistNN: Modular distances and DNN weight engines for U-statistics"
  )
  expect_false(any(vapply(tolower(stale_patterns), grepl, logical(1),
    x = tolower(source_text), fixed = TRUE
  )))
})
