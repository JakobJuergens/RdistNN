test_that("README simulation guidance excludes stale scientific claims", {
  root <- .rdistnn_source_root(required = c("README.Rmd", "README.md"))
  stale_phrases <- c(
    "use n_reps = 500 for paper results",
    "comparisons isolate the latter"
  )

  for (filename in c("README.Rmd", "README.md")) {
    text <- paste(
      readLines(file.path(root, filename), warn = FALSE),
      collapse = " "
    )
    text <- tolower(gsub("[[:space:]]+", " ", text))
    for (phrase in stale_phrases) {
      expect_false(grepl(phrase, text, fixed = TRUE), info = filename)
    }
  }
})
