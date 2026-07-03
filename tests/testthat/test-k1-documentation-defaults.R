test_that("documentation examples use k = 1 for DNN defaults", {
  root <- .rdistnn_source_root(required = c("README.md", "vignettes", "man", "R"))

  doc_files <- c(
    file.path(root, "README.md"),
    list.files(
      file.path(root, "vignettes"),
      pattern = "\\.Rmd$",
      full.names = TRUE
    ),
    list.files(
      file.path(root, "man"),
      pattern = "\\.Rd$",
      full.names = TRUE
    )
  )
  r_files <- list.files(file.path(root, "R"), pattern = "\\.R$", full.names = TRUE)

  doc_lines <- list()
  for (path in doc_files) {
    lines <- readLines(path, warn = FALSE)
    doc_lines[[path]] <- lines
  }
  for (path in r_files) {
    lines <- readLines(path, warn = FALSE)
    doc_lines[[path]] <- grep("^#'", lines, value = TRUE)
  }

  hits <- unlist(Map(function(path, lines) {
    matches <- gregexpr("\\bk\\s*=\\s*([0-9]+)L?\\b", lines, perl = TRUE)
    idx <- which(lengths(regmatches(lines, matches)) > 0L)
    if (length(idx) == 0L) {
      return(character())
    }
    offending <- idx[vapply(idx, function(i) {
      vals <- regmatches(lines[i], list(matches[[i]]))[[1L]]
      any(!grepl("\\bk\\s*=\\s*1L?\\b", vals, perl = TRUE))
    }, logical(1))]
    if (length(offending) == 0L) {
      return(character())
    }
    sprintf("%s:%d:%s", path, offending, lines[offending])
  }, names(doc_lines), doc_lines), use.names = FALSE)

  expect_equal(hits, character())
})
