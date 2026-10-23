# Select executable R examples, preserving source line numbers. Inline prose
# may explain a supported general-k boundary without changing a default call.
.rdistnn_documented_r_examples <- function(path, lines) {
  selected <- rep(FALSE, length(lines))
  active <- FALSE
  depth <- 0L
  for (i in seq_along(lines)) {
    line <- lines[i]
    if (grepl("\\.R$", path)) {
      if (grepl("^#' @examples", line)) active <- TRUE
      else if (grepl("^#' @", line)) active <- FALSE
      selected[i] <- active && grepl("^#'", line)
    } else if (grepl("\\.Rd$", path)) {
      if (grepl("^\\\\examples\\{", line)) active <- TRUE
      if (active) {
        selected[i] <- TRUE
        chars <- strsplit(line, "", fixed = TRUE)[[1L]]
        depth <- depth + sum(chars == "{") - sum(chars == "}")
        if (depth == 0L) active <- FALSE
      }
    } else {
      if (grepl("^```\\s*(\\{r(?:[ ,}]|$)|r\\s*$)", line, perl = TRUE)) {
        active <- TRUE
      } else if (grepl("^```\\s*$", line, perl = TRUE)) {
        active <- FALSE
      } else selected[i] <- active
    }
  }
  ifelse(selected, lines, "")
}

test_that("the default-example audit distinguishes executable calls from general-k prose", {
  prose <- "For `k = 2`, one closer observation and two boundary observations are averaged."
  markdown <- c(prose, "```{r example}", "weights_dnn(X, x0 = 0, k = 2L, s = 3L)", "```", prose)
  expect_identical(.rdistnn_documented_r_examples("example.Rmd", markdown),
    c("", "", markdown[3L], "", ""))
  roxygen <- c("#' For k = 2, ...", "#' @examples", "#' weights_dnn(X, x0 = 0, k = 2L, s = 3L)", "#' @export")
  expect_identical(.rdistnn_documented_r_examples("example.R", roxygen),
    c("", roxygen[2:3], ""))
  rd <- c("For k = 2, ...", "\\examples{", "weights_dnn(X, x0 = 0, k = 2L, s = 3L)", "}", "For k = 2, ...")
  expect_identical(.rdistnn_documented_r_examples("example.Rd", rd),
    c("", rd[2:4], ""))
})

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
    doc_lines[[path]] <- .rdistnn_documented_r_examples(path, lines)
  }
  for (path in r_files) {
    lines <- readLines(path, warn = FALSE)
    doc_lines[[path]] <- .rdistnn_documented_r_examples(path, lines)
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
