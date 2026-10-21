test_that("notation documents the current input and uncertainty fields", {
  root <- .rdistnn_source_root(required = c("vignettes"))
  text <- paste(readLines(file.path(root, "vignettes", "notation-and-conventions.Rmd"),
                          warn = FALSE), collapse = "\n")
  fields <- c("X", "x0", "X0", "s", "k", "n_eff", "theta_eval", "theta",
              "g", "mD", "density", "moment", "moment_se", "theta_se",
              "theta_vcov", "theta_omega")
  for (field in fields) {
    expect_true(grepl(paste0("`", field, "`"), text, fixed = TRUE), info = field)
  }
  for (estimator in c("np_dnn", "cate_dnn", "clate_dnn", "plqr_dnn", "gmm_dnn")) {
    expect_true(grepl(paste0(estimator, "()"), text, fixed = TRUE), info = estimator)
  }
})

test_that("articles are registered and their relative article links resolve", {
  root <- .rdistnn_source_root(required = c("vignettes", "_pkgdown.yml"))
  paths <- list.files(file.path(root, "vignettes"), pattern = "\\.Rmd$", full.names = TRUE)
  slugs <- sub("\\.Rmd$", "", basename(paths))
  config <- paste(readLines(file.path(root, "_pkgdown.yml"), warn = FALSE), collapse = "\n")
  for (i in seq_along(paths)) {
    expect_match(config, paste0("(?m)^\\s*-\\s+", slugs[i], "\\s*$"), perl = TRUE)
    text <- paste(readLines(paths[i], warn = FALSE), collapse = "\n")
    links <- regmatches(text, gregexpr("\\]\\([^)]*\\.html(?:#[^)]*)?\\)", text, perl = TRUE))[[1L]]
    targets <- sub("^\\]\\((.*)\\)$", "\\1", links)
    targets <- targets[!grepl("^(https?:)?//", targets)]
    destinations <- sub("\\.html(?:#.*)?$", "", basename(targets))
    expect_true(all(destinations %in% slugs), info = basename(paths[i]))
  }
})
