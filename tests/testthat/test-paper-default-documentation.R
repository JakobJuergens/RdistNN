test_that("public estimator docs state the finite-sample 0.4 rule", {
  root <- .rdistnn_source_root(required = c("R", "man"))
  stems <- c("npr", "cate", "plqr", "gmm")

  for (stem in stems) {
    source_text <- paste(readLines(
      file.path(root, "R", paste0("user-dnn-", stem, ".R")),
      warn = FALSE
    ), collapse = "\n")
    help_name <- switch(stem, npr = "np_dnn", paste0(stem, "_dnn"))
    help_text <- paste(readLines(
      file.path(root, "man", paste0(help_name, ".Rd")),
      warn = FALSE
    ), collapse = "\n")

    for (text in c(source_text, help_text)) {
      expect_true(grepl(
        "0.4 n^{(k_x+1)/(k_x+2)}",
        text,
        fixed = TRUE
      ))
      expect_false(grepl(
        "s = \\lfloor n^{(k_x+1)/(k_x+2)}",
        text,
        fixed = TRUE
      ))
    }
  }
})
