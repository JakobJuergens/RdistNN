test_that("notation vignette exposes current package and paper conventions", {
  root <- .rdistnn_source_root(required = c(
    "README.Rmd",
    "README.md",
    "_pkgdown.yml",
    "vignettes",
    "markdown_notes/rdistnn_roadmap.md"
  ))
  vignette_path <- file.path(root, "vignettes", "notation-and-conventions.Rmd")

  expect_true(file.exists(vignette_path))
  if (!file.exists(vignette_path)) {
    return(invisible())
  }

  vignette_text <- paste(readLines(vignette_path, warn = FALSE), collapse = "\n")
  required_sections <- c(
    "Localization and DNN conventions",
    "Estimator targets and signals",
    "CATE nuisance conventions",
    "PLQR target and nuisance conventions",
    "Moment diagnostics versus coefficient inference",
    "Reading fitted objects",
    "Matrix and data-frame interfaces"
  )
  required_names <- c(
    "`X`", "`x0`", "`X0`", "`s`", "`k`", "`n_eff`",
    "`theta_eval`", "`theta`", "`g`", "`mD`", "`density`",
    "`score`", "`moment`", "`moment_se`", "`theta_se`",
    "`theta_vcov`", "`theta_omega`"
  )

  for (section in required_sections) {
    expect_true(grepl(section, vignette_text, fixed = TRUE), info = section)
  }
  for (name in required_names) {
    expect_true(grepl(name, vignette_text, fixed = TRUE), info = name)
  }
  expect_true(grepl(
    "supplied externally",
    vignette_text,
    fixed = TRUE
  ))
  expect_true(grepl(
    "does not estimate a coefficient\\s+standard error",
    vignette_text,
    perl = TRUE
  ))

  pkgdown_text <- paste(
    readLines(file.path(root, "_pkgdown.yml"), warn = FALSE),
    collapse = "\n"
  )
  expect_true(grepl("- notation-and-conventions", pkgdown_text, fixed = TRUE))

  for (readme in c("README.Rmd", "README.md")) {
    readme_text <- paste(
      readLines(file.path(root, readme), warn = FALSE),
      collapse = "\n"
    )
    expect_true(grepl(
      'vignette("notation-and-conventions", package = "RdistNN")',
      readme_text,
      fixed = TRUE
    ))
  }

  roadmap_text <- paste(
    readLines(
      file.path(root, "markdown_notes", "rdistnn_roadmap.md"),
      warn = FALSE
    ),
    collapse = "\n"
  )
  expect_false(grepl("**Public notation page.**", roadmap_text, fixed = TRUE))
})

test_that("vignettes expose hub-and-spoke navigation", {
  root <- .rdistnn_source_root(required = c("vignettes"))
  vignette_dir <- file.path(root, "vignettes")
  hub_path <- file.path(vignette_dir, "notation-and-conventions.Rmd")
  hub_text <- paste(readLines(hub_path, warn = FALSE), collapse = "\n")

  hub_sections <- c(
    "How the vignette collection fits together",
    "Find information by task",
    "Where package objects come from",
    "Names that change meaning across estimators"
  )
  for (section in hub_sections) {
    expect_true(grepl(section, hub_text, fixed = TRUE), info = section)
  }

  spokes <- c(
    "getting-started",
    "conditional-z-estimation",
    "dnn-weights",
    "kernel-order-tuning",
    "diagnostics-reporting",
    "visualising-fits",
    "extending",
    "jackknife-gen-ustats",
    "local-iv-estimation",
    "paper-affine-replication"
  )
  for (spoke in spokes) {
    spoke_call <- sprintf('vignette("%s", package = "RdistNN")', spoke)
    expect_true(grepl(spoke_call, hub_text, fixed = TRUE), info = spoke_call)
    expect_true(grepl(
      sprintf("](%s.html)", spoke),
      hub_text,
      fixed = TRUE
    ), info = sprintf("clickable hub link to %s", spoke))

    spoke_path <- file.path(vignette_dir, paste0(spoke, ".Rmd"))
    spoke_text <- paste(readLines(spoke_path, warn = FALSE), collapse = "\n")
    expect_true(grepl(
      'vignette("notation-and-conventions", package = "RdistNN")',
      spoke_text,
      fixed = TRUE
    ), info = spoke_path)
    expect_true(grepl(
      "](notation-and-conventions.html)",
      spoke_text,
      fixed = TRUE
    ), info = sprintf("clickable notation link in %s", spoke_path))
  }

  continue_hits <- gregexpr("Continue with", hub_text, fixed = TRUE)[[1L]]
  continue_count <- if (identical(continue_hits, -1L)) 0L else length(continue_hits)
  expect_gte(continue_count, 5L)
})
