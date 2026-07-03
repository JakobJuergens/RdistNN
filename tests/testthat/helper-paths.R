.rdistnn_source_root <- function(required = character()) {
  test_root <- normalizePath(
    file.path(testthat::test_path(), "..", ".."),
    mustWork = FALSE
  )
  candidates <- unique(normalizePath(c(
    test_root,
    file.path(test_root, "00_pkg_src", "RdistNN"),
    file.path(test_root, "RdistNN"),
    getwd(),
    file.path(getwd(), ".."),
    file.path(getwd(), "..", ".."),
    file.path(getwd(), "..", "..", "00_pkg_src", "RdistNN")
  ), mustWork = FALSE))

  for (cand in candidates) {
    if (!file.exists(file.path(cand, "DESCRIPTION"))) {
      next
    }
    required_paths <- file.path(cand, required)
    if (all(file.exists(required_paths) | dir.exists(required_paths))) {
      return(cand)
    }
  }

  testthat::skip("RdistNN source tree files are unavailable in this test context.")
}
