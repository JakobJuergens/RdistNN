# The packaged affine study has explicit source and output paths. Do not carry
# the paper's superseded CATE manifest/orchestration infrastructure into it.
sim_find_project_root <- function(start = getwd()) {
  stop("Use the installed affine launcher with an explicit --output-dir; paper-checkout defaults are unavailable in the package.",
    call. = FALSE)
}
