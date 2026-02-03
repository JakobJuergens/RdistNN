args <- commandArgs(trailingOnly = TRUE)
error_on_lint <- if (length(args) > 0L) as.logical(args[[1L]]) else TRUE

lints <- lintr::lint_package()
if (length(lints) > 0L) {
  print(lints)
}
if (isTRUE(error_on_lint) && length(lints) > 0L) {
  quit(status = 1L)
}
