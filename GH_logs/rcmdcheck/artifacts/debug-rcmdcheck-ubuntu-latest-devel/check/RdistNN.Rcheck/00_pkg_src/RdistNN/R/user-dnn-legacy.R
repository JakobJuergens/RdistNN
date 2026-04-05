# Legacy compatibility methods.

#' Print a `cate_dnn_result` object
#'
#' Legacy print method retained for compatibility with older object classes.
#'
#' @inheritParams print.cate_result
#' @return `x`, invisibly.
#' @export
print.cate_dnn_result <- function(x, ...) {
  cat("<cate_dnn_result>\n", sep = "")
  cat("  data  : n=", x$meta$n, ", p=", x$meta$p, ", m=", x$meta$m, "\n", sep = "")
  cat("  dnn   : k=", x$meta$k, ", s=", x$meta$s, ", method=", x$meta$method, "\n", sep = "")
  cat("  rank  : ", x$meta$weight_name, "; distance: ", x$meta$distance_name, "\n", sep = "")
  cat("  components: ", paste(names(x), collapse = ", "), "\n", sep = "")
  invisible(x)
}
