# Methods for tdnn_np_result objects.

#' Print a `tdnn_np_result` object
#'
#' @param x An object returned by [np_tdnn()].
#' @param ... Not used.
#' @return `x`, invisibly.
#' @export
print.tdnn_np_result <- function(x, ...) {
  cat("<tdnn_np_result>\n", sep = "")
  cat("  data  : n=", x$meta$n, ", p=", x$meta$p, ", m=", x$meta$m, "\n", sep = "")
  cat("  TDNN  : method=", x$meta$method,
    ", s1=", x$meta$s1, ", s2=", x$meta$s2,
    ", dim=", x$meta$dim,
    ", distance=", x$meta$distance, "\n",
    sep = ""
  )
  if (!is.null(x$se)) {
    cat("  est   : mu[1:min(5,m)] = ",
      paste0(round(utils::head(x$mu, 5L), 4L), collapse = ", "),
      "  (SEs available)\n",
      sep = ""
    )
  } else {
    cat("  est   : mu[1:min(5,m)] = ",
      paste0(round(utils::head(x$mu, 5L), 4L), collapse = ", "),
      "  (SEs not requested)\n",
      sep = ""
    )
  }
  invisible(x)
}
