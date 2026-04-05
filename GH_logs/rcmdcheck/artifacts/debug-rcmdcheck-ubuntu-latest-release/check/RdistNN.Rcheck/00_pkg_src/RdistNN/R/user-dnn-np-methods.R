# Methods for np_result objects.

#' Print an `np_result` object
#'
#' @param x An object returned by [np_dnn()].
#' @param ... Not used.
#' @return `x`, invisibly.
#' @export
print.np_result <- function(x, ...) {
  cat("<np_result>\n", sep = "")
  cat("  data  : n=", x$meta$n, ", p=", x$meta$p, ", m=", x$meta$m, "\n", sep = "")
  cat("  DNN   : method=", x$meta$method, ", k=", x$meta$k, ", s=", x$meta$s,
    ", distance=", x$meta$distance, "\n",
    sep = ""
  )
  if (!is.null(x$se)) {
    cat("  est   : mu[1:min(5,m)] = ",
      paste0(round(utils::head(x$mu, 5), 4), collapse = ", "),
      "  (SEs available)\n",
      sep = ""
    )
  } else {
    cat("  est   : mu[1:min(5,m)] = ",
      paste0(round(utils::head(x$mu, 5), 4), collapse = ", "),
      "  (SEs not requested)\n",
      sep = ""
    )
  }
  invisible(x)
}
