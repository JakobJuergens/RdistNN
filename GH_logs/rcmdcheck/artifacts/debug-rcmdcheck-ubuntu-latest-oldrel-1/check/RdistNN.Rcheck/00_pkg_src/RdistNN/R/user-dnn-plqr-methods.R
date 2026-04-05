# Methods for plqr_result objects.

#' Print a `plqr_result` object
#'
#' @param x An object returned by [plqr_dnn()].
#' @param ... Not used.
#' @return `x`, invisibly.
#' @export
print.plqr_result <- function(x, ...) {
  cat("<plqr_result>\n", sep = "")
  cat("  data  : n=", x$meta$n, ", p=", x$meta$p, ", p_d=", x$meta$p_d,
    ", m=", x$meta$m, "\n",
    sep = ""
  )
  cat("  model : tau=", x$tau, ", beta=[",
    paste0(round(x$beta, 4), collapse = ", "),
    "]\n",
    sep = ""
  )
  cat("  DNN   : method=", x$meta$method, ", k=", x$meta$k, ", s=", x$meta$s,
    ", distance=", x$meta$distance, "\n",
    sep = ""
  )
  cat("  est   : moment[,1] = ",
    paste0(round(x$moment[, 1], 4), collapse = ", "),
    if (!is.null(x$se)) "  (SEs available)\n" else "  (SEs not requested)\n",
    sep = ""
  )
  invisible(x)
}
