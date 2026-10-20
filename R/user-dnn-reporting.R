# Descriptive localization labels shared by fit, summary and plot reporting.
.dnn_localization_label <- function(meta) {
  orders <- if (!is.null(meta$s1) && !is.null(meta$s2)) {
    c(s1 = meta$s1, s2 = meta$s2)
  } else if (!is.null(meta$s)) {
    c(s = meta$s)
  } else {
    return("")
  }
  n <- meta$n
  valid_n <- is.numeric(n) && length(n) == 1L && is.finite(n) && n > 0
  paste(vapply(names(orders), function(name) {
    s <- orders[[name]]
    fraction <- if (valid_n) format(signif(s / n, 3L), trim = TRUE) else "unavailable"
    sprintf("%s = %s, %s/n = %s", name, s, name, fraction)
  }, character(1L)), collapse = "; ")
}

.dnn_print_localization <- function(meta) {
  label <- .dnn_localization_label(meta)
  if (nzchar(label)) cat("  localization: ", label, "\n", sep = "")
  invisible(NULL)
}

.dnn_print_diagnostic_range <- function(x, field, label) {
  if (!field %in% names(x)) return(invisible(NULL))
  values <- x[[field]][is.finite(x[[field]])]
  text <- if (length(values)) {
    paste(format(signif(range(values), 3L), trim = TRUE), collapse = " to ")
  } else {
    "unavailable"
  }
  cat("  ", label, ": ", text, "\n", sep = "")
  invisible(NULL)
}
