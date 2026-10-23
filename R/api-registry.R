# Registry keys share a boundary without imposing syntactic-name restrictions.
.validate_registry_name <- function(name) {
  if (!is.character(name) || length(name) != 1L || is.na(name) || !nzchar(name)) {
    .abort_bad_arg("name", "must be a non-missing, non-empty character scalar")
  }
  name
}
