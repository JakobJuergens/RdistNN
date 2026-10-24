# Registry keys share a boundary without imposing syntactic-name restrictions.
.validate_registry_name <- function(name, argument = "name") {
  if (!is.character(name) || length(name) != 1L || is.na(name) || !nzchar(name)) {
    .abort_bad_arg(argument, "must be a non-missing, non-empty character scalar")
  }
  name
}

# Every accepted key is discoverable, including dot-prefixed names.
.registry_names <- function(registry) sort(ls(registry, all.names = TRUE))
