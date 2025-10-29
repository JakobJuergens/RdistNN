# Simple in-memory registry (package namespace)
.score_registry <- new.env(parent = emptyenv())

#' Register a score object by name
#' @export
register_score <- function(score) {
  stopifnot(is_score(score))
  assign(score$name, score, envir = .score_registry)
  invisible(score$name)
}

#' Retrieve a registered score
#' @export
get_score <- function(name) {
  if (!exists(name, envir = .score_registry, inherits = FALSE)) {
    stop("No score registered under '", name, "'.")
  }
  get(name, envir = .score_registry, inherits = FALSE)
}

#' List registered scores
#' @export
list_scores <- function() sort(ls(envir = .score_registry, all.names = FALSE))
