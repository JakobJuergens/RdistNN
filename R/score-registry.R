# Simple in-memory registry (package namespace)
.score_registry <- new.env(parent = emptyenv())

#' Register a score object by name
#' @param score a \code{ustat_score}
#' @return registered name (invisibly)
#' @export
register_score <- function(score) {
  stopifnot(is_score(score))
  assign(score$name, score, envir = .score_registry)
  invisible(score$name)
}

#' Retrieve a registered score
#' @param name character scalar
#' @return a \code{ustat_score}
#' @export
get_score <- function(name) {
  if (!exists(name, envir = .score_registry, inherits = FALSE)) {
    stop("No score registered under '", name, "'.")
  }
  get(name, envir = .score_registry, inherits = FALSE)
}

#' List registered scores
#' @return character vector of names
#' @export
list_scores <- function() sort(ls(envir = .score_registry, all.names = FALSE))
