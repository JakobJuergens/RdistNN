# R/rank-registry.R

.rank_registry <- new.env(parent = emptyenv())

#' Register a rank_weight by name
#' @export
register_rank <- function(rw) {
  stopifnot(is_rank_weight(rw))
  assign(rw$name, rw, envir = .rank_registry)
  invisible(rw$name)
}

#' Retrieve a registered rank_weight
#' @export
get_rank <- function(name) {
  if (!exists(name, envir = .rank_registry, inherits = FALSE)) {
    stop("No rank_weight registered under '", name, "'.")
  }
  get(name, envir = .rank_registry, inherits = FALSE)
}

#' List registered rank_weight names
#' @export
list_ranks <- function() sort(ls(envir = .rank_registry, all.names = FALSE))
