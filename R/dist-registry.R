# Simple in-memory registry (package namespace)
.distance_registry <- new.env(parent = emptyenv())

#' Register a distance object by name
#' @param distance a \code{dnn_distance}
#' @return registered name (invisibly)
#' @export
register_distance <- function(distance) {
  stopifnot(is_distance(distance))
  assign(distance$name, distance, envir = .distance_registry)
  invisible(distance$name)
}

#' Retrieve a registered distance
#' @param name character scalar
#' @return a \code{dnn_distance}
#' @export
get_distance <- function(name) {
  if (!exists(name, envir = .distance_registry, inherits = FALSE)) {
    stop("No distance registered under '", name, "'.")
  }
  get(name, envir = .distance_registry, inherits = FALSE)
}

#' List registered distances
#' @return character vector of names
#' @export
list_distances <- function() sort(ls(envir = .distance_registry, all.names = FALSE))
