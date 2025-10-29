# R/wdist-registry.R
# -----------------------------------------------------------------------------
# In-memory registry for distance-weight specs (dist_weight)
# Mirrors rank-registry.R
# -----------------------------------------------------------------------------

.wdist_registry <- new.env(parent = emptyenv())

#' Register a dist_weight by name
#' @param dw a \code{dist_weight}
#' @return registered name (invisibly)
#' @export
register_wdist <- function(dw) {
  stopifnot(is_dist_weight(dw))
  assign(dw$name, dw, envir = .wdist_registry)
  invisible(dw$name)
}

#' Retrieve a registered dist_weight
#' @param name character scalar
#' @return a \code{dist_weight}
#' @export
get_wdist <- function(name) {
  if (!exists(name, envir = .wdist_registry, inherits = FALSE)) {
    stop("No dist_weight registered under '", name, "'.")
  }
  get(name, envir = .wdist_registry, inherits = FALSE)
}

#' List registered dist_weight names
#' @return character vector of names
#' @export
list_wdists <- function() sort(ls(envir = .wdist_registry, all.names = FALSE))
