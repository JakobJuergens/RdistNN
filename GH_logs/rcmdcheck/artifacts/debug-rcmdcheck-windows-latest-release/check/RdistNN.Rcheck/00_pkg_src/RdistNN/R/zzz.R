#' Internal hooks
#' @keywords internal
#' @useDynLib RdistNN, .registration = TRUE
#' @importFrom Rcpp evalCpp

.onLoad <- function(libname, pkgname) {
  .register_builtin_distances()
}
