# R/zzz.R

#' Internal hooks
#' @keywords internal
#' @useDynLib RdistNN, .registration = TRUE
#' @importFrom Rcpp evalCpp

.onLoad <- function(libname, pkgname) {
  ns <- asNamespace(pkgname)

  # Coerce a raw list spec into a dnn_distance for robust registration
  as_dnn_spec <- function(spec) {
    if (inherits(spec, "dnn_distance")) {
      return(spec)
    }
    mk <- get0("new_distance", envir = ns, inherits = TRUE)
    if (!is.function(mk)) {
      return(spec)
    } # registry may accept raw specs
    nm <- tryCatch(spec$name, error = function(e) "unnamed")
    args <- list(name = nm, eval = spec$eval)
    if (!is.null(spec$eval_batch)) args$eval_batch <- spec$eval_batch
    if (!is.null(spec$setup)) args$setup <- spec$setup
    do.call(mk, args)
  }

  .register_builtin_distances()
}
