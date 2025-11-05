# R/zzz.R

#' Internal hooks
#' @keywords internal
#' @useDynLib RdistNN, .registration = TRUE
#' @importFrom Rcpp evalCpp

.onLoad <- function(libname, pkgname) {
  ns <- asNamespace(pkgname)

  # Coerce a raw list spec into a dnn_distance for robust registration
  as_dnn_spec <- function(spec) {
    if (inherits(spec, "dnn_distance")) return(spec)
    mk <- get0("new_distance", envir = ns, inherits = TRUE)
    if (!is.function(mk)) return(spec)  # registry may accept raw specs
    nm <- tryCatch(spec$name, error = function(e) "unnamed")
    args <- list(name = nm, eval = spec$eval)
    if (!is.null(spec$eval_batch)) args$eval_batch <- spec$eval_batch
    if (!is.null(spec$setup))      args$setup      <- spec$setup
    do.call(mk, args)
  }

  reg <- get0("register_distance", envir = ns, inherits = TRUE)
  reg_spec <- function(spec) {
    if (!is.function(reg)) return(invisible(NULL))
    spec <- as_dnn_spec(spec)
    nm   <- tryCatch(spec$name, error = function(e) NULL)
    fmls <- tryCatch(formals(args(reg)), error = function(e) NULL)
    n_args <- if (is.null(fmls)) 1L else length(fmls)
    if (n_args >= 2L && !is.null(nm)) try(reg(nm, spec), silent = TRUE) else try(reg(spec), silent = TRUE)
    invisible(NULL)
  }

  # Core
  if (is.function(get0("distance_l2", envir = ns))) {
    reg_spec(ns$distance_l2(squared = FALSE)) # "l2"
    reg_spec(ns$distance_l2(squared = TRUE))  # "l2_sq"
  }
  if (is.function(get0("distance_l1", envir = ns))) {
    reg_spec(ns$distance_l1()) # "l1"
  }

  # Scaled defaults
  if (is.function(get0("distance_l2_scaled", envir = ns))) {
    reg_spec(ns$distance_l2_scaled(rule = "sd",  squared = TRUE)) # "l2_scaled_sq[sd]"
    reg_spec(ns$distance_l2_scaled(rule = "mad", squared = TRUE)) # "l2_scaled_sq[mad]"
  }
  if (is.function(get0("distance_l1_scaled", envir = ns))) {
    reg_spec(ns$distance_l1_scaled(rule = "sd"))  # "l1_scaled[sd]"
    reg_spec(ns$distance_l1_scaled(rule = "mad")) # "l1_scaled[mad]"
  }
}
