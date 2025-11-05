# R/zzz.R

# Internal hooks & registrations -----------------------------------------------

#' @keywords internal
#' @useDynLib RdistNN, .registration = TRUE
#' @importFrom Rcpp evalCpp
#' @importFrom stats pbeta plogis rbinom rnorm runif
NULL

# Silence R CMD check for non-standard eval symbols used in examples/helpers
utils::globalVariables(c("xx"))

.onLoad <- function(libname, pkgname) {
  # ---- Register built-in distances (if registry present) ---------------------
  # These constructors do not require user parameters; setup() can compute any
  # data-dependent context later.
  if (exists("register_distance", mode = "function")) {
    # L2 (squared by default; ordering same as unsquared, but faster for kNN/kernels)
    try(register_distance(distance_l2(squared = TRUE)), silent = TRUE)
    # L1 (Manhattan)
    try(register_distance(distance_l1()), silent = TRUE)
    # Normalized L2

    # Scaled / normalized variants — safe to auto-register because scaling
    # is learned from data later in distance_setup().
    if (exists("distance_l2_scaled", mode = "function")) {
      # Default (sd) with squared distances pairs nicely with Gaussian kernels
      try(register_distance(distance_l2_scaled(rule = "sd", squared = TRUE)), silent = TRUE)
      # A robust default many folks like
      try(register_distance(distance_l2_scaled(rule = "mad", squared = TRUE)), silent = TRUE)
    }
    if (exists("distance_l1_scaled", mode = "function")) {
      try(register_distance(distance_l1_scaled(rule = "sd")), silent = TRUE)
      try(register_distance(distance_l1_scaled(rule = "mad")), silent = TRUE)
    }
  }

  # (Optional) Pre-register distance-based weight rules that need no params.
  # Most weighting rules (e.g., Gaussian kernel) need a bandwidth h, so we
  # do NOT auto-register them here.
}

.onUnload <- function(libpath) {
  # Ensure native symbols are cleanly unloaded
  library.dynam.unload("RdistNN", libpath)
}

# (Optional) package attach message; keep silent by default
# .onAttach <- function(libname, pkgname) {
#   if (exists("list_distances", mode = "function")) {
#     ds <- try(list_distances(), silent = TRUE)
#     if (is.character(ds) && length(ds)) {
#       packageStartupMessage("RdistNN distances: ", paste(ds, collapse = ", "))
#     }
#   }
# }
