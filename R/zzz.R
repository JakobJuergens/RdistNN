# R/zzz.R

#' @useDynLib RdistNN, .registration = TRUE
#' @importFrom Rcpp evalCpp
NULL

.onLoad <- function(libname, pkgname) {
  # ---- Register built-in distances ----
  # These wrap compiled fast paths if available (else R fallbacks).
  # They are safe to call here because they take no required args.
  if (exists("register_distance", mode = "function")) {
    # L2 (default: squared for kNN/kernels)
    try(register_distance(distance_l2(squared = TRUE)), silent = TRUE)
    # L1
    try(register_distance(distance_l1()), silent = TRUE)
  }

  # ---- Optionally register pre-implemented scores ----
  # Some scores may require parameters (e.g., tau), so we only auto-register
  # those whose constructors take *no required arguments*. Others should be
  # registered by the user explicitly (or via a small helper you expose).
  if (exists("register_score", mode = "function") &&
    exists("is_score", mode = "function")) {
    safe_register_score <- function(fname) {
      if (!exists(fname, mode = "function")) {
        return(invisible(FALSE))
      }
      f <- get(fname, inherits = TRUE)
      fm <- formals(f)
      # Count required (non-default) arguments
      required <- if (is.null(fm)) {
        0L
      } else {
        sum(vapply(fm, function(x) identical(x, quote(expr = )), logical(1L)))
      }
      if (required == 0L) {
        obj <- try(f(), silent = TRUE)
        if (inherits(obj, "ustat_score")) {
          try(register_score(obj), silent = TRUE)
          return(invisible(TRUE))
        }
      }
      invisible(FALSE)
    }

    # Try common zero-arg constructors if you provide them:
    # e.g., `score_cate()` with sensible defaults.
    safe_register_score("score_cate")
    safe_register_score("score_plqr") # only if zero-arg exists
    safe_register_score("score_plqr_component") # only if zero-arg exists
  }
}

# (Optional) A lightweight attach message; comment out if you prefer silence.
# .onAttach <- function(libname, pkgname) {
#  if (exists("list_distances", mode = "function")) {
#    ds <- try(list_distances(), silent = TRUE)
#    if (is.character(ds) && length(ds)) {
#      packageStartupMessage("RdistNN distances: ", paste(ds, collapse = ", "))
#    }
#  }
# }
