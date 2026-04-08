# Robustly extract registry names from distance_names() output
registry_names <- function(ds) {
  if (is.character(ds)) {
    return(ds)
  }
  if (is.data.frame(ds) && "name" %in% names(ds)) {
    return(as.character(ds$name))
  }
  if (is.list(ds) && !is.null(names(ds))) {
    return(names(ds))
  }
  if (is.list(ds) && length(ds) &&
    all(vapply(ds, function(x) is.list(x) && !is.null(x$name), logical(1)))) {
    return(vapply(ds, `[[`, character(1), "name"))
  }
  character()
}

# Coerce a raw spec (list) into a proper dnn_distance using new_distance()
as_dnn_spec <- function(spec, fallback_name = "unnamed") {
  if (inherits(spec, "dnn_distance")) {
    return(spec)
  }
  mk <- get0("new_distance", mode = "function", inherits = TRUE)
  if (!is.function(mk)) {
    testthat::skip("new_distance() not available to build dnn_distance specs.")
  }
  nm <- tryCatch(spec$name, error = function(e) NULL)
  if (is.null(nm) || !nzchar(nm)) nm <- fallback_name

  args <- list(name = nm)

  # Required: eval
  ev <- tryCatch(spec$eval, error = function(e) NULL)
  if (!is.null(ev)) {
    args$eval <- ev
  } else {
    testthat::fail(sprintf("Spec '%s' lacks an 'eval' function.", nm))
  }

  # Optional: eval_batch, setup
  eb <- tryCatch(spec$eval_batch, error = function(e) NULL)
  if (!is.null(eb)) args$eval_batch <- eb
  st <- tryCatch(spec$setup, error = function(e) NULL)
  if (!is.null(st)) args$setup <- st

  do.call(mk, args)
}

# Register helper:
#  - finds registry API
#  - coerces to dnn_distance
#  - registers via register_distance(name, spec)
#  - verifies success via distance_names()
ensure_registered <- function(named_specs) {
  reg <- get0("register_distance", mode = "function", inherits = TRUE)
  if (!is.function(reg)) {
    reg <- tryCatch(getFromNamespace("register_distance", "RdistNN"),
      error = function(e) NULL
    )
  }
  lister <- get0("distance_names", mode = "function", inherits = TRUE)
  if (!is.function(lister)) {
    lister <- tryCatch(getFromNamespace("distance_names", "RdistNN"),
      error = function(e) NULL
    )
  }
  if (!is.function(reg) || !is.function(lister)) {
    testthat::skip("Registry API not available in this test run.")
  }

  for (nm in names(named_specs)) {
    raw <- named_specs[[nm]]
    spec <- as_dnn_spec(raw, fallback_name = nm)
    name <- spec$name

    success <- FALSE
    try(
      {
        reg(name, spec)
        if (name %in% registry_names(lister())) success <- TRUE
      },
      silent = TRUE
    )

    if (!success) {
      testthat::fail(sprintf("Failed to register '%s' via register_distance(name, spec).", name))
    }
  }
  invisible(NULL)
}
