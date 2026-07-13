# Internal compiled-engine dispatch.

.rdistnn_compiled_enabled <- function() {
  isTRUE(getOption("RdistNN.use_compiled", TRUE))
}

.rdistnn_compiled_fun <- function(name) {
  if (!.rdistnn_compiled_enabled()) {
    return(NULL)
  }
  fun <- get0(name, mode = "function", inherits = TRUE)
  if (is.null(fun)) {
    return(NULL)
  }
  symbol <- paste0("_RdistNN_", name)
  if (!isTRUE(is.loaded(symbol, PACKAGE = "RdistNN"))) {
    return(NULL)
  }
  fun
}
