# Observe namespace calls without replacing constructor environments or attached
# functions. The helper frame owns cleanup on success and error.
.count_package_calls <- function(fun, expr) {
  calls <- 0L
  count <- function() {
    calls <<- calls + 1L
  }
  namespace <- asNamespace("RdistNN")
  # Pass the closure itself; a symbolic tracer name is looked up in the
  # observed function's frame rather than this helper's frame.
  suppressMessages(invisible(do.call(trace,
    list(what = fun, tracer = count, where = namespace, print = FALSE))))
  on.exit(suppressMessages(untrace(fun, where = namespace)), add = TRUE)
  force(expr)
  calls
}

# Public geometry callbacks with independently observable setup/evaluation state.
.counted_l1_distance <- function(on_setup = function() NULL,
                                 on_eval = function() NULL) {
  counts <- new.env(parent = emptyenv())
  counts$setup <- counts$eval <- 0L
  metric <- new_distance("counted_l1",
    setup = function(X) {
      counts$setup <- counts$setup + 1L
      on_setup()
      list()
    },
    eval_batch = function(X, X0, ctx = NULL) {
      counts$eval <- counts$eval + 1L
      on_eval()
      abs(outer(X[, 1L], X0[, 1L], "-"))
    })
  list(metric = metric, counts = counts)
}
