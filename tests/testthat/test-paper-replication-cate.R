.cate_repl_path <- function(...) {
  path <- system.file("paper-replication", "cate", ..., package = "RdistNN")
  if (!nzchar(path)) {
    stop("Could not locate installed CATE paper-replication files.", call. = FALSE)
  }
  path
}

.source_cate_benchmarks <- function() {
  env <- new.env(parent = globalenv())
  sys.source(.cate_repl_path("benchmarks.R"), envir = env)
  env
}

.collect_named_calls <- function(expr, target) {
  found <- list()
  visit <- function(node) {
    if (missing(node)) {
      return(invisible(NULL))
    }
    if (!is.call(node)) {
      return(invisible(NULL))
    }
    label <- paste(deparse(node[[1L]], width.cutoff = 500L), collapse = "")
    if (identical(label, target)) {
      found[[length(found) + 1L]] <<- node
    }
    if (length(node) >= 2L) {
      for (i in 2L:length(node)) {
        if (identical(paste(deparse(node[[i]]), collapse = ""), "")) next
        visit(node[[i]])
      }
    }
    invisible(NULL)
  }
  visit(expr)
  found
}

.call_argument_names <- function(call) {
  args <- as.list(call)[-1L]
  nm <- names(args)
  if (is.null(nm)) rep("", length(args)) else nm
}

test_that("installed CATE local-linear GRF benchmark returns finite inference", {
  skip_if_not_installed("grf")
  env <- .source_cate_benchmarks()
  args <- list(
    n = 160L, kx = 2L, design = "dgp1_linear", alpha = 0.05,
    seed = 817L, sigma = 0.5, clip = 1e-3,
    num.trees = 200L, num.threads = 1L
  )

  out <- do.call(env$one_rep_llf, args)
  draw <- do.call(env$.bench_draw, args[c(
    "n", "kx", "design", "sigma", "clip", "seed"
  )])
  tau_hat <- unname(out[["bias"]]) + draw$tau0
  variance_hat <- unname(out[["se"]])^2

  expect_true(is.finite(tau_hat))
  expect_true(is.finite(variance_hat))
  expect_gte(variance_hat, 0)
  expect_identical(
    attr(out, "method_variant"),
    "grf_local_linear_internal_nuisance"
  )
  expect_identical(attr(out, "forest_nuisance"), "internal")
  expect_false(attr(out, "y_hat_supplied"))
  expect_false(attr(out, "w_hat_supplied"))
})

test_that("local-linear GRF options are prediction-time only", {
  env <- .source_cate_benchmarks()
  replication_body <- body(env$one_rep_llf)
  maintained_body <- body(getFromNamespace(".method_llf", "RdistNN"))

  replication_train <- .collect_named_calls(
    replication_body, "grf::causal_forest"
  )
  replication_predict <- .collect_named_calls(
    replication_body, "stats::predict"
  )
  maintained_train <- .collect_named_calls(
    maintained_body, "grf::causal_forest"
  )
  maintained_predict <- .collect_named_calls(
    maintained_body, "stats::predict"
  )

  expect_length(replication_train, 1L)
  expect_length(replication_predict, 1L)
  expect_length(maintained_train, 1L)
  expect_length(maintained_predict, 1L)

  forbidden_training <- c("ll.lambda", "ll.weight.penalty")
  intended_prediction <- c(
    "ll.lambda", "ll.weight.penalty", "linear.correction.variables"
  )
  for (training_call in c(replication_train, maintained_train)) {
    expect_false(any(
      forbidden_training %in% .call_argument_names(training_call)
    ))
  }
  for (prediction_call in c(replication_predict, maintained_predict)) {
    expect_true(all(
      intended_prediction %in% .call_argument_names(prediction_call)
    ))
  }
})
