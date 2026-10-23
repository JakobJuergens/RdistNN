test_that("shared argument errors carry bounded facts and ordinary error semantics", {
  err <- tryCatch(RdistNN:::.validate_scalar_int(0, "k"), error = identity)
  expect_s3_class(err, "rdistnn_input_error")
  expect_identical(err$code, "argument_invalid")
  expect_identical(err$argument, "k")
  expect_match(err$actual, "0", fixed = TRUE)
  expect_null(err$call)
  expect_match(conditionMessage(err), "positive integer", fixed = TRUE)
  err <- tryCatch(RdistNN:::.abort_bad_arg("X0", "has incompatible columns",
    expected = "3 columns", actual = "2 columns"), error = identity)
  expect_identical(conditionMessage(err),
    "`X0` has incompatible columns. Expected: 3 columns. Got: 2 columns.")
  expect_identical(err$expected, "3 columns")
  expect_identical(err$actual, "2 columns")
})

test_that("unknown extras retain their remedy with a bounded literal list", {
  unknown <- paste0("unknown_", seq_len(200))
  unknown[1] <- paste0("brace{%}\n", paste(rep("z", 150), collapse = ""))
  got <- NULL
  withCallingHandlers(RdistNN:::.warn_unknown_return(unknown, c("W", "dx", "ranks")),
    warning = function(w) { got <<- w; invokeRestart("muffleWarning") })
  expect_s3_class(got, "rdistnn_unknown_return_warning")
  expect_equal(got$affected, 200L)
  expect_length(got$examples, 5L)
  expect_true(all(nchar(got$examples) <= 80L))
  expect_lt(nchar(conditionMessage(got)), 800L)
  expect_match(conditionMessage(got), "195 more", fixed = TRUE)
  expect_match(conditionMessage(got), "Valid options: W, dx, ranks", fixed = TRUE)
  expect_false(grepl("\n", got$examples[1], fixed = TRUE))
})

test_that("owned warning identities follow the base warning protocol", {
  for (warn in c(-1, 0, 1, 2)) {
    withr::local_options(warn = warn)
    got <- NULL
    withCallingHandlers(RdistNN:::.dnn_warn("partial_result", "Results are unavailable.",
      affected = 6L, examples = 1:6, hint = "Inspect status."), warning = function(w) {
        got <<- w
        invokeRestart("muffleWarning")
      })
    expect_s3_class(got, "rdistnn_partial_result_warning")
    expect_null(got$call)
    expect_identical(got$examples, as.character(1:5))
    expect_match(conditionMessage(got), "1 more", fixed = TRUE)
    expect_match(conditionMessage(got), "Inspect status.", fixed = TRUE)
  }
  withr::local_options(warn = 2)
  expect_error(RdistNN:::.dnn_warn("partial_result", "Results unavailable."),
    "Results unavailable.", fixed = TRUE)
})

test_that("warning families are scoped to a public operation and do not leak", {
  operation <- function() {
    RdistNN:::.dnn_distance_scope()
    RdistNN:::.dnn_warn("partial_result", "Partial.")
    RdistNN:::.dnn_warn("partial_result", "Partial again.")
    RdistNN:::.dnn_warn("weak_localization", "Weak.")
  }
  for (i in 1:2) {
    seen <- character()
    withCallingHandlers(operation(), warning = function(w) {
      seen <<- c(seen, w$code)
      invokeRestart("muffleWarning")
    })
    expect_identical(seen, c("partial_result", "weak_localization"))
  }
  expect_false(exists("..dnn_distance_warning_scope", .GlobalEnv, inherits = FALSE))
})

test_that("condition families distinguish state, extras, dependencies and invariants", {
  mapping <- c(state_invalid = "rdistnn_state_error", extra_unavailable = "rdistnn_extra_unavailable",
    dependency_missing = "rdistnn_dependency_error", internal_invariant = "rdistnn_internal_error")
  for (code in names(mapping)) {
    err <- tryCatch(RdistNN:::.dnn_abort(code, "Problem.", hint = "Action."), error = identity)
    expect_s3_class(err, mapping[[code]])
    expect_match(conditionMessage(err), "Action.", fixed = TRUE)
  }
})

test_that("ordinary valid fits do not construct conditions", {
  local_mocked_bindings(.dnn_abort = function(...) stop("unexpected error construction"),
    .dnn_warn = function(...) stop("unexpected warning construction"), .package = "RdistNN")
  fit <- np_dnn(X = matrix(1:10), x0 = 0, Y = 1:10, s = 3L, se = FALSE)
  expect_true(all(is.finite(fit$mu)))
})
