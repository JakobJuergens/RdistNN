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
  for (value in list(TRUE, list("W"), NA_character_)) {
    err <- tryCatch(RdistNN:::.warn_unknown_return(value, "W"), error = identity)
    expect_s3_class(err, "rdistnn_input_error")
    expect_identical(err$argument, "return")
  }
  expect_no_condition(RdistNN:::.warn_unknown_return(NULL, "W"))
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
    dependency_missing = "rdistnn_dependency_error", internal_invariant = "rdistnn_internal_error",
    argument_invalid = "rdistnn_input_error", engine_unavailable = "rdistnn_engine_error",
    numerical_range = "rdistnn_numerical_error", callback_invalid = "rdistnn_callback_error",
    candidate_failure = "rdistnn_candidate_failure")
  for (code in names(mapping)) {
    err <- tryCatch(RdistNN:::.dnn_abort(code, "Problem.", hint = "Action."), error = identity)
    expect_s3_class(err, mapping[[code]])
    expect_match(conditionMessage(err), "Action.", fixed = TRUE)
  }
  parent <- structure(list(message = strrep("problem ", 1000),
    call = quote(callback(private_data)), payload = raw(100000)),
    class = c(strrep("x", 200), "rdistnn_method_failure", "error", "condition"))
  err <- tryCatch(RdistNN:::.dnn_abort("callback_invalid", "Invalid result.", parent = parent), error = identity)
  expect_s3_class(err, "rdistnn_callback_error")
  expect_true(all(nchar(err$parent$class) <= 80L))
  expect_lte(nchar(err$parent$message), 240L)
  expect_named(err$parent, c("message", "class"))
  expect_lt(as.numeric(object.size(err)), 8192)
})

test_that("ordinary valid fits do not construct conditions", {
  local_mocked_bindings(.dnn_abort = function(...) stop("unexpected error construction"),
    .dnn_warn = function(...) stop("unexpected warning construction"),
    .dnn_message_examples = function(...) stop("unexpected message formatting"), .package = "RdistNN")
  fit <- np_dnn(X = matrix(1:10), x0 = 0, Y = 1:10, s = 3L, se = FALSE)
  expect_true(all(is.finite(fit$mu)))
})

test_that("choice validation preserves base matching and gives typed public facts", {
  select <- function(method = c("exact", "asymp")) RdistNN:::.dnn_match_arg(method)
  expect_identical(select(), "exact")
  expect_identical(select(NULL), "exact")
  expect_identical(select("asy"), "asymp")
  for (value in list("excat", NA_character_, character(), 1, c("exact", "bad"))) {
    err <- tryCatch(select(value), error = identity)
    expect_s3_class(err, "rdistnn_input_error")
    expect_identical(err$code, "argument_invalid")
    expect_identical(err$argument, "method")
    expect_match(err$expected, "exact", fixed = TRUE)
    expect_match(err$expected, "asymp", fixed = TRUE)
  }
  expect_identical(RdistNN:::.dnn_match_arg(c("asy", "bad", "exa"),
    c("exact", "asymp"), several.ok = TRUE), c("asymp", "exact"))
})

test_that("owned public conditions identify the operation without retaining its call", {
  err <- tryCatch(np_dnn(X = matrix(1:8), Y = 1:8, x0 = 0, s = 3, k = 0), error = identity)
  expect_s3_class(err, "rdistnn_input_error")
  expect_identical(err$operation, "np_dnn")
  expect_null(conditionCall(err))
  expect_match(conditionMessage(err), "np_dnn():", fixed = TRUE)
  expect_false(any(vapply(err, is.environment, logical(1L))))
  args <- list(X = matrix(1:8), Y = 1:8, x0 = 0, s = 3, k = 0)
  .alias <- np_dnn
  for (operation in list(function() do.call(np_dnn, args),
      function() do.call("np_dnn", args), function() do.call(.alias, args))) {
    err <- tryCatch(operation(), error = identity)
    expect_identical(err$operation, "np_dnn")
    expect_null(conditionCall(err))
  }
  fit <- np_dnn(X = matrix(1:8), Y = 1:8, x0 = 0, s = 3, se = FALSE)
  err <- tryCatch(predict(fit, x0 = "bad"), error = identity)
  expect_identical(err$operation, "predict.rdistnn_np_result")
  expect_null(conditionCall(err))
})

test_that("warning locations identify targets rather than bare numbers", {
  got <- NULL
  withCallingHandlers(RdistNN:::.dnn_warn("partial_result", "Unavailable.",
    affected = 1L, examples = 2L, location = "target"), warning = function(w) {
      got <<- w
      invokeRestart("muffleWarning")
    })
  expect_s3_class(got, "rdistnn_partial_result_warning")
  expect_identical(got$location, "target")
  expect_match(conditionMessage(got), "target 2", fixed = TRUE)
})

test_that("localization facts remain accessible when warnings and extras are absent", {
  X <- matrix(c(-1, 1, 2, 3, 4))
  fit <- suppressWarnings(np_dnn(X = X, X0 = matrix(c(0, 10)), Y = 1:5,
    s = 3, se = FALSE, return = character()))
  facts <- fit$meta$target_notices
  expect_s3_class(facts, "data.frame")
  expect_identical(facts$target_id, 1:2)
  expect_identical(facts$distance_ties, c(TRUE, FALSE))
  expect_identical(facts$weak_localization, c(FALSE, FALSE))
  expect_identical(fit$meta$n_tied_targets, 1L)
  expect_identical(fit$meta$n_weak_targets, 0L)
  expect_null(fit$W)
  restored <- unserialize(serialize(fit, NULL))
  expect_equal(unname(coef(restored)), fit$mu)
  broken <- restored
  broken$meta$target_notices$distance_ties[] <- FALSE
  expect_error(coef(broken), class = "rdistnn_state_error")
  broken <- restored
  broken$meta$n_tied_targets <- 0L
  expect_error(coef(broken), class = "rdistnn_state_error")
  broken <- restored
  broken$meta$target_notices <- NULL
  expect_error(coef(broken), class = "rdistnn_state_error")
})

test_that("candidate notice facts reflect candidate weights without regrouping", {
  fit <- suppressWarnings(np_dnn(X = matrix(1:12), x0 = 0, Y = 1:12,
    tuning = new_tuning("notice_grid", eval = function(...) list(s = 3L),
      eval_grid = function(...) data.frame(s = c(1L, 3L))), tuning_mode = "grid",
    se = FALSE, return = character()))
  candidate <- suppressWarnings(candidate_fit(fit, 1L))
  rows <- fit$tuning_grid$estimates
  chosen <- fit$tuning_grid$candidates$candidate_id[[1L]]
  expect_identical(candidate$meta$n_weak_targets,
    as.integer(candidate$meta$s == 1L))
  expect_identical(candidate$meta$n_tied_targets, 0L)
  expect_equal(unname(coef(candidate)), rows$estimate[rows$candidate_id == chosen])
})
