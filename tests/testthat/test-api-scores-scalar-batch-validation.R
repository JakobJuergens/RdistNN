test_that("parameter-aware scalar callbacks reject invalid batch results before assignment", {
  data <- 1:4
  subs <- matrix(c(1L, 3L, 1L))
  for (route in c("eval_param", "eval")) {
    for (bad in list(c(1, 2), TRUE, NA_real_, Inf, numeric(), "1", list(1))) {
      callback <- local({
        value <- bad
        function(idxs, data, theta, ctx = NULL, rng = NULL) {
          if (idxs == 1L) 7L else value
        }
      })
      args <- list(name = "invalid_param_scalar", order = 1L)
      args[[route]] <- callback
      path <- score_diff_path(do.call(new_score, args), data)
      for (evaluate in list(function() path$value(3L, 0),
                            function() path$value_batch(subs, 0),
                            function() path$grad_batch(subs, 0))) {
        condition <- tryCatch(suppressWarnings(evaluate()), error = identity)
        expect_s3_class(condition, "error")
        if (inherits(condition, "error")) {
          expect_match(conditionMessage(condition), "score\\$.*finite scalar numeric")
        }
      }
    }
  }
})

test_that("parameter-aware scalar batches preserve integer values, contexts and order", {
  data <- c(2L, 4L, 7L)
  subs <- matrix(c(3L, 1L, 3L))
  for (route in c("eval_param", "eval")) {
    setups <- 0L
    args <- list(name = "integer_param_scalar", order = 1L,
      setup = function(data) {
        setups <<- setups + 1L
        list(offset = 5L)
      })
    args[[route]] <- function(idxs, data, theta, ctx = NULL, rng = NULL) {
      as.integer(data[idxs] + theta[1L] + ctx$offset)
    }
    path <- score_diff_path(do.call(new_score, args), data)
    actual <- path$value_batch(subs, 2)
    expect_type(actual, "double")
    expect_identical(actual, c(14, 9, 14))
    expect_equal(path$value(1L, 2), 9, tolerance = 0)
    expect_identical(setups, 1L)
    restored <- unserialize(serialize(path, NULL))
    expect_identical(restored$value_batch(subs[3:1, , drop = FALSE], 2), c(14, 9, 14))
    single <- path$value_batch(matrix(1L), 2)
    expect_identical(single, 9)
  }
})

test_that("differentiable paths bind supplied data when constructed in a loop", {
  datasets <- list(c(2, 3), c(20, 30))
  subs <- matrix(c(2L, 1L, 2L))
  expected_values <- list(c(5.5, 5, 5.5), c(19, 14, 19))
  expected_gradients <- list(cbind(c(3, 2, 3), 4), cbind(c(30, 20, 30), 4))
  theta <- c(.5, 2)
  for (supplied_ctx in c(FALSE, TRUE)) for (method in c("analytic", "fd")) {
    score <- new_score("bound_path_data", 1L,
      setup = if (supplied_ctx) function(data) stop("supplied context must be reused") else NULL,
      eval_param = function(idxs, data, theta, ctx = NULL) {
        theta[1L] * data[idxs] + theta[2L]^2 + if (is.null(ctx)) 0 else ctx$offset
      },
      grad = function(idxs, data, theta) c(data[idxs], 2 * theta[2L]))
    ctx <- if (supplied_ctx) list(offset = 7) else NULL
    paths <- list()
    for (datum in datasets) {
      paths[[length(paths) + 1L]] <- score_diff_path(score, datum, ctx = ctx, method = method)
    }
    datum <- c(200, 300)
    restored <- unserialize(serialize(paths, NULL))
    for (i in 1:2) {
      expected <- expected_values[[i]] + if (supplied_ctx) 7 else 0
      for (path in list(paths[[i]], restored[[i]])) {
        expect_equal(path$value_batch(subs, theta), expected, tolerance = 1e-12)
        expect_equal(path$value(1L, theta), expected[2L], tolerance = 1e-12)
        expect_equal(path$grad_batch(subs, theta), expected_gradients[[i]], tolerance = 1e-8)
      }
    }
  }
})

test_that("differentiable paths evaluate data once after validating options", {
  score <- new_score("path_evaluation_time", 1L,
    eval_param = function(idxs, data, theta) data[idxs] + theta[1L])
  evaluations <- 0L
  path <- score_diff_path(score, {
    evaluations <- evaluations + 1L
    c(2, 3)
  })
  expect_identical(evaluations, 1L)
  expect_equal(path$value(2L, 1), 4)
  expect_identical(evaluations, 1L)
  expect_error(score_diff_path(score, stop("data must stay unevaluated"), fd_step = 0),
    "fd_step.*positive")
})

test_that("analytic batch gradients reject higher-dimensional and non-real outputs", {
  data <- c(2, 4, 7)
  bad_outputs <- list(
    function(K, p) array(1, dim = c(K, p, 1L)),
    function(K, p) array(seq_len(K * p * 2L), dim = c(K, p, 2L)),
    function(K, p) matrix(TRUE, K, p),
    function(K, p) matrix(1 + 1i, K, p),
    function(K, p) matrix(NA_real_, K, p),
    function(K, p) matrix(Inf, K, p))
  for (K in c(1L, 3L)) for (p in c(1L, 2L)) for (output in bad_outputs) {
    score <- new_score("invalid_gradient_shape", 1L,
      eval_param = function(idxs, data, theta) data[idxs] * sum(theta),
      grad_batch = function(subs, data, theta) output(nrow(subs), length(theta)))
    path <- score_diff_path(score, data, method = "analytic")
    subs <- matrix(rep(1L, K))
    expect_error(suppressWarnings(path$grad_batch(subs, seq_len(p))),
      "score\\$grad_batch.*finite matrix")
    expect_error(suppressWarnings(path$grad(1L, seq_len(p))),
      "score\\$grad_batch.*finite matrix")
  }
  for (output in list(function(K) rep(TRUE, K), function(K) rep(1 + 1i, K))) {
    score <- new_score("invalid_gradient_vector", 1L,
      eval_param = function(idxs, data, theta) data[idxs] * theta,
      grad_batch = function(subs, data, theta) output(nrow(subs)))
    path <- score_diff_path(score, data, method = "analytic")
    expect_error(suppressWarnings(path$grad_batch(matrix(1:3), 1)), "finite matrix")
    expect_error(suppressWarnings(path$grad(1L, 1)), "finite matrix")
  }
})

test_that("analytic batch gradients preserve matrix and one-parameter vector outputs", {
  data <- c(2L, 4L, 7L)
  for (K in c(1L, 3L)) for (p in c(1L, 2L)) {
    score <- new_score("valid_gradient_shape", 1L,
      eval_param = function(idxs, data, theta) data[idxs] * sum(seq_along(theta) * theta),
      grad_batch = function(subs, data, theta) outer(data[subs[, 1L]], seq_along(theta), "*"))
    path <- score_diff_path(score, data, method = "analytic")
    subs <- matrix(if (K == 1L) 3L else c(1L, 3L, 1L))
    expected <- if (K == 1L) matrix(c(7, 14)[seq_len(p)], 1L) else {
      matrix(c(2, 7, 2, 4, 14, 4)[seq_len(K * p)], K, p)
    }
    actual <- path$grad_batch(subs, seq_len(p))
    expect_identical(actual, expected)
    expect_identical(path$grad(3L, seq_len(p)), c(7, 14)[seq_len(p)])
  }
  score <- new_score("valid_gradient_vector", 1L,
    eval_param = function(idxs, data, theta) data[idxs] * theta,
    grad_batch = function(subs, data, theta) data[subs[, 1L]])
  path <- score_diff_path(score, data, method = "analytic")
  expect_identical(path$grad_batch(matrix(c(3L, 1L, 3L)), 1), matrix(c(7, 2, 7)))
  expect_identical(path$grad_batch(matrix(1L), 1), matrix(2))
  expect_identical(path$grad(1L, 1), 2)
})

test_that("invalid score return values carry callback identity", {
  bad <- new_score("invalid_value", order = 1L,
    eval = function(idxs, data, ctx, rng) "invalid")
  error <- tryCatch(score_eval(bad, 1L, 1:4), error = identity)
  expect_s3_class(error, "rdistnn_callback_error")
  expect_identical(error$code, "callback_invalid")
  expect_identical(error$argument, "score$eval")
  expect_match(conditionMessage(error), "Correct the callback")
})

test_that("score batches preserve zero-row shapes through value and gradient paths", {
  data <- c(2, 4, 7)
  ordinary <- new_score("empty_ordinary_batch", 2L,
    eval = function(idxs, data, ctx, rng) sum(data[idxs]))
  ordinary_batch <- new_score("empty_explicit_batch", 2L,
    eval_batch = function(subs, data, ctx, rng) {
      rowSums(matrix(data[subs], nrow = nrow(subs), ncol = ncol(subs)))
    })
  for (empty in list(matrix(integer(), 0L, 2L), matrix(numeric(), 0L, 2L))) {
    expect_identical(score_eval_batch(ordinary, empty, data), numeric())
    expect_identical(score_eval_batch(ordinary_batch, empty, data), numeric())
    for (route in c("eval_param", "eval", "eval_param_batch", "eval_batch")) {
      args <- list(name = "empty_parameter_batch", order = 2L)
      args[[route]] <- if (grepl("batch$", route)) {
        function(subs, data, theta) {
          rowSums(matrix(data[subs], nrow = nrow(subs), ncol = ncol(subs))) * sum(theta)
        }
      } else {
        function(idxs, data, theta) sum(data[idxs]) * sum(theta)
      }
      for (derivative in c("grad", "grad_batch", "forward", "central")) {
        local_args <- args
        method <- if (derivative %in% c("grad", "grad_batch")) "analytic" else "fd"
        if (derivative == "grad") {
          local_args$grad <- function(idxs, data, theta) rep(sum(data[idxs]), length(theta))
        } else if (derivative == "grad_batch") {
          local_args$grad_batch <- function(subs, data, theta) {
            outer(rowSums(matrix(data[subs], nrow = nrow(subs), ncol = ncol(subs))),
              rep(1, length(theta)))
          }
        }
        path <- score_diff_path(do.call(new_score, local_args), data, method = method,
          fd_scheme = if (method == "fd") derivative else "central")
        for (theta in list(1, c(1, 2))) {
          expect_identical(path$value_batch(empty, theta), numeric())
          expect_identical(path$grad_batch(empty, theta), matrix(numeric(), 0L, length(theta)))
          # Keep a nonempty control beside every evaluation/gradient route.
          expect_equal(path$value_batch(rbind(c(1L, 3L)), theta), 9 * sum(theta))
          expect_equal(path$grad_batch(rbind(c(1L, 3L)), theta),
            matrix(rep(9, length(theta)), 1L), tolerance = 1e-7)
        }
      }
    }
  }
})

test_that("empty score batches retain input and callback validation", {
  data <- c(2, 4, 7)
  score <- new_score("empty_subset_validation", 2L,
    eval = function(idxs, data, theta) sum(data[idxs]) * sum(theta),
    grad = function(idxs, data, theta) rep(sum(data[idxs]), length(theta)))
  path <- score_diff_path(score, data, method = "analytic")
  for (bad in list(matrix(integer(), 0L, 1L), matrix(logical(), 0L, 2L),
                   matrix(complex(), 0L, 2L), matrix(character(), 0L, 2L))) {
    for (evaluate in list(function() score_eval_batch(score, bad, data),
                          function() path$value_batch(bad, 1),
                          function() path$grad_batch(bad, 1))) {
      error <- expect_error(evaluate(), "columns|numeric or integer matrix",
        class = "rdistnn_input_error")
      expect_identical(error$code, "argument_invalid")
      expect_identical(error$argument, "subs")
    }
  }
  scalar_error <- expect_error(path$value(integer(), 1), "must have 2 elements",
    class = "rdistnn_input_error")
  expect_identical(scalar_error$argument, "idxs")

  empty <- matrix(integer(), 0L, 2L)
  for (bad in list(1, logical(), character())) {
    invalid <- new_score("empty_invalid_callback", 2L,
      eval_param_batch = function(subs, data, theta) bad,
      grad_batch = function(subs, data, theta) matrix(TRUE, 0L, length(theta)))
    invalid_path <- score_diff_path(invalid, data, method = "analytic")
    value_error <- expect_error(invalid_path$value_batch(empty, 1),
      class = "rdistnn_callback_error")
    expect_identical(value_error$code, "callback_invalid")
    expect_identical(value_error$argument, "score$eval_param_batch")
    expect_match(conditionMessage(value_error), "numeric\\(K\\).*finite entries")
    expect_match(value_error$hint, "Correct the callback")
    gradient_error <- expect_error(invalid_path$grad_batch(empty, c(1, 2)),
      class = "rdistnn_callback_error")
    expect_identical(gradient_error$code, "callback_invalid")
    expect_identical(gradient_error$argument, "score$grad_batch")
    expect_match(conditionMessage(gradient_error), "finite matrix")
    invalid_ordinary <- new_score("empty_invalid_ordinary_callback", 2L,
      eval_batch = function(subs, data, ctx, rng) bad)
    ordinary_error <- expect_error(score_eval_batch(invalid_ordinary, empty, data),
      class = "rdistnn_callback_error")
    expect_identical(ordinary_error$code, "callback_invalid")
    expect_identical(ordinary_error$argument, "score$eval_batch")
    expect_match(ordinary_error$hint, "Correct the callback")
  }
})
