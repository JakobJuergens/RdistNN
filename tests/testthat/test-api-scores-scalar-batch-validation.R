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
