test_that("scalar rank callbacks preserve singleton and repeated target columns", {
  ctx <- list(scale = 2)
  for (make in list(
    function() new_rank_weight("scalar_rank_shape",
      eval = function(ranks, ctx) ctx$scale * ranks,
      setup = function(...) stop("rank setup is manual")),
    function() as_rank_weight(function(ranks, ctx) ctx$scale * ranks))) {
    spec <- make()
    for (n in c(1L, 3L)) for (m in c(1L, 2L, 3L)) {
      ranks <- matrix(rep(seq_len(n), m), n, m)
      if (m > 1L) ranks[, 2L] <- rev(seq_len(n))
      expected <- matrix(2 * ranks, n, m)
      actual <- rank_weight_eval_batch(spec, ranks, ctx)
      expect_type(actual, "double")
      expect_identical(dim(actual), c(n, m))
      expect_equal(actual, expected, tolerance = 0)
      for (j in seq_len(m)) {
        expect_equal(rank_weight_eval(spec, ranks[, j], ctx), expected[, j], tolerance = 0)
      }
    }
    fit <- weights_matrix(matrix(0), X0 = matrix(c(1, 2, 1)),
      ranking = "l2", weight = spec, ctx = ctx)
    expect_equal(fit$W, matrix(2, 1L, 3L), tolerance = 0)
  }
})

test_that("scalar rank callbacks validate each column before assembly", {
  ranks <- cbind(1:3, 3:1, 1:3)
  for (bad in list(TRUE, c(1, 2), NA_real_, numeric(), "1")) {
    spec <- new_rank_weight("invalid_scalar_rank", eval = local({
      value <- bad
      function(ranks, ctx) if (ranks[1L] == 1L) as.numeric(ranks) else value
    }))
    expect_error(suppressWarnings(rank_weight_eval_batch(spec, ranks)), "rank_weight")
  }
  integer_spec <- new_rank_weight("integer_scalar_rank", eval = function(ranks, ctx) ranks)
  expect_equal(rank_weight_eval_batch(integer_spec, ranks), matrix(as.numeric(ranks), 3L), tolerance = 0)
})
