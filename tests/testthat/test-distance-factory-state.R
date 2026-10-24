test_that("scaled distance factories capture configuration before first use", {
  old <- options(RdistNN.use_compiled = TRUE)
  on.exit(options(old), add = TRUE)
  dispatch <- .rdistnn_compiled_fun
  native_calls <- character()
  local_mocked_bindings(.rdistnn_compiled_fun = function(name) {
    fun <- dispatch(name)
    if (is.null(fun)) return(NULL)
    function(...) {
      native_calls <<- c(native_calls, name)
      fun(...)
    }
  }, .package = "RdistNN")
  X <- matrix(c(0, 1, 2), ncol = 1L)
  X0 <- matrix(c(0, 1), ncol = 1L)
  expected <- cbind(c(0, .5, 1), c(.5, 0, .5))

  for (compiled in c(FALSE, TRUE)) {
    options(RdistNN.use_compiled = compiled)
    native_calls <- character()
    for (factory in list(distance_l1_scaled, distance_l2_scaled)) {
      for (clamped in c(FALSE, TRUE)) {
        for (batch_first in c(FALSE, TRUE)) {
          rule <- "custom"
          scale <- if (clamped) 0 else 2
          eps <- if (clamped) 2 else 1e-12
          spec <- factory(rule = rule, scale_vec = scale, eps = eps)
          rule <- "range"
          scale <- 8
          eps <- 8
          if (batch_first) {
            expect_equal(distance_eval_batch(spec, X, X0), expected, tolerance = 0)
          } else {
            expect_equal(distance_eval(spec, X, 0), expected[, 1L], tolerance = 0)
          }
          restored <- unserialize(serialize(spec, NULL))
          expect_equal(distance_eval_batch(restored, X, X0), expected, tolerance = 0)
        }
      }
    }
    if (compiled) {
      expect_setequal(unique(native_calls), c("rustat_dist_l1_scaled",
        "rustat_dist_l1_scaled_batch", "rustat_dist_l2_scaled", "rustat_dist_l2_scaled_batch"))
    } else {
      expect_identical(native_calls, character())
    }
  }
})

test_that("fixed custom scales retain configuration without the training matrix", {
  small <- matrix(1:6, ncol = 2L)
  large <- matrix(seq_len(20000L), ncol = 2L)
  for (factory in list(distance_l1_scaled, distance_l2_scaled)) {
    spec <- factory("custom", scale_vec = c(2, 4))
    small_ctx <- spec$setup(small)
    large_ctx <- spec$setup(large)
    expect_identical(large_ctx, small_ctx)
    expect_equal(distance_eval_batch(spec, small, small, large_ctx),
      distance_eval_batch(spec, small, small), tolerance = 0)
  }
})
