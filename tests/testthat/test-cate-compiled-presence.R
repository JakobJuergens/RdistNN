test_that("compiled symbols present when built", {
  testthat::skip_if_not(has_cate_cpp(), "C++ CATE functions not available")
  expect_true(is.function(rustat_cate_ctx_create))
  expect_true(is.function(rustat_cate_eval_batch))
})
