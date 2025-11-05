test_that("L1 matches basic expectations", {
  X <- gen_X()
  x0 <- gen_x0(p = ncol(X))
  d <- dist_l1(X, x0)
  # non-negativity and zero iff x0 is a row
  expect_true(all(d >= 0))
  expect_equal(min(dist_l1(rbind(x0, X), x0)), 0)
})

run_distance_battery(
  metric_name = "l1",
  mk_spec = distance_l1,
  orderer = order_by_l1,
  args = list(),
)
