# Public CATE demo schema consumed by nuisance and score examples.
test_that("CATE demo data includes arm nuisances and the AIPW signal", {
  d <- cate_demo_data(n = 20, seed = 9)$data
  expect_true(all(c("m0_oracle", "m1_oracle", "aipw_oracle") %in% names(d)))
})
