test_that("All orderers return consistent structure and lengths", {
  X <- gen_X(n = 17, p = 5)
  x0 <- gen_x0(p = ncol(X))
  inv1 <- rep(1, ncol(X))

  o1 <- order_by_l1(X, x0)
  o2 <- order_by_l2(X, x0, squared = FALSE)
  o2s <- order_by_l2(X, x0, squared = TRUE)
  os1 <- order_by_l1_scaled(X, x0, rule = "custom", scale_vec = inv1)
  os2 <- order_by_l2_scaled(X, x0, rule = "custom", scale_vec = inv1, squared = FALSE)
  os2s <- order_by_l2_scaled(X, x0, rule = "custom", scale_vec = inv1, squared = TRUE)

  for (ob in list(o1, o2, o2s, os1, os2, os2s)) {
    expect_true(is.list(ob) && all(c("dx", "ord") %in% names(ob)))
    expect_length(ob$dx, nrow(X))
    expect_length(ob$ord, nrow(X))
    expect_true(all(sort(ob$ord) == seq_len(nrow(X))))
  }

  # With unit scale, scaled and unscaled orderings must match
  expect_identical(o1$ord, os1$ord)
  expect_identical(o2$ord, os2$ord)
  expect_identical(o2s$ord, os2s$ord)

  # Squared vs unsquared L2 order identical
  expect_identical(o2$ord, o2s$ord)
})
