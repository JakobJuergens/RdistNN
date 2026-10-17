test_that("single-observation rank and weight matrices retain every target", {
  X <- matrix(c(1, 2), 1L)
  X0 <- rbind(c(0, 0), c(1, 3), c(-2, 4))
  for (ranking in list("l2", list(ord = matrix(1L, 1L, 3L)),
                      list(ranks = matrix(1L, 1L, 3L)))) {
    z <- weights_matrix(X, X0 = X0, ranking = ranking,
                         weight = rank_weight_dnn_exact(k = 1L, s = 1L),
                         return = c("W", "ranks", "ord"))
    expect_identical(dim(z$W), c(1L, 3L))
    expect_equal(unname(z$W), matrix(1, 1L, 3L))
    expect_identical(dim(z$ranks), c(1L, 3L))
  }
})
