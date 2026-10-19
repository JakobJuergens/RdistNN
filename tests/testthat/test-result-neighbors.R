test_that("neighbors have their own internal rank dependency", {
  for (se in c(FALSE, TRUE)) for (jk in c("none", "fast", "loo")) {
    a <- np_dnn(X = matrix(1:8), Y = 1:8, X0 = matrix(c(0, 9)), k = 2, s = 3,
      se = se, jackknife = jk, return = "neighbors")
    expect_equal(a$neighbors, matrix(c(1L, 2L, 8L, 7L), 2))
    expect_null(a$ranks)
    b <- cate_dnn(X = matrix(1:8), pseudo_outcome = 1:8, x0 = 0, s = 3,
      se = se, jackknife = jk, return = "neighbors")
    expect_equal(b$neighbors, matrix(1L))
    c <- np_tdnn(X = matrix(1:8), Y = 1:8, x0 = 0, s1 = 2, s2 = 3,
      se = se, jackknife = jk, return = "neighbors")
    expect_equal(c$neighbors, matrix(1L))
  }
})
