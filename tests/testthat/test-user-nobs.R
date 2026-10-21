test_that("nobs reports training observations independently of targets and uncertainty", {
  X <- matrix(seq_len(12L), ncol = 1L)
  Y <- sin(X[, 1L])
  X0 <- matrix(c(2, 6, 10), ncol = 1L)
  for (se in c(FALSE, TRUE)) {
    fits <- list(
      np_dnn(X = X, Y = Y, X0 = X0, s = 3L, se = se, return = character()),
      np_tdnn(X = X, Y = Y, X0 = X0, s1 = 2L, s2 = 4L,
        se = se, return = character()),
      cate_dnn(X = X, pseudo_outcome = Y, X0 = X0, s = 3L,
        se = se, return = character())
    )
    for (fit in fits) {
      expect_identical(stats::nobs(fit), 12L)
      expect_identical(stats::nobs(unserialize(serialize(fit, NULL))), 12L)
    }
  }
})
