test_that("convenience weights validate integers before coercion", {
  X <- matrix(1:8)
  for (fun in list(weights_dnn, weights_dnn_exact, weights_dnn_asymp)) {
    expect_error(fun(X = X, x0 = 0, k = 1.5, s = 3), "k")
    expect_error(fun(X = X, x0 = 0, k = 1, s = 3.5), "s")
    for (value in c(NA_real_, Inf, -.5, 2^32)) {
      expect_error(fun(X = X, x0 = 0, k = value, s = 3), "k")
    }
  }
  for (fun in list(weights_tdnn, weights_tdnn_exact, weights_tdnn_asymp)) {
    for (arg in c("s1", "s2", "dim")) {
      args <- list(X = X, x0 = 0, s1 = 2, s2 = 4, dim = 1)
      args[[arg]] <- args[[arg]] + .5
      expect_error(do.call(fun, args), arg)
    }
  }
})

test_that("normalize FALSE preserves the raw asymptotic rank table", {
  raw <- dnn_weights_asymp_full_from_ranks(1:80, 1, 2, FALSE)$w_rank
  spec <- rank_weight_dnn_asymp(1, 2, normalize = FALSE)
  expect_identical(rank_weight_eval(spec, 1:80), as.numeric(raw))
  expect_identical(as.numeric(weights_dnn_asymp(X = matrix(1:80), x0 = 0,
    k = 1, s = 2, normalize = FALSE)$W), as.numeric(raw))
})
