test_that("NP and CATE score factories capture target configuration", {
  dat <- data.frame(X = 1:3, Y = c(2, 4, 6))
  for (kind in c("np", "cate")) {
    x0 <- 0
    cols <- "X"
    score <- if (kind == "np") {
      score_np_y(x0 = x0, x_cols = cols, y_col = "Y")
    } else {
      score_cate_rho(x0 = x0, x_cols = cols, rho_col = "Y", s = 1L)
    }
    x0 <- c(0, 0)
    cols <- "missing"
    ctx <- score_setup(score, dat)
    expect_identical(score_eval_batch(score, matrix(1:3, ncol = 1L), dat, ctx), dat$Y)
  }
})

test_that("grid inference factories capture the prepared target groups", {
  X <- matrix(c(-1, 1, 2, 3, 4), ncol = 1L)
  signal <- c(1, 4, 9, 16, 25)
  ranks <- matrix(1:5, ncol = 1L)
  for (kind in c("np", "cate", "plqr")) {
    groups <- matrix(c(1L, 1L, 2L, 3L, 4L), ncol = 1L)
    callback <- switch(kind,
      np = .np_dnn_grid_se_fn(ranks, signal, X, matrix(0), "l2", "exact", 1L,
        TRUE, NULL, "fast", groups = groups),
      cate = .cate_dnn_grid_se_fn(ranks, signal, X, matrix(0), "l2", "exact", 1L,
        TRUE, NULL, "fast", groups = groups),
      plqr = .plqr_dnn_grid_se_fn(ranks, matrix(signal), X, matrix(0), "l2", "exact", 1L,
        TRUE, NULL, "fast", groups = groups))
    groups <- matrix(1:5, ncol = 1L)
    # Independently enumerate the symmetric nearest-neighbor subset kernel.
    deleted <- vapply(seq_len(nrow(X)), function(i) {
      ids <- setdiff(seq_len(nrow(X)), i)
      subsets <- combn(ids, 3L)
      mean(apply(subsets, 2L, function(S) {
        nearest <- S[abs(X[S, 1L]) == min(abs(X[S, 1L]))]
        mean(signal[nearest])
      }))
    }, numeric(1L))
    expected <- sqrt(4 / 5 * sum((deleted - mean(deleted))^2))
    for (table in list(NULL, c(.75, .25, 0, 0))) {
      out <- callback(data.frame(s = 3L), table)
      expect_equal(out$se, expected, tolerance = 1e-14)
    }
  }
})
