.subset_demo_helpers <- function() {
  list(np = np_demo_subsets, cate = cate_demo_subsets, plqr = plqr_demo_subsets)
}

test_that("demo subset helpers return integer K by s matrices at boundaries", {
  cases <- list(c(4L, 1L), c(1L, 3L), c(1L, 1L), c(4L, 3L), c(7L, 4L))
  for (helper in .subset_demo_helpers()) {
    for (case in cases) {
      K <- case[[1L]]
      s <- case[[2L]]
      draws <- helper(n = 12L, K = K, s = s, seed = 74L)
      expect_identical(dim(draws), c(K, s))
      expect_type(draws, "integer")
      expect_identical(draws, helper(n = 12L, K = K, s = s, seed = 74L))
      expect_true(all(draws >= 1L & draws <= 12L))
      if (!identical(dim(draws), c(K, s))) next
      for (row in seq_len(K)) {
        expect_identical(draws[row, ], sort(draws[row, ]))
        expect_identical(length(unique(draws[row, ])), s)
      }
    }
  }
})

test_that("demo subset helpers preserve row draws and caller RNG consumption", {
  .test_rng_scope({
    for (kind in c("Mersenne-Twister", "L'Ecuyer-CMRG")) {
      RNGkind(kind, normal.kind = "Inversion", sample.kind = "Rejection")
      configured_kind <- RNGkind()
      for (helper in .subset_demo_helpers()) {
        for (s in c(1L, 3L)) {
          set.seed(74L)
          expected <- matrix(NA_integer_, nrow = 4L, ncol = s)
          for (row in seq_len(4L)) expected[row, ] <- sort(sample.int(12L, s))
          after_draws <- .Random.seed

          set.seed(271L)
          before <- .Random.seed
          expect_identical(as.vector(helper(12L, K = 4L, s = s, seed = 74L)),
            as.vector(expected))
          expect_identical(.Random.seed, before)
          expect_identical(RNGkind(), configured_kind)
          rm(".Random.seed", envir = .GlobalEnv)
          expect_identical(as.vector(helper(12L, K = 4L, s = s, seed = 74L)),
            as.vector(expected))
          expect_false(exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))

          set.seed(74L)
          actual <- helper(12L, K = 4L, s = s, seed = NULL)
          expect_identical(as.vector(actual), as.vector(expected))
          expect_identical(.Random.seed, after_draws)
          if (s > 1L) {
            set.seed(74L)
            previous <- t(replicate(4L, sort(sample.int(12L, s))))
            expect_identical(expected, previous)
          }
        }
      }
    }
  })
})

test_that("singleton demo subsets feed batch scores directly", {
  demo <- np_demo_data(n = 12L, seed = 73L)
  score <- score_np_y(s = 1L, x0 = demo$x0,
    x_cols = c("X1", "X2"), y_col = "Y")
  ctx <- score_setup(score, demo$data)
  for (helper in .subset_demo_helpers()) {
    draws <- helper(n = 12L, K = 4L, s = 1L, seed = 74L)
    expect_identical(score_eval_batch(score, draws, demo$data, ctx),
      demo$data$Y[as.vector(draws)])
  }
})

test_that("demo subset helpers reject invalid sampling controls before drawing", {
  .test_rng_scope({
    for (helper in .subset_demo_helpers()) {
      set.seed(271L)
      before <- .Random.seed
      expect_error(helper(n = 3L, K = 2L, s = 4L), "s")
      expect_error(helper(n = 12L, K = 0L), "K")
      expect_error(helper(n = 0L), "n")
      expect_identical(.Random.seed, before)
    }
  })
})
