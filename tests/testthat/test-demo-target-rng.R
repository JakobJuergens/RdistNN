test_that("invalid demo targets neither reset nor create the caller's RNG state", {
  .test_rng_scope({
    for (demo in list(np_demo_data, cate_demo_data, plqr_demo_data)) {
      for (x0 in list(0, c(0, NA_real_), c(0, Inf))) {
        set.seed(271L)
        before <- .Random.seed
        expect_error(demo(n = 12L, p = 2L, x0 = x0, seed = 7L), "x0")
        expect_identical(.Random.seed, before)
        rm(".Random.seed", envir = .GlobalEnv)
        expect_error(demo(n = 12L, p = 2L, x0 = x0, seed = 7L), "x0")
        expect_false(exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
      }
    }
  })
})
