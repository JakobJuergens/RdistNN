test_that("np_demo_data returns its documented schema reproducibly", {
  first <- np_demo_data(
    n = 24L, p = 3L, design = "sinusoidal", seed = 71L
  )
  second <- np_demo_data(
    n = 24L, p = 3L, design = "sinusoidal", seed = 71L
  )

  expect_named(first, c("data", "x0", "truth", "m_x0"))
  expect_named(first$data, c("X1", "X2", "X3", "Y", "m_oracle"))
  expect_equal(first$data, second$data)
  expect_equal(first$x0, rep(0, 3L))
  expect_equal(first$m_x0, second$m_x0)
})

test_that("np_demo_data zero-noise outcomes equal the oracle for every design", {
  x0 <- c(0.2, -0.1, 0.4)
  for (design in c("linear", "sinusoidal", "bump", "constant")) {
    demo <- np_demo_data(
      n = 30L, p = 3L, design = design, sigma = 0,
      seed = 72L, x0 = x0
    )
    X <- as.matrix(demo$data[paste0("X", 1:3)])

    expect_equal(demo$data$Y, demo$data$m_oracle)
    expect_equal(demo$truth$m_fn(X), demo$data$m_oracle)
    expect_equal(
      demo$m_x0,
      demo$truth$m_fn(matrix(x0, nrow = 1L))
    )
  }
})

test_that("np_demo_data supports diagnostic injection and validates it", {
  na_demo <- np_demo_data(
    n = 12L, p = 2L, seed = 73L,
    inject = "na", idx_inject = 4L
  )
  inf_demo <- np_demo_data(
    n = 12L, p = 2L, seed = 73L,
    inject = "inf", idx_inject = 5L
  )

  expect_true(is.na(na_demo$data$Y[4L]))
  expect_true(is.infinite(inf_demo$data$Y[5L]))
  expect_error(np_demo_data(sigma = -1), "sigma")
  expect_error(np_demo_data(p = 2L, x0 = 0), "x0")
  expect_error(np_demo_data(inject = "bad"), "inject")
  expect_error(
    np_demo_data(n = 5L, inject = "na", idx_inject = 6L),
    "idx_inject"
  )
})

test_that("one-coordinate NP bumps match scalar truth and restore seeded RNG", {
  .test_rng_scope({
    for (seeded in c(TRUE, FALSE)) {
      set.seed(276L)
      before <- .Random.seed
      if (!seeded) rm(".Random.seed", envir = .GlobalEnv)
      demo <- np_demo_data(n = 8L, p = 1L, design = "bump", sigma = 0,
        x0 = .2, seed = 72L)
      expected <- .2 + .8 * exp(-.5 * (demo$data$X1 - .5)^2)
      expect_equal(demo$data$m_oracle, expected, tolerance = 1e-15)
      expect_equal(demo$data$Y, expected, tolerance = 1e-15)
      expect_equal(demo$truth$m_fn(matrix(c(-.1, .2, .5))),
        .2 + .8 * exp(-.5 * (c(-.1, .2, .5) - .5)^2), tolerance = 1e-15)
      expect_equal(demo$m_x0, .2 + .8 * exp(-.5 * (.2 - .5)^2),
        tolerance = 1e-15)
      if (seeded) expect_identical(.Random.seed, before) else
        expect_false(exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
    }
  })
})
