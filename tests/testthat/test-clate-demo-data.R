test_that("bounded CLATE data preserve structural outcomes and oracle regressions", {
  z <- clate_demo_data(n = 80L, seed = 127L, keep_latent = TRUE)
  d <- z$data
  l <- z$latent
  expect_equal(nrow(d), 80L)
  expect_identical(d$row_id, 1:80)
  expect_true(all(abs(d$Y) <= 1))
  expect_true(all(d$e_oracle >= 3 / 8 & d$e_oracle <= 5 / 8))
  s <- sin(2 * pi * as.matrix(d[paste0("Z", 1:4)]))
  g <- (s[, 1L] + s[, 4L] + d$X1 * s[, 2L] + d$X2 * s[, 3L]) / 8
  tau <- (d$X1 + d$X2) / 8
  a <- 1 / 4 + s[, 1L] / 32
  c <- 1 / 3 + s[, 2L] / 32
  expect_equal(d$m_d0_oracle, a)
  expect_equal(d$m_d1_oracle, a + c)
  expect_equal(d$m_y0_oracle, g + tau * a)
  expect_equal(d$m_y1_oracle, g + tau * (a + c))
  expect_true(all(l$D1 >= l$D0))
  expect_equal(d$D, ifelse(d$A == 1, l$D1, l$D0))
  expect_equal(d$Y, g + tau * d$D + (l$H - .5) / 8 + l$epsilon)
  expect_equal(d$Y, ifelse(d$D == 1, l$Y1, l$Y0))
  expect_equal(z$truth$theta(c(.5, .5)), 1 / 8)
  expect_equal(z$truth$first_stage(rbind(c(.2, .3), c(.7, .8))), rep(1 / 3, 2L))
  expect_equal(z$truth$sigma, 13 / (48 * sqrt(15)))
  expect_equal(z$truth$omega, 39 / (16 * sqrt(15)))
  expect_equal(z$data, clate_demo_data(n = 80L, seed = 127L)$data)
  expect_null(clate_demo_data(n = 5L, seed = 1L)$latent)
})

test_that("heterogeneous bounded data distinguish the complier target from CATE", {
  z <- clate_demo_data(n = 40L, seed = 12L, heterogeneous = TRUE, keep_latent = TRUE)
  x <- rbind(c(.2, .3), c(.6, .7))
  expect_equal(z$truth$theta(x) - z$truth$cate(x), rep(3 / 1024, 2L))
  d <- z$data
  tau <- (d$X1 + d$X2) / 8 + sin(2 * pi * d$Z2) / 16
  expect_equal(z$latent$Y1 - z$latent$Y0, tau)
  expect_true(all(abs(d$Y) <= 1))
  expect_true(is.na(z$truth$sigma))
  expect_true(is.na(z$truth$omega))
  expect_error(z$truth$theta(c(1, 2, 3)), "X")
})

test_that("bounded CLATE generator validates controls and runs with the estimator", {
  for (n in list(0, -1, 1.5, NA_real_, Inf, "10")) {
    expect_error(clate_demo_data(n = n), "n")
  }
  for (seed in list(-1, NA_real_, Inf, c(1, 2), "1")) {
    expect_error(clate_demo_data(seed = seed), "seed")
  }
  expect_error(clate_demo_data(heterogeneous = 1), "heterogeneous")
  expect_error(clate_demo_data(keep_latent = 1), "keep_latent")
  z <- clate_demo_data(n = 36L, seed = 0L)
  f <- clate_dnn(
    data = z$data, x_cols = c("X1", "X2"), y_col = "Y", a_col = "A", d_col = "D",
    e_col = "e_oracle", m_y0_col = "m_y0_oracle", m_y1_col = "m_y1_oracle",
    m_d0_col = "m_d0_oracle", m_d1_col = "m_d1_oracle",
    x0 = z$x0, s = 6L, training = "oracle"
  )
  expect_true(all(is.finite(f$theta)))
  expect_true(all(is.finite(f$se)))
})

test_that("latent interval integration recovers the benchmark means and variance", {
  z <- clate_demo_data(n = 9L, seed = 412L)
  for (i in seq_len(nrow(z$data))) {
    d <- z$data[i, ]
    a <- d$m_d0_oracle
    b <- d$m_d1_oracle
    cuts <- c(0, a, b, 1)
    # Two-point Gauss integration is exact for the quadratic score moments.
    H <- unlist(lapply(seq_len(3L), function(j) {
      mean(cuts[j + 0:1]) + c(-1, 1) * diff(cuts[j + 0:1]) / (2 * sqrt(3))
    }))
    mass <- rep(diff(cuts) / 2, each = 2L)
    h <- rep(H, 2L)
    w <- rep(mass / 2, 2L)
    epsilon <- rep(c(-1 / 8, 1 / 8), each = length(H))
    tau <- (d$X1 + d$X2) / 8
    g <- d$m_y0_oracle - tau * a
    for (A in 0:1) {
      D <- as.numeric(h <= if (A == 0) a else b)
      Y <- g + tau * D + (h - .5) / 8 + epsilon
      expect_equal(sum(w * D), if (A == 0) a else b, tolerance = 1e-14)
      expect_equal(sum(w * Y), d[[if (A == 0) "m_y0_oracle" else "m_y1_oracle"]],
                   tolerance = 1e-14)
    }
    expect_equal(sum(w * ((h - .5) / 8 + epsilon)^2), 13 / 768, tolerance = 1e-14)
  }
  u <- (seq_len(4096L) - .5) / 4096
  e <- .5 + sin(2 * pi * u) / 8
  integrated <- (13 / 768) * mean(1 / (e * (1 - e)))
  expect_equal(integrated, z$truth$sigma, tolerance = 1e-13)
  expect_equal(integrated / (1 / 3)^2, z$truth$omega, tolerance = 1e-13)
})
