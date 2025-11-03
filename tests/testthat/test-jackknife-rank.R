# jackknife (rank L-stat): shapes, parity, compiled path & fallback
# -------------------------------------------------------------------
# Tags:
#   [engine]   shape/finite & normalize toggle; s={1,>1}
#   [compiled] C++ parity vs slow reference (matrix phi)
#   [fallback] R fallback (phi(j,t)) parity vs slow reference
#   [ties]     deterministic under duplicate rows
# -------------------------------------------------------------------

test_that("[engine] jackknife_slow_rank returns correct shapes & finite values", {
  set.seed(1)
  n <- 60
  p <- 4
  s <- 6
  q <- 3
  X <- matrix(rnorm(n * p), n, p)
  ranks <- knn_ranks(X, k = s + 1L, exclude_self = TRUE)

  PHI <- cbind(rnorm(n), rnorm(n), rnorm(n)) # matrix phi (n x q)
  a <- runif(s)

  out1 <- jackknife_slow_rank(ranks, a, PHI, normalize = TRUE)
  out2 <- jackknife_slow_rank(ranks, a, PHI, normalize = FALSE)

  expect_length(out1$theta0, q)
  expect_true(all(is.finite(out1$theta0)))
  expect_true(all(is.finite(out1$theta_loo)))
  expect_equal(dim(out1$theta_loo), c(n, q))

  expect_length(out2$theta0, q)
  expect_equal(dim(out2$theta_loo), c(n, q))
})

test_that("[engine] boundary case s=1 and scalar phi works", {
  set.seed(2)
  n <- 40
  p <- 3
  s <- 1
  q <- 1
  X <- matrix(rnorm(n * p), n, p)
  ranks <- knn_ranks(X, k = s + 1L, exclude_self = TRUE)
  y <- rnorm(n)
  a <- 2.5 # nonunit to exercise normalize flag

  outT <- jackknife_slow_rank(ranks, a, cbind(y), normalize = TRUE)
  outF <- jackknife_slow_rank(ranks, a, cbind(y), normalize = FALSE)

  expect_length(outT$theta0, q)
  expect_equal(dim(outT$theta_loo), c(n, q))
  expect_true(all(is.finite(outT$theta_loo)))
  expect_true(all(is.finite(outF$theta_loo)))
})

test_that("[compiled] jackknife_fast_rank (matrix phi) matches slow reference", {
  set.seed(3)
  skip_if_not(
    exists("rustat_jackknife_fast_rank", mode = "function"),
    "C++ symbol unavailable; skipping compiled parity test"
  )

  n <- 80
  p <- 5
  s <- 7
  q <- 2
  X <- matrix(rnorm(n * p), n, p)
  ranks <- knn_ranks(X, k = s + 1L, exclude_self = TRUE)
  PHI <- cbind(rnorm(n), rnorm(n))
  a <- runif(s)
  a <- a / mean(a)

  ref <- jackknife_slow_rank(ranks, a, PHI, normalize = TRUE)
  fast <- jackknife_fast_rank(ranks, a, PHI, normalize = TRUE)

  expect_equal(fast$theta0, ref$theta0, tolerance = 1e-12)
  expect_equal(fast$theta_loo, ref$theta_loo, tolerance = 1e-12)

  # also check normalize = FALSE
  ref2 <- jackknife_slow_rank(ranks, a, PHI, normalize = FALSE)
  fast2 <- jackknife_fast_rank(ranks, a, PHI, normalize = FALSE)
  expect_equal(fast2$theta0, ref2$theta0, tolerance = 1e-12)
  expect_equal(fast2$theta_loo, ref2$theta_loo, tolerance = 1e-12)
})

test_that("[fallback] jackknife_fast_rank (phi as function) matches slow reference", {
  set.seed(4)
  n <- 70
  p <- 3
  s <- 5
  X <- matrix(rnorm(n * p), n, p)
  ranks <- knn_ranks(X, k = s + 1L, exclude_self = TRUE)
  y <- rnorm(n)

  # target-dependent, vector-valued phi
  phi <- function(j, t) c(y[j], y[j] * t, 1.0)
  a <- runif(s)

  ref <- jackknife_slow_rank(ranks, a, phi, normalize = TRUE)
  out <- jackknife_fast_rank(ranks, a, phi, normalize = TRUE) # uses R fallback path

  expect_equal(out$theta0, ref$theta0, tolerance = 1e-12)
  expect_equal(out$theta_loo, ref$theta_loo, tolerance = 1e-12)
})

test_that("[ties] duplicate rows yield deterministic outputs", {
  set.seed(5)
  n <- 50
  p <- 2
  s <- 4
  base <- matrix(rnorm((n / 2) * p), n / 2, p)
  X <- rbind(base, base) # induce ties
  ranks1 <- knn_ranks(X, k = s + 1L, exclude_self = TRUE)
  ranks2 <- knn_ranks(X, k = s + 1L, exclude_self = TRUE)
  expect_identical(ranks1, ranks2)

  PHI <- cbind(rnorm(n))
  a <- rep(1, s)

  out1 <- jackknife_fast_rank(ranks1, a, PHI)
  out2 <- jackknife_fast_rank(ranks2, a, PHI)
  expect_equal(out1$theta0, out2$theta0, tolerance = 0)
  expect_equal(out1$theta_loo, out2$theta_loo, tolerance = 0)
})
