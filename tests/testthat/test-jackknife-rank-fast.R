# jackknife (rank L-stat): shapes, parity, compiled path & fallback
# -------------------------------------------------------------------
# Tags:
#   [engine]   shape/finite & normalize toggle; s={1,>1}
#   [compiled] C++ parity vs slow reference (matrix phi)
#   [fallback] R fallback (phi(j,t)) parity vs slow reference
#   [ties]     deterministic under duplicate rows
# -------------------------------------------------------------------

test_that("[engine] jackknife_slow_rank returns correct shapes & finite values", {
  withr::local_preserve_seed()
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
  withr::local_preserve_seed()
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
  withr::local_preserve_seed()
  withr::local_options(RdistNN.use_compiled = TRUE)
  set.seed(3)
  native <- .rdistnn_compiled_fun("rustat_jackknife_fast_rank")
  native_calls <- 0L
  local_mocked_bindings(rustat_jackknife_fast_rank = function(...) {
    native_calls <<- native_calls + 1L
    native(...)
  }, .package = "RdistNN")

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
  expect_identical(native_calls, 2L)
})

test_that("[fallback] jackknife_fast_rank (phi as function) matches slow reference", {
  withr::local_preserve_seed()
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
  withr::local_preserve_seed()
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

test_that("[compiled] averages over T, not hits[i]", {
  withr::local_preserve_seed()
  set.seed(123)
  Tn <- 20L
  n_ref <- 50L
  s <- 4L
  q <- 2L
  ranks <- t(replicate(Tn, sample(n_ref, s + 1L, replace = FALSE)))
  a <- c(1, 2, 3, 4)
  PHI <- cbind(rnorm(n_ref), runif(n_ref))

  ref <- jackknife_slow_rank(ranks, a, PHI, normalize = TRUE)
  fast <- jackknife_fast_rank(ranks, a, PHI, normalize = TRUE)

  # if someone divides by hits[i], the max abs diff will be clearly > 1e-8
  expect_lt(max(abs(fast$theta_loo - ref$theta_loo)), 1e-10)
})

test_that("[golden-matrix] small matrix-phi case matches fixed values", {
  # Setup: s=2, T=2, n_ref=5, q=1, normalize=TRUE
  # ranks: row t=1: (2,4 | 5); row t=2: (1,3 | 2)
  ranks <- rbind(
    c(2L, 4L, 5L),
    c(1L, 3L, 2L)
  )
  a <- c(1, 3) # Z = 4
  PHI <- matrix(1:5, nrow = 5, ncol = 1) # phi(j) = j

  out <- jackknife_fast_rank(ranks, a, PHI, normalize = TRUE)

  # Hand-computed:
  # theta0 = average_t [(a1*phi(j1)+a2*phi(j2))/Z] = 3.0
  theta0_exp <- 3.0

  # theta_loo (n_ref x 1):
  # i=1: 2.875; i=2: 3.625; i=3: 2.625; i=4: 3.375; i=5: 3.0
  theta_loo_exp <- c(2.875, 3.625, 2.625, 3.375, 3.000)

  expect_equal(as.numeric(out$theta0), theta0_exp, tolerance = 1e-12)
  expect_equal(drop(out$theta_loo[, 1, drop = TRUE]), theta_loo_exp, tolerance = 1e-12)
})

test_that("[info] output list has correct names and info fields", {
  withr::local_preserve_seed()
  set.seed(10)
  Tn <- 6L; n_ref <- 15L; s <- 3L; q <- 2L
  ranks <- t(replicate(Tn, sample(n_ref, s + 1L, replace = FALSE)))
  a <- runif(s)
  PHI <- cbind(rnorm(n_ref), rnorm(n_ref))

  out <- jackknife_fast_rank(ranks, a, PHI, normalize = TRUE)
  expect_type(out, "list")
  expect_named(out, c("theta0", "theta_loo", "info"), ignore.order = TRUE)
  expect_length(out$theta0, q)
  expect_identical(dim(out$theta_loo), c(n_ref, q))
  expect_true(is.list(out$info))
  expect_identical(out$info$s, s)
  expect_identical(out$info$T, Tn)
  expect_identical(out$info$q, q)
  expect_identical(out$info$normalize, TRUE)

  # normalize = FALSE propagated in info
  out_f <- jackknife_fast_rank(ranks, a, PHI, normalize = FALSE)
  expect_identical(out_f$info$normalize, FALSE)
})

test_that("[errors] K < s+1 raises an error", {
  withr::local_preserve_seed()
  set.seed(11)
  Tn <- 4L; n_ref <- 10L; s <- 3L
  # Only s cols — need s+1; force R fallback via function phi
  ranks_short <- t(replicate(Tn, sample(n_ref, s, replace = FALSE)))
  a <- runif(s)
  PHI <- cbind(rnorm(n_ref))
  phi_fn <- function(j, t) PHI[j, ]
  expect_error(jackknife_fast_rank(ranks_short, a, phi_fn), regexp = "s\\+1")
})

test_that("[errors] ranks diagnostics distinguish NA vs non-finite", {
  withr::local_preserve_seed()
  set.seed(111)
  Tn <- 4L; n_ref <- 10L; s <- 3L
  ranks <- t(replicate(Tn, sample(n_ref, s + 1L, replace = FALSE)))
  a <- runif(s)
  PHI <- cbind(rnorm(n_ref))
  phi_fn <- function(j, t) PHI[j, ]

  ranks_na <- ranks
  ranks_na[1, 1] <- NA_integer_
  expect_error(jackknife_fast_rank(ranks_na, a, phi_fn), regexp = "ranks.*NA")

  ranks_inf <- ranks
  storage.mode(ranks_inf) <- "double"
  ranks_inf[1, 1] <- Inf
  expect_error(jackknife_fast_rank(ranks_inf, a, phi_fn), regexp = "ranks.*non-finite")

  ranks_frac <- ranks
  storage.mode(ranks_frac) <- "double"
  ranks_frac[1, 1] <- 1.5
  expect_error(jackknife_fast_rank(ranks_frac, a, PHI), regexp = "ranks.*integer indices")
})

test_that("[errors] finite-support input validation runs before compiled dispatch", {
  withr::local_preserve_seed()
  set.seed(112)
  Tn <- 4L; n_ref <- 10L; s <- 3L
  ranks <- t(replicate(Tn, sample(n_ref, s + 1L, replace = FALSE)))
  a <- runif(s)
  PHI <- cbind(rnorm(n_ref), rnorm(n_ref))

  expect_error(
    jackknife_fast_rank(ranks, a, PHI, normalize = NA),
    "normalize.*TRUE or FALSE"
  )
  expect_error(
    jackknife_fast_rank(ranks, a, PHI, normalize = "yes"),
    "normalize.*TRUE or FALSE"
  )
  expect_error(
    jackknife_fast_rank(ranks, a, PHI, normalize = c(TRUE, FALSE)),
    "normalize.*TRUE or FALSE"
  )

  a_na <- a
  a_na[1L] <- NA_real_
  expect_error(jackknife_fast_rank(ranks, a_na, PHI), "a.*NA")

  a_inf <- a
  a_inf[1L] <- Inf
  expect_error(jackknife_fast_rank(ranks, a_inf, PHI), "a.*non-finite")

  PHI_bad <- PHI
  PHI_bad[1L, 1L] <- NA_real_
  expect_error(jackknife_fast_rank(ranks, a, PHI_bad), "phi.*NA")

  PHI_bad[1L, 1L] <- Inf
  expect_error(jackknife_fast_rank(ranks, a, PHI_bad), "phi.*non-finite")

  phi_na <- function(j, t) c(NA_real_, 1)
  expect_error(jackknife_fast_rank(ranks, a, phi_na), "phi\\(j, t\\).*NA")

  phi_bad_length <- function(j, t) {
    if (t == 1L && j == ranks[1L, 1L]) c(1, 2) else 1
  }
  expect_error(jackknife_fast_rank(ranks, a, phi_bad_length), "consistent length")
})

test_that("[errors] sum(a)=0 with normalize=TRUE raises an error", {
  withr::local_preserve_seed()
  set.seed(12)
  Tn <- 3L; n_ref <- 8L; s <- 2L
  ranks <- t(replicate(Tn, sample(n_ref, s + 1L, replace = FALSE)))
  a_zero <- c(1, -1)  # sum = 0
  PHI <- cbind(rnorm(n_ref))
  phi_fn <- function(j, t) PHI[j, ]
  expect_error(jackknife_fast_rank(ranks, a_zero, phi_fn, normalize = TRUE),
               regexp = "sum\\(a\\)")
})

test_that("[normalize] normalize=FALSE is a pure rescaling vs normalize=TRUE", {
  withr::local_preserve_seed()
  set.seed(14)
  Tn <- 8L; n_ref <- 20L; s <- 4L
  ranks <- t(replicate(Tn, sample(n_ref, s + 1L, replace = FALSE)))
  a <- runif(s) + 0.1
  PHI <- cbind(rnorm(n_ref), rnorm(n_ref))

  out_t <- jackknife_fast_rank(ranks, a, PHI, normalize = TRUE)
  out_f <- jackknife_fast_rank(ranks, a, PHI, normalize = FALSE)
  Z <- sum(a)
  expect_equal(out_f$theta0,    Z * out_t$theta0,    tolerance = 1e-12)
  expect_equal(out_f$theta_loo, Z * out_t$theta_loo, tolerance = 1e-12)
})

test_that("[golden-fun] small function-phi (t-dependent) matches fixed values", {
  # Setup: s=2, T=2, n_ref=4, q=2, normalize=TRUE
  # ranks: row t=1: (2,3 | 4); row t=2: (1,4 | 2)
  ranks <- rbind(
    c(2L, 3L, 4L),
    c(1L, 4L, 2L)
  )
  a <- c(2, 1) # Z = 3

  # t-dependent phi: phi(j,t) = (j, j + t)
  phi_fun <- function(j, t) c(j, j + t)

  out <- jackknife_fast_rank(ranks, a, phi_fun, normalize = TRUE)

  # Hand-computed:
  # theta0 = (13/6, 11/3)
  theta0_exp <- c(13 / 6, 11 / 3)

  # theta_loo rows (i=1..4):
  # i=1: (17/6, 13/3)
  # i=2: ( 8/3, 25/6)
  # i=3: ( 7/3, 23/6)
  # i=4: (11/6, 10/3)
  theta_loo_exp <- rbind(
    c(17 / 6, 13 / 3),
    c(8 / 3, 25 / 6),
    c(7 / 3, 23 / 6),
    c(11 / 6, 10 / 3)
  )

  expect_equal(as.numeric(out$theta0), theta0_exp, tolerance = 1e-12)
  expect_equal(as.matrix(out$theta_loo), theta_loo_exp, tolerance = 1e-12)
})
