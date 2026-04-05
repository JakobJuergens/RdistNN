# tests/testthat/test-jackknife-slow-rank.R
# -------------------------------------------------------------------
# jackknife_slow_rank(): reference delete-one jackknife for fixed ranks
#
# We verify:
#  - [shapes] correct output dims/info fields
#  - [baseline] theta0 equals average of per-target normalized sums
#  - [bf-mat] parity with a brute-force oracle (phi as matrix)
#  - [bf-fun] parity with a brute-force oracle (phi as function)
#  - [tdep]   t-dependent phi path is respected
#  - [s=1]    edge case s = 1 matches hand formula
#  - [s+1only] if i only occurs at s+1, delete-one has no effect
#  - [norm]   normalize=FALSE is a pure rescaling by sum(a)
#  - [errors] invalid inputs raise the documented errors
#  - [det]    determinism
# -------------------------------------------------------------------

test_that("[shapes][baseline] theta0 and dims match the definition", {
  set.seed(1)
  Tn <- 5L
  n_ref <- 7L
  s <- 3L
  q <- 2L
  # Make ranks with s+1 unique neighbors per target
  ranks <- t(replicate(Tn, sample(n_ref, s + 1L, replace = FALSE)))
  a <- c(2, 1, 3)
  PHI <- cbind(1:n_ref, (1:n_ref)^2) # q = 2

  out <- jackknife_slow_rank(ranks, a, PHI, normalize = TRUE)
  expect_type(out, "list")
  expect_named(out, c("theta0", "theta_loo", "info"))
  expect_length(out$theta0, q)
  expect_identical(dim(out$theta_loo), c(n_ref, q))
  expect_true(is.list(out$info))
  expect_identical(out$info$s, s)
  expect_identical(out$info$T, Tn)
  expect_identical(out$info$q, q)
  expect_identical(out$info$normalize, TRUE)

  # Baseline recomputation (definition): theta0 = mean_t (1/Z) sum_k a_k phi(j_k(t))
  Z <- sum(a)
  S_t_divZ <- matrix(0, Tn, q)
  for (k in seq_len(s)) {
    S_t_divZ <- S_t_divZ + a[k] * PHI[ranks[, k], , drop = FALSE]
  }
  S_t_divZ <- S_t_divZ / Z
  theta0_bf <- colMeans(S_t_divZ)
  expect_equal(out$theta0, as.numeric(theta0_bf))
})

# ---------- helpers: brute-force oracle for delete-one with fixed ranks -------
# Works for phi matrix (n_ref x q)
bf_jackknife_matrix <- function(ranks, a, PHI, normalize = TRUE) {
  ranks <- as.matrix(ranks)
  PHI <- as.matrix(PHI)
  Tn <- nrow(ranks)
  s <- length(a)
  q <- ncol(PHI)
  n_ref <- nrow(PHI)
  if (ncol(ranks) < s + 1L) stop("need s+1 columns")
  Z <- if (normalize) sum(a) else 1
  if (normalize && Z == 0) stop("sum(a)==0")
  # baseline
  S_t_divZ <- matrix(0, Tn, q)
  for (k in seq_len(s)) S_t_divZ <- S_t_divZ + a[k] * PHI[ranks[, k], , drop = FALSE]
  S_t_divZ <- S_t_divZ / Z
  theta0 <- colMeans(S_t_divZ)

  theta_loo <- matrix(rep(theta0, each = n_ref), n_ref, q)
  for (i in seq_len(n_ref)) {
    acc <- numeric(q)
    for (t in seq_len(Tn)) {
      nb <- ranks[t, 1:(s + 1L)]
      rpos <- which(nb == i)
      if (length(rpos) && rpos[1L] <= s) {
        r <- rpos[1L]
        # remove nb[r] and pull in nb[s+1] to fill length s
        new_first_s <- c(nb[seq_len(r - 1L)], nb[(r + 1L):(s + 1L)])[seq_len(s)]
      } else {
        new_first_s <- nb[seq_len(s)]
      }
      Snew <- colSums(PHI[new_first_s, , drop = FALSE] * a) / Z
      Sold <- colSums(PHI[nb[seq_len(s)], , drop = FALSE] * a) / Z
      acc <- acc + (Snew - Sold)
    }
    theta_loo[i, ] <- theta0 + acc / Tn
  }
  list(theta0 = as.numeric(theta0), theta_loo = theta_loo)
}

# Works for phi function phi(j, t)
bf_jackknife_function <- function(ranks, a, phi, normalize = TRUE) {
  ranks <- as.matrix(ranks)
  Tn <- nrow(ranks)
  s <- length(a)
  K <- ncol(ranks)
  if (K < s + 1L) stop("need s+1 columns")
  Z <- if (normalize) sum(a) else 1
  if (normalize && Z == 0) stop("sum(a)==0")

  q <- length(phi(ranks[1L, 1L], 1L))
  n_ref <- max(ranks)

  y_tk <- array(0, dim = c(Tn, s + 1L, q))
  for (t in seq_len(Tn)) for (k in seq_len(s + 1L)) y_tk[t, k, ] <- phi(ranks[t, k], t)

  S_t_divZ <- matrix(0, Tn, q)
  for (k in seq_len(s)) S_t_divZ <- S_t_divZ + a[k] * matrix(y_tk[, k, ], Tn, q)
  S_t_divZ <- S_t_divZ / Z
  theta0 <- colMeans(S_t_divZ)

  theta_loo <- matrix(rep(theta0, each = n_ref), n_ref, q)
  for (i in seq_len(n_ref)) {
    acc <- numeric(q)
    for (t in seq_len(Tn)) {
      nb <- ranks[t, 1:(s + 1L)]
      yS <- matrix(y_tk[t, 1:s, ], s, q)
      rpos <- which(nb == i)
      if (length(rpos) && rpos[1L] <= s) {
        r <- rpos[1L]
        y_shift <- matrix(0, s, q)
        if (r > 1L) y_shift[1:(r - 1L), ] <- matrix(y_tk[t, 1:(r - 1L), ], r - 1L, q)
        if (r <= s - 1L) y_shift[r:(s - 1L), ] <- matrix(y_tk[t, (r + 1L):s, ], s - r, q)
        y_shift[s, ] <- y_tk[t, s + 1L, ]
        y_new <- y_shift
      } else {
        y_new <- yS
      }
      acc <- acc + (colSums(y_new * a) - colSums(yS * a)) / Z
    }
    theta_loo[i, ] <- theta0 + acc / Tn
  }
  list(theta0 = as.numeric(theta0), theta_loo = theta_loo)
}

test_that("[bf-mat] matches brute-force oracle for phi as matrix", {
  set.seed(2)
  Tn <- 8L
  n_ref <- 20L
  s <- 4L
  q <- 3L
  ranks <- t(replicate(Tn, sample(n_ref, s + 1L, replace = FALSE)))
  a <- c(1.5, 0.5, 2, 1.25)
  PHI <- cbind(seq_len(n_ref), rnorm(n_ref), runif(n_ref))

  ref <- bf_jackknife_matrix(ranks, a, PHI, normalize = TRUE)
  out <- jackknife_slow_rank(ranks, a, PHI, normalize = TRUE)

  expect_equal(out$theta0, ref$theta0, tolerance = 1e-12)
  expect_equal(out$theta_loo, ref$theta_loo, tolerance = 1e-12)
})

test_that("[bf-fun] matches brute-force oracle for phi as function", {
  set.seed(3)
  Tn <- 6L
  n_ref <- 15L
  s <- 3L
  ranks <- t(replicate(Tn, sample(n_ref, s + 1L, replace = FALSE)))
  a <- c(2, 1, 3)
  # t-dependent phi: q = 2
  phi_fun <- function(j, t) c(j + 0.1 * t, j^2 - 0.2 * t)

  ref <- bf_jackknife_function(ranks, a, phi_fun, normalize = TRUE)
  out <- jackknife_slow_rank(ranks, a, phi_fun, normalize = TRUE)

  expect_equal(out$theta0, ref$theta0, tolerance = 1e-12)
  expect_equal(out$theta_loo, ref$theta_loo, tolerance = 1e-12)
})

test_that("[tdep] t-dependent path actually differs from t-independent", {
  set.seed(4)
  Tn <- 5L
  n_ref <- 12L
  s <- 3L
  ranks <- t(replicate(Tn, sample(n_ref, s + 1L, replace = FALSE)))
  a <- c(1, 2, 1)
  PHI <- cbind(1:n_ref) # q=1
  phi_const <- function(j, t) PHI[j, ]
  phi_tdep <- function(j, t) PHI[j, ] + 0.01 * t

  out_const <- jackknife_slow_rank(ranks, a, phi_const, normalize = TRUE)
  out_tdep <- jackknife_slow_rank(ranks, a, phi_tdep, normalize = TRUE)
  # They should differ unless degenerate ranks/a
  expect_false(isTRUE(all.equal(out_const$theta0, out_tdep$theta0)))
  expect_false(isTRUE(all.equal(out_const$theta_loo, out_tdep$theta_loo)))
})

test_that("[s=1] edge case matches hand formula", {
  set.seed(5)
  Tn <- 7L
  n_ref <- 18L
  s <- 1L
  ranks <- t(replicate(Tn, sample(n_ref, s + 1L, replace = FALSE)))
  a <- 2.5
  PHI <- cbind(runif(n_ref)) # q = 1
  Z <- sum(a)

  out <- jackknife_slow_rank(ranks, a, PHI, normalize = TRUE)
  # hand-check: when i appears as first neighbor at target t, replace y_{t,1} by y_{t,2}
  theta0 <- out$theta0
  theta_loo_hand <- matrix(theta0, n_ref, 1)
  for (i in seq_len(n_ref)) {
    delta <- 0
    for (t in seq_len(Tn)) {
      nb <- ranks[t, 1:2]
      if (nb[1] == i) {
        delta <- delta + (a * PHI[nb[2], ] / Z - a * PHI[nb[1], ] / Z)
      }
    }
    theta_loo_hand[i, ] <- theta0 + delta / Tn
  }
  expect_equal(out$theta_loo, theta_loo_hand, tolerance = 1e-12)
})

test_that("[s+1only] if i only appears at s+1, deletion leaves theta unchanged", {
  set.seed(6)
  Tn <- 5L
  n_ref <- 10L
  s <- 2L
  # Construct ranks so that index 'i_star' never appears in first s positions
  i_star <- 7L
  ranks <- matrix(NA_integer_, Tn, s + 1L)
  for (t in seq_len(Tn)) {
    firsts <- sample(setdiff(1:n_ref, i_star), s, replace = FALSE)
    last <- i_star
    ranks[t, ] <- c(firsts, last)
  }
  a <- c(1, 3)
  PHI <- cbind(rnorm(n_ref), runif(n_ref))
  out <- jackknife_slow_rank(ranks, a, PHI, normalize = TRUE)
  expect_equal(out$theta_loo[i_star, ], out$theta0, tolerance = 0)
})

test_that("[norm] normalize=FALSE is a pure rescaling by sum(a)", {
  set.seed(7)
  Tn <- 6L
  n_ref <- 14L
  s <- 3L
  q <- 2L
  ranks <- t(replicate(Tn, sample(n_ref, s + 1L, replace = FALSE)))
  a <- c(2, 1, 4)
  PHI <- cbind(rnorm(n_ref), rpois(n_ref, 3))
  out_true <- jackknife_slow_rank(ranks, a, PHI, normalize = TRUE)
  out_false <- jackknife_slow_rank(ranks, a, PHI, normalize = FALSE)
  Z <- sum(a)
  expect_equal(out_false$theta0, Z * out_true$theta0, tolerance = 1e-12)
  expect_equal(out_false$theta_loo, Z * out_true$theta_loo, tolerance = 1e-12)
})

test_that("[errors] invalid inputs raise clear errors", {
  set.seed(8)
  Tn <- 3L
  n_ref <- 6L
  s <- 2L
  ranks <- t(replicate(Tn, sample(n_ref, s + 1L, replace = FALSE)))
  a <- c(1, -1)
  PHI <- matrix(rnorm(n_ref), n_ref, 1)

  # K < s+1
  ranks_bad <- ranks[, 1:s, drop = FALSE]
  expect_error(jackknife_slow_rank(ranks_bad, a, PHI), "at least s\\+1")

  # sum(a) == 0 with normalize=TRUE
  a_zero <- c(1, -1)
  expect_error(jackknife_slow_rank(ranks, a_zero, PHI, normalize = TRUE), "sum\\(a\\) == 0")

  # ranks out of PHI range
  ranks_oob <- ranks
  ranks_oob[1, 1] <- n_ref + 1L
  expect_error(jackknife_slow_rank(ranks_oob, a, PHI), "PHI range")
})

test_that("[det] determinism", {
  set.seed(9)
  Tn <- 5L
  n_ref <- 12L
  s <- 3L
  ranks <- t(replicate(Tn, sample(n_ref, s + 1L, replace = FALSE)))
  a <- c(1, 2, 1)
  PHI <- cbind(rnorm(n_ref), runif(n_ref))
  o1 <- jackknife_slow_rank(ranks, a, PHI)
  o2 <- jackknife_slow_rank(ranks, a, PHI)
  expect_identical(o1$theta0, o2$theta0)
  expect_identical(o1$theta_loo, o2$theta_loo)
  expect_identical(o1$info, o2$info)
})

test_that("[sign] shift block has the correct sign (deterministic toy)", {
  Tn <- 1L
  n_ref <- 8L
  s <- 3L
  q <- 1L
  # first s = (2,5,6), s+1 = 4
  ranks <- matrix(c(2, 5, 6, 4), nrow = Tn)
  a <- c(1, 2, 4)
  PHI <- matrix(c(10, 20, 30, 40, 50, 60, 70, 80), nrow = n_ref, ncol = q)

  ref <- bf_jackknife_matrix(ranks, a, PHI, normalize = TRUE)
  out <- jackknife_slow_rank(ranks, a, PHI, normalize = TRUE)

  expect_equal(out$theta0, ref$theta0, tolerance = 1e-12)
  expect_equal(out$theta_loo, ref$theta_loo, tolerance = 1e-12)

  changed <- which(out$theta_loo[, 1] != out$theta0[1])
  expect_setequal(changed, ranks[1, 1:s]) # must be {2,5,6}
})

# -------------------------------------------------------------------
# "Golden value" tests for jackknife_slow_rank()
# These compare against fixed, hand-computed values (no reliance on
# the fast implementation).
# -------------------------------------------------------------------

test_that("[golden-matrix][slow] small matrix-phi case matches fixed values", {
  # Setup: s=2, T=2, n_ref=5, q=1, normalize=TRUE
  # ranks: row t=1: (2,4 | 5); row t=2: (1,3 | 2)
  ranks <- rbind(
    c(2L, 4L, 5L),
    c(1L, 3L, 2L)
  )
  a <- c(1, 3) # Z = 4
  PHI <- matrix(1:5, nrow = 5, ncol = 1) # phi(j) = j

  out <- jackknife_slow_rank(ranks, a, PHI, normalize = TRUE)

  # Hand-computed:
  # theta0 = average_t [(a1*phi(j1)+a2*phi(j2))/Z] = 3.0
  theta0_exp <- 3.0

  # theta_loo (n_ref x 1):
  # i=1: 2.875; i=2: 3.625; i=3: 2.625; i=4: 3.375; i=5: 3.0
  theta_loo_exp <- c(2.875, 3.625, 2.625, 3.375, 3.000)

  expect_equal(as.numeric(out$theta0), theta0_exp, tolerance = 1e-12)
  expect_equal(drop(out$theta_loo[, 1, drop = TRUE]), theta_loo_exp, tolerance = 1e-12)
})

test_that("[golden-fun][slow] small function-phi (t-dependent) matches fixed values", {
  # Setup: s=2, T=2, n_ref=4, q=2, normalize=TRUE
  # ranks: row t=1: (2,3 | 4); row t=2: (1,4 | 2)
  ranks <- rbind(
    c(2L, 3L, 4L),
    c(1L, 4L, 2L)
  )
  a <- c(2, 1) # Z = 3

  # t-dependent phi: phi(j,t) = (j, j + t)
  phi_fun <- function(j, t) c(j, j + t)

  out <- jackknife_slow_rank(ranks, a, phi_fun, normalize = TRUE)

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
