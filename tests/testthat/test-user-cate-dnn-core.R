# tests/testthat/test-user-cate-dnn-core.R

skip_if_not(
  all(vapply(
    c("cate_dnn", "weights_dnn", "as_rank_weight"),
    exists, logical(1),
    mode = "function"
  )),
  "cate_dnn() or required helpers not found; skipping file"
)

set.seed(123)

test_that("basic run: returns tau (and se if requested) with expected shapes", {
  dat <- .make_data()
  psi <- dat$Y

  out <- cate_dnn(
    X = dat$X, X0 = dat$X0,
    pseudo_outcome = psi,
    distance = "l2",
    k = 3L, s = 5L,
    method = "exact",
    se = TRUE,
    return = c("W")
  )

  expect_true(is.list(out))
  expect_true(is.numeric(out$tau))
  expect_length(out$tau, nrow(dat$X0))
  expect_true(is.null(out$se) || (is.numeric(out$se) && length(out$se) == nrow(dat$X0)))
  if (!is.null(out$se)) {
    expect_true(all(is.finite(out$se)))
  }

  expect_true(is.list(out$meta))
  expect_identical(out$meta$n, nrow(dat$X))
  expect_identical(out$meta$m, nrow(dat$X0))
  expect_identical(out$meta$p, ncol(dat$X))
  expect_true(out$meta$k == 3L && out$meta$s == 5L)
  expect_true(out$meta$method %in% c("exact", "asymp"))

  skip_if_not("W" %in% names(out), "W not returned to verify tau")
  expect_equal(
    as.numeric(colSums(out$W * matrix(psi, nrow(out$W), ncol(out$W)))),
    as.numeric(out$tau),
    tolerance = 1e-12
  )
})

test_that("jackknife modes: none omits SE, loo and fast both return finite SEs", {
  dat <- .make_data(n = 24, p = 2, m = 3)
  psi <- dat$Y

  out_default <- cate_dnn(
    X = dat$X, X0 = dat$X0,
    pseudo_outcome = psi,
    distance = "l2", k = 2L, s = 4L, method = "exact",
    se = TRUE
  )

  out_none <- cate_dnn(
    X = dat$X, X0 = dat$X0,
    pseudo_outcome = psi,
    distance = "l2", k = 2L, s = 4L, method = "exact",
    se = TRUE, jackknife = "none"
  )
  expect_true(is.null(out_none$se))

  out_loo <- cate_dnn(
    X = dat$X, X0 = dat$X0,
    pseudo_outcome = psi,
    distance = "l2", k = 2L, s = 4L, method = "exact",
    se = TRUE, jackknife = "loo"
  )
  expect_true(is.numeric(out_loo$se))
  expect_length(out_loo$se, nrow(dat$X0))
  expect_true(all(is.finite(out_loo$se)))
  expect_true(all(out_loo$se >= 0))

  out_fast <- cate_dnn(
    X = dat$X, X0 = dat$X0,
    pseudo_outcome = psi,
    distance = "l2", k = 2L, s = 4L, method = "exact",
    se = TRUE, jackknife = "fast"
  )
  expect_true(is.numeric(out_fast$se))
  expect_length(out_fast$se, nrow(dat$X0))
  expect_true(all(is.finite(out_fast$se)))
  expect_equal(out_default$se, out_fast$se, tolerance = 1e-12)
  expect_equal(out_fast$se, out_loo$se, tolerance = 1e-12)
})

test_that("delete-1 jackknife requires s < n when se is requested", {
  n <- 10L
  p <- 2L
  X <- matrix(rnorm(n * p), n, p)
  X0 <- matrix(rnorm(3L * p), 3L, p)
  psi <- rnorm(n)

  expect_error(
    cate_dnn(
      X = X, X0 = X0, pseudo_outcome = psi,
      distance = "l2", k = 2L, s = n, method = "exact",
      se = TRUE, jackknife = "loo"
    ),
    regexp = "s < n"
  )
})

test_that("se = FALSE bypasses delete-1 jackknife feasibility check", {
  n <- 10L
  p <- 2L
  X <- matrix(rnorm(n * p), n, p)
  X0 <- matrix(rnorm(2L * p), 2L, p)
  psi <- rnorm(n)

  expect_no_error(
    cate_dnn(
      X = X, X0 = X0, pseudo_outcome = psi,
      distance = "l2", k = 2L, s = n, method = "exact",
      se = FALSE, jackknife = "loo"
    )
  )
})

test_that("constant pseudo_outcome implies zero jackknife SE", {
  n <- 18L
  p <- 2L
  m <- 4L
  X <- matrix(rnorm(n * p), n, p)
  X0 <- matrix(rnorm(m * p), m, p)
  psi <- rep(3.5, n)

  out <- cate_dnn(
    X = X, X0 = X0, pseudo_outcome = psi,
    distance = "l2", k = 2L, s = 4L, method = "exact",
    se = TRUE, jackknife = "loo", return = "W"
  )

  expect_true(all(abs(out$tau - 3.5) < 1e-12))
  expect_true(all(is.finite(out$se)))
  expect_equal(out$se, rep(0, m), tolerance = 1e-10)
})

test_that("data-frame API with pseudo_outcome supports jackknife SEs", {
  dat <- .make_data(n = 26L, p = 3L, m = 4L)
  df <- data.frame(
    Y = rnorm(nrow(dat$X)),
    T = rbinom(nrow(dat$X), 1, 0.5),
    X1 = dat$X[, 1],
    X2 = dat$X[, 2],
    X3 = dat$X[, 3]
  )
  psi <- dat$Y

  out <- cate_dnn(
    data = df, x_cols = c("X1", "X2", "X3"),
    y_col = "Y", t_col = "T",
    pseudo_outcome = psi,
    X0 = dat$X0,
    distance = "l2", k = 2L, s = 4L, method = "asymp",
    se = TRUE, jackknife = "loo"
  )

  expect_true(is.numeric(out$se))
  expect_length(out$se, nrow(dat$X0))
  expect_true(all(is.finite(out$se)))
  expect_true(all(out$se >= 0))
})

test_that("invalid jackknife argument errors clearly", {
  dat <- .make_data(n = 20L, p = 2L, m = 3L)

  expect_error(
    cate_dnn(
      X = dat$X, X0 = dat$X0, pseudo_outcome = dat$Y,
      distance = "l2", k = 2L, s = 3L,
      se = TRUE, jackknife = "bad_option"
    ),
    regexp = "one of"
  )
})

test_that("k and s are validated at the cate_dnn entry point", {
  dat <- .make_data(n = 20L, p = 2L, m = 3L)
  psi <- dat$Y

  # non-positive
  expect_error(
    cate_dnn(X = dat$X, X0 = dat$X0, pseudo_outcome = psi, k = 0L, s = 3L, se = FALSE),
    "positive integer"
  )
  expect_error(
    cate_dnn(X = dat$X, X0 = dat$X0, pseudo_outcome = psi, k = 1L, s = 0L, se = FALSE),
    "positive integer"
  )
  # cross-parameter constraint
  expect_error(
    cate_dnn(X = dat$X, X0 = dat$X0, pseudo_outcome = psi, k = 3L, s = 2L, se = FALSE),
    "k <= s"
  )
  # non-integer scalars (coercion trap)
  expect_error(
    cate_dnn(X = dat$X, X0 = dat$X0, pseudo_outcome = psi, k = 1.5, s = 3L, se = FALSE),
    "integer scalar"
  )
  expect_error(
    cate_dnn(X = dat$X, X0 = dat$X0, pseudo_outcome = psi, k = "foo", s = 3L, se = FALSE),
    "integer scalar"
  )
})

test_that("x0 convenience equals single-row X0 path", {
  dat <- .make_data(n = 20, p = 2, m = 3)
  psi <- dat$Y
  j <- 2L
  x0 <- dat$X0[j, ]

  out_X0 <- cate_dnn(
    X = dat$X, X0 = dat$X0,
    pseudo_outcome = psi,
    distance = "l2", k = 2L, s = 3L,
    method = "exact",
    se = FALSE
  )
  out_x0 <- cate_dnn(
    X = dat$X, x0 = x0,
    pseudo_outcome = psi,
    distance = "l2", k = 2L, s = 3L,
    method = "exact",
    se = FALSE
  )

  expect_equal(out_x0$tau, out_X0$tau[j], tolerance = 1e-12)
})

test_that("return= extras: ranks/ord/dx, psi, ess, neighbors are populated correctly", {
  dat <- .make_data(n = 30, p = 3, m = 5)
  psi <- dat$Y

  out <- cate_dnn(
    X = dat$X, X0 = dat$X0,
    pseudo_outcome = psi,
    distance = "l2",
    k = 4L, s = 5L, method = "exact", se = FALSE,
    return = c("W", "ranks", "ord", "dx", "psi", "ess", "neighbors")
  )

  expect_equal(dim(out$W), c(nrow(dat$X), nrow(dat$X0)))
  if (!is.null(out$ranks)) expect_equal(dim(out$ranks), c(nrow(dat$X), nrow(dat$X0)))
  if (!is.null(out$ord)) expect_equal(dim(out$ord), c(nrow(dat$X), nrow(dat$X0)))
  if (!is.null(out$dx)) expect_equal(dim(out$dx), c(nrow(dat$X), nrow(dat$X0)))
  expect_true(is.numeric(out$psi))
  expect_equal(out$psi, psi)
  expect_true(is.numeric(out$ess))
  expect_length(out$ess, nrow(dat$X0))
  expect_true(all(is.finite(out$ess)))
  expect_true(all(out$ess > 0))

  skip_if_not("neighbors" %in% names(out), "neighbors not returned")
  expect_equal(dim(out$neighbors)[2], nrow(dat$X0))
  expect_true(all(out$neighbors >= 1 & out$neighbors <= nrow(dat$X)))
})

test_that("error handling: both X0 and x0 provided; bad shapes; missing pseudo/nuisance", {
  dat <- .make_data(n = 20, p = 2, m = 3)
  psi <- dat$Y
  x0 <- dat$X0[1, ]

  expect_error(
    cate_dnn(
      X = dat$X, X0 = dat$X0, x0 = x0,
      pseudo_outcome = psi, distance = "l2", k = 2L, s = 3L, method = "exact"
    ),
    regexp = "either `X0` or `x0`", ignore.case = TRUE
  )

  expect_error(
    cate_dnn(
      X = dat$X, X0 = matrix(0, 2, ncol(dat$X) + 1),
      pseudo_outcome = psi, distance = "l2", k = 2L, s = 3L
    ),
    regexp = "ncol\\(X0\\) .* ncol\\(X\\)", ignore.case = TRUE
  )

  expect_error(
    cate_dnn(
      X = dat$X, X0 = dat$X0,
      distance = "l2", k = 2L, s = 3L
    )
  )
})

test_that("metadata fields are stable and informative", {
  dat <- .make_data(n = 22, p = 3, m = 5)
  psi <- dat$Y

  out <- cate_dnn(
    X = dat$X, X0 = dat$X0, pseudo_outcome = psi,
    distance = "l2", k = 3L, s = 4L, method = "asymp",
    se = FALSE
  )

  expect_true(is.list(out$meta))
  expect_identical(out$meta$n, nrow(dat$X))
  expect_identical(out$meta$m, nrow(dat$X0))
  expect_identical(out$meta$p, ncol(dat$X))
  expect_identical(out$meta$k, 3L)
  expect_identical(out$meta$s, 4L)
  expect_identical(out$meta$method, "asymp")
  expect_true(is.character(out$meta$distance))
})
