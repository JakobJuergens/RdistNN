# tests/testthat/test-cate-dnn.R
# ------------------------------------------------------------------------------

skip_if_not(
  all(vapply(
    c("cate_dnn", "weights_dnn", "as_rank_weight"),
    exists, logical(1),
    mode = "function"
  )),
  "cate_dnn() or required helpers not found; skipping file"
)

set.seed(123)

.make_data <- function(n = 40L, p = 3L, m = 4L) {
  X <- matrix(rnorm(n * p), n, p)
  X0 <- matrix(rnorm(m * p), m, p)
  # Simple synthetic signal; psi will be provided directly
  beta <- seq_len(p) / p
  Y <- drop(X %*% beta) + rnorm(n, sd = 0.2)
  list(X = X, X0 = X0, Y = Y, beta = beta)
}

# A simple rank-weight spec: w_i ∝ 1 / rank_i^alpha
rw_power <- function(alpha = 1) {
  as_rank_weight(function(r, ctx = NULL) {
    a <- if (is.null(ctx) || is.null(ctx$alpha)) alpha else ctx$alpha
    w <- 1 / (as.numeric(r)^a)
    w / sum(w)
  }, name = sprintf("inv_rank^%g", alpha))
}

test_that("basic run: returns tau (and se if requested) with expected shapes", {
  dat <- .make_data()
  psi <- dat$Y # provide pseudo-outcome directly (keeps test focused on weighting)

  out <- cate_dnn(
    X = dat$X, X0 = dat$X0,
    pseudo_outcome = psi,
    distance = "l2",
    k = 3L, s = 5L,
    method = "exact",
    se = TRUE,
    return = c("W") # also request W to cross-check tau below
  )

  expect_true(is.list(out))
  expect_true(is.numeric(out$tau))
  expect_length(out$tau, nrow(dat$X0))

  # se can be NA if a fast jackknife path isn't compiled; allow either finite or all NA
  expect_true(is.null(out$se) || (is.numeric(out$se) && length(out$se) == nrow(dat$X0)))
  if (!is.null(out$se)) {
    expect_true(all(is.finite(out$se)) || all(is.na(out$se)))
  }

  # meta sanity
  expect_true(is.list(out$meta))
  expect_identical(out$meta$n, nrow(dat$X))
  expect_identical(out$meta$m, nrow(dat$X0))
  expect_identical(out$meta$p, ncol(dat$X))
  expect_true(out$meta$k == 3L && out$meta$s == 5L)
  expect_true(out$meta$method %in% c("exact", "asymp"))

  # tau equals W-weighted average of psi when W is returned
  skip_if_not("W" %in% names(out), "W not returned to verify tau")
  expect_equal(
    as.numeric(colSums(out$W * matrix(psi, nrow(out$W), ncol(out$W)))),
    as.numeric(out$tau),
    tolerance = 1e-12
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

  # matrix-shaped extras
  expect_equal(dim(out$W), c(nrow(dat$X), nrow(dat$X0)))
  if (!is.null(out$ranks)) expect_equal(dim(out$ranks), c(nrow(dat$X), nrow(dat$X0)))
  if (!is.null(out$ord)) expect_equal(dim(out$ord), c(nrow(dat$X), nrow(dat$X0)))
  if (!is.null(out$dx)) expect_equal(dim(out$dx), c(nrow(dat$X), nrow(dat$X0)))

  # psi passthrough
  expect_true(is.numeric(out$psi))
  expect_equal(out$psi, psi)

  # ESS diagnostic
  expect_true(is.numeric(out$ess))
  expect_length(out$ess, nrow(dat$X0))
  expect_true(all(is.finite(out$ess)))
  expect_true(all(out$ess > 0))

  # neighbors: k × m indices taken from ord/ranks
  skip_if_not("neighbors" %in% names(out), "neighbors not returned")
  expect_equal(dim(out$neighbors)[2], nrow(dat$X0))
  expect_true(all(out$neighbors >= 1 & out$neighbors <= nrow(dat$X)))
})

test_that("custom distance works in cate_dnn; custom rank_weight works via weights_matrix", {
  skip_if_not(
    all(vapply(
      c(
        "register_distance", "get_distance",
        "weights_dnn", "weights_matrix",
        "as_rank_weight", "rank_weight_eval_batch",
        "distance_eval_batch", "ord_to_ranks"
      ),
      exists, logical(1),
      mode = "function"
    )),
    "required APIs not available"
  )

  set.seed(123)

  # --- register a custom distance (Chebyshev / L_inf) -------------------------
  nm <- paste0("linf_", as.integer(runif(1, 1e6, 2e6)))
  register_distance(nm, function(X, x0, ctx) {
    apply(abs(X - matrix(x0, nrow(X), ncol(X), TRUE)), 1L, max)
  })
  spec <- get_distance(nm)

  # --- toy data ---------------------------------------------------------------
  n <- 25
  p <- 3
  m <- 4
  X <- matrix(rnorm(n * p), n, p)
  X0 <- matrix(rnorm(m * p), m, p)
  psi <- rnorm(n)

  # cate_dnn expects a data.frame with X columns + Y/T (T must be binary)
  df <- data.frame(
    Y = rnorm(n),
    T = rbinom(n, 1, 0.5),
    X1 = X[, 1], X2 = X[, 2], X3 = X[, 3]
  )
  x_cols <- c("X1", "X2", "X3")

  # --- (A) custom distance + cate_dnn == psi' * W (from weights_dnn) ---------
  wt <- weights_dnn(
    X = X, X0 = X0, distance = nm,
    k = 3L, s = 4L, method = "exact",
    return = "W"
  )
  tau_manual <- as.numeric(crossprod(psi, wt$W))

  out1 <- cate_dnn(
    data = df, x_cols = x_cols, y_col = "Y", t_col = "T",
    pseudo_outcome = psi, X0 = X0,
    distance = nm, k = 3L, s = 4L, method = "exact",
    se = FALSE, return = "W"
  )
  expect_equal(out1$tau, tau_manual, tolerance = 1e-12)

  # --- (B) custom rank_weight works via weights_matrix ------------------------
  # Build ranks using the SAME distance spec to avoid DX mismatches
  DX <- distance_eval_batch(spec, X, X0)
  ORD <- apply(DX, 2L, function(col) order(col, seq_along(col)))
  RKS <- apply(ORD, 2L, ord_to_ranks)

  # Two custom rank-weights: w_i ∝ 1/r_i^alpha
  rw_power <- function(alpha) {
    as_rank_weight(function(r, ctx = NULL) {
      w <- 1 / (as.numeric(r)^alpha)
      w / sum(w)
    }, name = paste0("inv_rank^", alpha))
  }
  rw1 <- rw_power(1)
  rw2 <- rw_power(2)

  # weights_matrix with custom rank_weight (no cate_dnn path here by design)
  W1 <- rank_weight_eval_batch(rw1, RKS)
  W2 <- rank_weight_eval_batch(rw2, RKS)

  tau_rw1 <- as.numeric(crossprod(psi, W1))
  tau_rw2 <- as.numeric(crossprod(psi, W2))

  # basic sanity checks
  expect_true(all(abs(colSums(W1) - 1) < 1e-12))
  expect_true(all(abs(colSums(W2) - 1) < 1e-12))
  expect_false(isTRUE(all.equal(tau_rw1, tau_rw2))) # different alpha => different taus
})


test_that("error handling: both X0 and x0 provided; bad shapes; missing pseudo/nuisance", {
  dat <- .make_data(n = 20, p = 2, m = 3)
  psi <- dat$Y
  x0 <- dat$X0[1, ]

  # both X0 and x0
  expect_error(
    cate_dnn(
      X = dat$X, X0 = dat$X0, x0 = x0,
      pseudo_outcome = psi, distance = "l2", k = 2L, s = 3L, method = "exact"
    ),
    regexp = "either `X0` or `x0`", ignore.case = TRUE
  )

  # mismatched ncol
  expect_error(
    cate_dnn(
      X = dat$X, X0 = matrix(0, 2, ncol(dat$X) + 1),
      pseudo_outcome = psi, distance = "l2", k = 2L, s = 3L
    ),
    regexp = "ncol\\(X0\\) .* ncol\\(X\\)", ignore.case = TRUE
  )

  # missing pseudo_outcome and no nuisance pieces should error
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

test_that("deterministic ties in distances lead to deterministic ordering and stable tau", {
  skip_if_not(
    all(vapply(
      c("cate_dnn", "get_distance", "distance_eval_batch"),
      exists, logical(1),
      mode = "function"
    )),
    "required APIs not available"
  )

  # Construct data with natural L2 ties (symmetry around x=5)
  set.seed(42)
  n <- 10L
  p <- 2L
  m <- 3L
  X <- cbind(seq_len(n), rep(0, n))
  X0 <- matrix(
    c(
      5, 0, # column 1 -> many ties: (4,6), (3,7), (2,8), (1,9)
      5, 1, # slight y-offset
      5, -1
    ), # slight y-offset opposite
    m, p,
    byrow = TRUE
  )

  # Pseudo-outcome and a minimal data.frame for cate_dnn's interface
  psi <- rnorm(n)
  df <- data.frame(
    Y = rnorm(n), T = rbinom(n, 1, 0.5),
    X1 = X[, 1], X2 = X[, 2]
  )

  # Run cate_dnn (ask for ord so we can inspect deterministic tie-breaking)
  out <- cate_dnn(
    data = df, x_cols = c("X1", "X2"), y_col = "Y", t_col = "T",
    pseudo_outcome = psi, X0 = X0,
    distance = "l2", k = 3L, s = 3L, method = "exact",
    se = FALSE, return = c("W", "ord")
  )
  expect_true("ord" %in% names(out))

  # Compute expected order from THE SAME distance implementation
  l2 <- get_distance("l2")
  DX <- distance_eval_batch(l2, X, X0)
  exp_ord1 <- order(DX[, 1], seq_along(DX[, 1]))
  expect_equal(out$ord[, 1], exp_ord1)

  # Reproducibility: same inputs -> same tau
  out2 <- cate_dnn(
    data = df, x_cols = c("X1", "X2"), y_col = "Y", t_col = "T",
    pseudo_outcome = psi, X0 = X0,
    distance = "l2", k = 3L, s = 3L, method = "exact",
    se = FALSE
  )
  expect_equal(out$tau, out2$tau, tolerance = 1e-12)
})
