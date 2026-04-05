# Helper utilities shared by all distance tests
# --------------------------------------------

gen_X <- function(n = 20L, p = 6L, seed = 1L) {
  set.seed(seed)
  matrix(rnorm(n * p), n, p)
}
gen_X0 <- function(m = 5L, p = 6L, seed = 2L) {
  set.seed(seed)
  matrix(rnorm(m * p), m, p)
}
gen_x0 <- function(p = 6L, seed = 3L) {
  set.seed(seed)
  rnorm(p)
}

# Return BOTH an X with a constant column and an inv scale that triggers eps-floor
gen_constcol_and_inv <- function(X, eps = 1e-12) {
  X2 <- X
  p <- ncol(X2)
  if (p >= 1L) X2[, 1L] <- 5
  inv <- rep(1, p)
  if (p >= 1L) inv[1L] <- 1 / eps
  list(X = X2, inv = inv)
}

# Deterministic tie-case matrix for order tests
make_tiecase <- function() {
  X <- rbind(c(1, 0), c(1, 0), c(0, 1), c(2, 0))
  x0 <- c(0, 0) # rows 1 and 2 tie exactly; 1 should appear before 2
  list(X = X, x0 = x0)
}

# Battery for a given spec factory & orderer
#   mk_spec: function(...) -> list(name, setup, eval, eval_batch)
#   orderer: function(X, x0, ...)
#   args:    named list forwarded to mk_spec
#   expect_sq: NULL (no check), TRUE (name should have 'sq'), or FALSE (should not)
run_distance_battery <- function(metric_name, mk_spec, orderer, args = list(),
                                 expect_sq = NULL) {
  test_that(paste0(metric_name, ": spec name is a scalar string (and squared flag if requested)"), {
    sp <- do.call(mk_spec, args)
    nm <- sp$name
    expect_true(is.character(nm) && length(nm) == 1L)
    if (!is.null(expect_sq)) {
      if (isTRUE(expect_sq)) expect_true(grepl("sq", nm)) else expect_false(grepl("sq", nm))
    }
  })

  test_that(paste0(metric_name, ": eval equals column 1 of batched"), {
    X <- gen_X()
    X0 <- gen_X0(p = ncol(X))
    x0 <- X0[1L, ]

    sp <- do.call(mk_spec, args)
    ctx <- if (!is.null(sp$setup)) sp$setup(X) else NULL

    d1 <- sp$eval(X, x0, ctx)
    D2 <- sp$eval_batch(X, X0, ctx)
    expect_equal(as.numeric(d1), as.numeric(D2[, 1L]))
  })

  test_that(paste0(metric_name, ": batched equals loop over eval()"), {
    X <- gen_X()
    X0 <- gen_X0(p = ncol(X))

    sp <- do.call(mk_spec, args)
    ctx <- if (!is.null(sp$setup)) sp$setup(X) else NULL

    D1 <- sp$eval_batch(X, X0, ctx)
    D2 <- sapply(seq_len(nrow(X0)), function(j) sp$eval(X, X0[j, ], ctx))
    expect_equal(D1, as.matrix(D2))
  })

  test_that(paste0(metric_name, ": orderer is deterministic on ties"), {
    tc <- make_tiecase()
    ob <- do.call(orderer, c(list(tc$X, tc$x0), args))
    expect_true(is.list(ob) && all(c("dx", "ord") %in% names(ob)))
    i1 <- which(ob$ord == 1L)
    i2 <- which(ob$ord == 2L)
    expect_lt(i1, i2)
  })

  test_that(paste0(metric_name, ": NA policy errors consistently"), {
    X <- gen_X()
    x0 <- gen_x0(p = ncol(X))
    X[1, 1] <- NA_real_
    sp <- do.call(mk_spec, args)
    ctx <- NULL
    expect_error(sp$eval(X, x0, ctx))
    expect_error(sp$eval_batch(X, matrix(x0, 1, byrow = TRUE), ctx))
  })
}
