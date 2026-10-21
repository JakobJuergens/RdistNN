# TDNN rank-weight shapes, signs, dispatch, extras and constant identities.

.tdnn_weights_X <- function(n, p, seed = 1L) {
  .test_rng_scope(matrix(rnorm(n * p), n, p), seed = seed)
}

test_that("[exact] weights_tdnn_exact returns W with correct shape", {
  .test_rng_scope({
    n <- 20L; p <- 2L; m <- 3L
    X  <- .tdnn_weights_X(n, p, seed = 1L)
    X0 <- .make_X0(m, p, seed = 2L)
    out <- weights_tdnn_exact(X, X0 = X0, s1 = 3L, s2 = 8L)
    expect_true(is.list(out))
    expect_true(!is.null(out$W))
    expect_equal(dim(out$W), c(n, m))
  }, seed = 30L)
})

test_that("[exact] weights_tdnn_exact colSums = 1", {
  .test_rng_scope({
    n <- 25L; p <- 2L
    X  <- .tdnn_weights_X(n, p, seed = 3L)
    x0 <- c(0, 0)
    out <- weights_tdnn_exact(X, x0 = x0, s1 = 4L, s2 = 10L)
    expect_equal(colSums(out$W), 1, tolerance = 1e-10)
  }, seed = 30L)
})

test_that("[exact] weights_tdnn_exact has negative weights (bias-correction)", {
  .test_rng_scope({
    n <- 30L; p <- 2L
    X  <- .tdnn_weights_X(n, p, seed = 4L)
    x0 <- c(0.5, -0.5)
    out <- weights_tdnn_exact(X, x0 = x0, s1 = 3L, s2 = 9L)
    expect_true(any(out$W < 0))
  }, seed = 30L)
})

test_that("[exact] weights_tdnn_exact x0 and X0 interfaces agree", {
  .test_rng_scope({
    n <- 20L; p <- 2L
    X  <- .tdnn_weights_X(n, p, seed = 5L)
    x0 <- c(0.1, -0.2)
    X0 <- matrix(x0, nrow = 1L)
    out_x0 <- weights_tdnn_exact(X, x0 = x0, s1 = 3L, s2 = 7L)
    out_X0 <- weights_tdnn_exact(X, X0 = X0, s1 = 3L, s2 = 7L)
    expect_equal(out_x0$W, out_X0$W, tolerance = 0)
  }, seed = 30L)
})

test_that("[exact] weights_tdnn_exact dim defaults to ncol(X)", {
  .test_rng_scope({
    n <- 18L; p <- 3L
    X  <- .tdnn_weights_X(n, p, seed = 6L)
    x0 <- rep(0, p)
    out_default  <- weights_tdnn_exact(X, x0 = x0, s1 = 2L, s2 = 6L)
    out_explicit <- weights_tdnn_exact(X, x0 = x0, s1 = 2L, s2 = 6L, dim = p)
    expect_equal(out_default$W, out_explicit$W, tolerance = 0)
  }, seed = 30L)
})

test_that("[exact] weights_tdnn_exact with NULL X and NULL dim gives a clear error", {
  .test_rng_scope({
    # Before the fix, NULL X with dim = NULL produced integer(0) for dim (via
    # ncol(NULL) = NULL -> as.integer(NULL) = integer(0)), which cascaded into a
    # confusing sprintf warning or deferred failure deep in eval_batch.
    # After the fix, an explicit NULL guard fires immediately.
    expect_error(
      weights_tdnn_exact(NULL, x0 = c(0, 0), s1 = 2L, s2 = 6L),
      "`dim`"
    )
  }, seed = 30L)
})

test_that("[exact] weights_tdnn_exact return='ranks' includes rank matrix", {
  .test_rng_scope({
    n <- 15L; p <- 2L; m <- 2L
    X  <- .tdnn_weights_X(n, p, seed = 7L)
    X0 <- .make_X0(m, p, seed = 8L)
    out <- weights_tdnn_exact(X, X0 = X0, s1 = 2L, s2 = 5L, return = c("W", "ranks"))
    expect_equal(dim(out$ranks), c(n, m))
    expect_true(is.integer(out$ranks))
  }, seed = 30L)
})

test_that("[exact] weights_tdnn_exact return='ord' and 'dx' included when requested", {
  .test_rng_scope({
    n <- 15L; p <- 2L
    X  <- .tdnn_weights_X(n, p, seed = 9L)
    x0 <- c(0, 0)
    out <- weights_tdnn_exact(X, x0 = x0, s1 = 2L, s2 = 6L, return = c("W", "ord", "dx"))
    expect_false(is.null(out$ord))
    expect_false(is.null(out$dx))
    expect_equal(dim(out$ord), c(n, 1L))
    expect_equal(dim(out$dx),  c(n, 1L))
  }, seed = 30L)
})

test_that("[asymp] weights_tdnn_asymp returns W with correct shape", {
  .test_rng_scope({
    n <- 20L; p <- 2L; m <- 3L
    X  <- .tdnn_weights_X(n, p, seed = 10L)
    X0 <- .make_X0(m, p, seed = 11L)
    out <- weights_tdnn_asymp(X, X0 = X0, s1 = 3L, s2 = 8L)
    expect_equal(dim(out$W), c(n, m))
  }, seed = 30L)
})

test_that("[asymp] weights_tdnn_asymp colSums = 1 (normalize=TRUE)", {
  .test_rng_scope({
    n <- 25L; p <- 2L
    X  <- .tdnn_weights_X(n, p, seed = 12L)
    x0 <- c(0, 0)
    out <- weights_tdnn_asymp(X, x0 = x0, s1 = 4L, s2 = 10L)
    expect_equal(colSums(out$W), 1, tolerance = 1e-10)
  }, seed = 30L)
})

test_that("[asymp] weights_tdnn_asymp has negative weights", {
  .test_rng_scope({
    n <- 30L; p <- 2L
    X  <- .tdnn_weights_X(n, p, seed = 13L)
    x0 <- c(0.5, -0.5)
    out <- weights_tdnn_asymp(X, x0 = x0, s1 = 3L, s2 = 9L)
    expect_true(any(out$W < 0))
  }, seed = 30L)
})

test_that("[asymp] weights_tdnn_asymp normalize=FALSE sign structure preserved", {
  .test_rng_scope({
    n <- 20L; p <- 2L
    X  <- .tdnn_weights_X(n, p, seed = 14L)
    x0 <- c(0, 0)
    out_norm <- weights_tdnn_asymp(X, x0 = x0, s1 = 3L, s2 = 7L, normalize = TRUE)
    out_raw  <- weights_tdnn_asymp(X, x0 = x0, s1 = 3L, s2 = 7L, normalize = FALSE)
    expect_equal(sign(out_norm$W), sign(out_raw$W))
    expect_equal(colSums(out_norm$W), 1, tolerance = 1e-10)
  }, seed = 30L)
})

test_that("[asymp] weights_tdnn_asymp x0 and X0 interfaces agree", {
  .test_rng_scope({
    n <- 20L; p <- 2L
    X  <- .tdnn_weights_X(n, p, seed = 15L)
    x0 <- c(0.1, -0.2)
    X0 <- matrix(x0, nrow = 1L)
    out_x0 <- weights_tdnn_asymp(X, x0 = x0, s1 = 3L, s2 = 7L)
    out_X0 <- weights_tdnn_asymp(X, X0 = X0, s1 = 3L, s2 = 7L)
    expect_equal(out_x0$W, out_X0$W, tolerance = 0)
  }, seed = 30L)
})

test_that("[asymp] weights_tdnn_asymp dim defaults to ncol(X)", {
  .test_rng_scope({
    n <- 18L; p <- 3L
    X  <- .tdnn_weights_X(n, p, seed = 16L)
    x0 <- rep(0, p)
    out_default  <- weights_tdnn_asymp(X, x0 = x0, s1 = 2L, s2 = 6L)
    out_explicit <- weights_tdnn_asymp(X, x0 = x0, s1 = 2L, s2 = 6L, dim = p)
    expect_equal(out_default$W, out_explicit$W, tolerance = 0)
  }, seed = 30L)
})

test_that("[asymp] weights_tdnn_asymp with NULL X and NULL dim gives a clear error", {
  .test_rng_scope({
    expect_error(
      weights_tdnn_asymp(NULL, x0 = c(0, 0), s1 = 2L, s2 = 6L),
      "`dim`"
    )
  }, seed = 30L)
})

test_that("[dispatch] weights_tdnn(method='exact') equals weights_tdnn_exact", {
  .test_rng_scope({
    n <- 20L; p <- 2L; m <- 2L
    X  <- .tdnn_weights_X(n, p, seed = 20L)
    X0 <- .make_X0(m, p, seed = 21L)
    out_direct <- weights_tdnn_exact(X, X0 = X0, s1 = 3L, s2 = 8L)
    out_unified <- weights_tdnn(X, X0 = X0, s1 = 3L, s2 = 8L, method = "exact")
    expect_equal(out_unified$W, out_direct$W, tolerance = 0)
  }, seed = 30L)
})

test_that("[dispatch] weights_tdnn(method='asymp') equals weights_tdnn_asymp", {
  .test_rng_scope({
    n <- 20L; p <- 2L; m <- 2L
    X  <- .tdnn_weights_X(n, p, seed = 22L)
    X0 <- .make_X0(m, p, seed = 23L)
    out_direct  <- weights_tdnn_asymp(X, X0 = X0, s1 = 3L, s2 = 8L)
    out_unified <- weights_tdnn(X, X0 = X0, s1 = 3L, s2 = 8L, method = "asymp")
    expect_equal(out_unified$W, out_direct$W, tolerance = 0)
  }, seed = 30L)
})

test_that("[dispatch] weights_tdnn default method is 'exact'", {
  .test_rng_scope({
    n <- 20L; p <- 2L
    X  <- .tdnn_weights_X(n, p, seed = 24L)
    x0 <- c(0, 0)
    out_default <- weights_tdnn(X, x0 = x0, s1 = 3L, s2 = 8L)
    out_exact   <- weights_tdnn(X, x0 = x0, s1 = 3L, s2 = 8L, method = "exact")
    expect_equal(out_default$W, out_exact$W, tolerance = 0)
  }, seed = 30L)
})

test_that("[constant] constant signal: exact weighted mean = constant", {
  .test_rng_scope({
    set.seed(30L)
    n <- 20L; p <- 2L; m <- 3L; c0 <- 7.5
    X     <- .tdnn_weights_X(n, p, seed = 30L)
    X0    <- .make_X0(m, p, seed = 31L)
    out   <- weights_tdnn_exact(X, X0 = X0, s1 = 3L, s2 = 8L)
    psi   <- rep(c0, n)
    ests  <- as.numeric(crossprod(psi, out$W))   # W^T psi (length m)
    expect_equal(ests, rep(c0, m), tolerance = 1e-10)
  }, seed = 30L)
})

test_that("[constant] constant signal: asymp weighted mean = constant", {
  .test_rng_scope({
    n <- 20L; p <- 2L; m <- 3L; c0 <- -2.1
    X     <- .tdnn_weights_X(n, p, seed = 32L)
    X0    <- .make_X0(m, p, seed = 33L)
    out   <- weights_tdnn_asymp(X, X0 = X0, s1 = 3L, s2 = 8L)
    psi   <- rep(c0, n)
    ests  <- as.numeric(crossprod(psi, out$W))
    expect_equal(ests, rep(c0, m), tolerance = 1e-10)
  }, seed = 30L)
})
