.cate_vcov_data <- function() list(
  X = cbind(a = c(-2, -1, 0, 0, 1, 2, 4), b = c(3, 1, 0, 0, 3, 2, -1)),
  X0 = rbind(c(-.2, .1), c(1.7, 1), c(0, 0)),
  psi = c(2, -1, 3, 7, -2, 4, 1)
)

.cate_vcov_subset <- function(X, X0, psi, s, k = 1L) {
  subsets <- utils::combn(nrow(X), s)
  vapply(seq_len(nrow(X0)), function(j) {
    distance <- rowSums(sweep(X, 2L, X0[j, ], "-")^2)
    mean(vapply(seq_len(ncol(subsets)), function(b) {
      ids <- subsets[, b]
      mean(psi[ids[.tie_reference_order(distance[ids], X[ids, , drop = FALSE], ids)[seq_len(k)]]])
    }, numeric(1L)))
  }, numeric(1L))
}

.cate_vcov_reference <- function(d, s = 3L, k = 2L, refit = FALSE) {
  n <- nrow(d$X)
  L <- t(vapply(seq_len(n), function(i) {
    X <- d$X[-i, , drop = FALSE]
    scale <- if (refit) apply(X, 2L, stats::sd) else rep(1, ncol(X))
    .cate_vcov_subset(sweep(X, 2L, scale, "/"), sweep(d$X0, 2L, scale, "/"), d$psi[-i], s, k)
  }, numeric(nrow(d$X0))))
  C <- sweep(L, 2L, colMeans(L), "-")
  (n - 1) / n * crossprod(C)
}

test_that("CATE vcov returns independent common-deletion target covariance", {
  d <- .cate_vcov_data()
  for (jk in c("fast", "loo")) {
    f <- cate_dnn(X = d$X, X0 = d$X0, pseudo_outcome = d$psi,
      s = 3, k = 2, jackknife = jk, return = character())
    V <- vcov(f)
    ref <- .cate_vcov_reference(d)
    expect_equal(unname(V), ref, tolerance = 1e-12)
    expect_equal(dimnames(V), rep(list(names(coef(f))), 2L))
    expect_equal(unname(diag(V)), f$se^2, tolerance = 1e-12)
    expect_gt(abs(V[1, 2]), .01)
    a <- c(1, -1, 0)
    expect_equal(drop(crossprod(a, V %*% a)), drop(crossprod(a, ref %*% a)), tolerance = 1e-12)
    expect_equal(vcov(unserialize(serialize(f, NULL))), V)
  }
})

test_that("CATE vcov preserves deletion-refitted geometry and frozen candidates", {
  d <- .cate_vcov_data()
  f <- cate_dnn(X = d$X, X0 = d$X0, pseudo_outcome = d$psi,
    s = 3, k = 2, distance = distance_l2_scaled("sd"), jackknife = "loo")
  expect_equal(unname(vcov(f)), .cate_vcov_reference(d, refit = TRUE), tolerance = 1e-12)
  tun <- new_tuning("covariance", eval = function(...) list(s = 3L),
    eval_grid = function(...) data.frame(s = c(2L, 3L, 4L)))
  f <- cate_dnn(X = d$X, X0 = d$X0, pseudo_outcome = d$psi,
    tuning = tun, tuning_mode = "grid", return = character())
  for (i in 1:3) {
    candidate <- candidate_fit(f, i)
    expect_equal(unname(vcov(candidate)), .cate_vcov_reference(d, s = i + 1L, k = 1L), tolerance = 1e-12)
    expect_equal(unname(diag(vcov(candidate))), candidate$se^2, tolerance = 1e-12)
  }
})

test_that("CATE vcov distinguishes zero covariance, duplicates and absent inference", {
  d <- .cate_vcov_data()
  f <- cate_dnn(X = d$X, X0 = d$X0, pseudo_outcome = rep(5, 7), s = 3, return = character())
  expect_equal(unname(vcov(f)), matrix(0, 3, 3))
  f <- cate_dnn(X = d$X, X0 = d$X0[c(1, 1), ], pseudo_outcome = d$psi, s = 3)
  expect_equal(vcov(f)[1, ], unname(vcov(f)[2, ]), ignore_attr = TRUE)
  f <- cate_dnn(X = d$X, x0 = d$X0[1, ], pseudo_outcome = d$psi, s = 3)
  expect_equal(dim(vcov(f)), c(1L, 1L))
  f0 <- cate_dnn(X = d$X, X0 = d$X0, pseudo_outcome = d$psi, s = 3, se = FALSE)
  expect_error(vcov(f0), "se = TRUE")
  f$.model <- NULL
  expect_error(vcov(f), "fitted model state")
})

test_that("CATE covariance follows supplied AIPW signals and asymptotic normalization", {
  d <- cate_demo_data(n = 25, seed = 72)$data
  X0 <- rbind(c(0, 0), c(.2, -.1))
  f <- cate_dnn(data = d, x_cols = c("X1", "X2"), X0 = X0,
    y_col = "Y", t_col = "D", e_col = "e_oracle",
    m0_col = "m0_oracle", m1_col = "m1_oracle", s = 4, return = "psi")
  explicit <- cate_dnn(X = as.matrix(d[c("X1", "X2")]), X0 = X0,
    pseudo_outcome = f$psi, s = 4, return = character())
  expect_equal(vcov(f), vcov(explicit), tolerance = 1e-12)
  d <- .cate_vcov_data()
  for (normalize in c(TRUE, FALSE)) {
    f <- cate_dnn(X = d$X, X0 = d$X0, pseudo_outcome = d$psi,
      method = "asymp", normalize = normalize, s = 3, k = 2)
    L <- t(vapply(1:7, function(i) {
      W <- weights_dnn(X = d$X[-i, , drop = FALSE], X0 = d$X0,
        method = "asymp", normalize = normalize, s = 3, k = 2)$W
      as.numeric(crossprod(d$psi[-i], W))
    }, numeric(3L)))
    C <- sweep(L, 2L, colMeans(L), "-")
    expect_equal(unname(vcov(f)), 6 / 7 * crossprod(C), tolerance = 1e-12)
    expect_equal(unname(diag(vcov(f))), f$se^2, tolerance = 1e-12)
  }
})

test_that("CATE covariance retains registered geometry and validates inference state", {
  d <- .cate_vcov_data()
  name <- "stage1_covariance_metric"
  register_distance(name, distance_l2_scaled("sd"))
  on.exit(rm(list = name, envir = RdistNN:::.dist_registry), add = TRUE)
  f <- cate_dnn(X = d$X, X0 = d$X0, pseudo_outcome = d$psi,
    distance = name, s = 3, jackknife = "loo", return = character())
  V <- vcov(f)
  register_distance(name, distance_l1())
  expect_equal(vcov(f), V)
  f$.model["inference_ctx"] <- NULL
  expect_error(vcov(f), "fitted geometry.*Refit")
  f <- cate_dnn(X = d$X, X0 = d$X0, pseudo_outcome = 2^40 + d$psi, s = 3)
  centered <- cate_dnn(X = d$X, X0 = d$X0, pseudo_outcome = d$psi, s = 3)
  expect_equal(vcov(f), vcov(centered), tolerance = 1e-12)
})

test_that("CATE joint covariance exercises native and forced-R deletions", {
  original <- RdistNN:::.rdistnn_compiled_fun
  native_calls <- 0L
  testthat::local_mocked_bindings(.rdistnn_compiled_fun = function(name) {
    fn <- original(name)
    if (name != "rustat_jackknife_rank_full" || is.null(fn)) return(fn)
    function(...) {
      result <- fn(...)
      native_calls <<- native_calls + 1L
      result
    }
  }, .package = "RdistNN")
  old <- options(RdistNN.use_compiled = TRUE)
  on.exit(options(old), add = TRUE)
  d <- .cate_vcov_data()
  native <- cate_dnn(X = d$X, X0 = d$X0, pseudo_outcome = d$psi, s = 3, k = 2)
  V <- vcov(native)
  expect_gt(native_calls, 0L)
  count <- native_calls
  options(RdistNN.use_compiled = FALSE)
  fallback <- cate_dnn(X = d$X, X0 = d$X0, pseudo_outcome = d$psi, s = 3, k = 2)
  expect_equal(vcov(fallback), V, tolerance = 1e-12)
  expect_identical(native_calls, count)
})
