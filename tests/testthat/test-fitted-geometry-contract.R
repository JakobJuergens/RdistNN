test_that("saved core fits retain the resolved metric and fitted context", {
  metric <- "core_audit_saved_metric"
  on.exit(rm(list = metric, envir = .dist_registry), add = TRUE)
  calls <- 0L
  spec <- new_distance(metric,
    setup = function(X) { calls <<- calls + 1L; list(scale = max(X)) },
    eval = function(X, x0, ctx) abs(X[, 1] - x0[1]) / ctx$scale)
  register_distance(metric, spec)
  X <- matrix(1:8)
  tun <- new_tuning("core_audit", eval = function(...) list(s = 3L),
    eval_grid = function(...) data.frame(s = c(3L, 4L)))
  fits <- list(
    np = np_dnn(X = X, Y = 1:8, x0 = 0, distance = metric, tuning = tun,
      tuning_mode = "grid"),
    cate = cate_dnn(X = X, pseudo_outcome = 1:8, x0 = 0, distance = metric,
      tuning = tun, tuning_mode = "grid"),
    ratio = affine_ratio_dnn(X = X, numerator = 1:8, denominator = rep(1, 8),
      x0 = 0, distance = metric, tuning = tun, tuning_mode = "grid"))
  expect_identical(calls, 3L)
  tdnn <- np_tdnn(X = X, Y = 1:8, x0 = 0, distance = metric, s1 = 2, s2 = 4)
  A <- rep(0:1, 4)
  clate <- clate_dnn(X = X, Y = 2 * A + (1:8) / 10, D = A, A = A,
    e = rep(.5, 8), m_y0 = rep(0, 8), m_y1 = rep(2, 8),
    m_d0 = rep(0, 8), m_d1 = rep(1, 8), x0 = 0, distance = metric, s = 3)
  expect_identical(calls, 5L)
  views <- lapply(fits, candidate_fit, 2L)
  before <- lapply(c(fits, views, list(tdnn, clate)), predict, x0 = 0, se.fit = TRUE)
  register_distance(metric, new_distance(metric,
    eval = function(X, x0, ctx = NULL) 9 - abs(X[, 1] - x0[1])))
  after <- lapply(c(fits, views, list(tdnn, clate)), function(fit) {
    predict(unserialize(serialize(fit, NULL)), x0 = 0, se.fit = TRUE)
  })
  expect_equal(after, before, tolerance = 0)
  expect_identical(calls, 5L)
})

test_that("CATE LOO refits geometry and explicit contexts bypass setup", {
  X <- cbind(c(0, 2, 3, 8, 4, 7, 1), c(0, 7, 1, 2, 9, 3, 5))
  targets <- rbind(c(0, 0), c(3, 4))
  signal <- c(1, -2, 3, 4, -1, 6, 2)
  metric <- distance_l2_scaled("sd")
  explicit <- t(vapply(seq_len(nrow(X)), function(i) {
    w <- weights_dnn(X = X[-i, , drop = FALSE], X0 = targets,
      distance = metric, k = 1, s = 3)$W
    as.numeric(crossprod(signal[-i], w))
  }, numeric(nrow(targets))))
  ref <- 6/7 * colSums(sweep(explicit, 2, colMeans(explicit), "-")^2)
  fit <- cate_dnn(X = X, X0 = targets, pseudo_outcome = signal,
    distance = metric, s = 3, jackknife = "loo")
  expect_equal(fit$se^2, ref, tolerance = 1e-12)
  metric$setup <- function(X) stop("setup must not run")
  ctx <- list(inv_s = c(1, .5))
  # Use a custom metric so the context contract is explicit and independent of
  # built-in scaled-distance representation.
  metric <- new_distance("explicit_context", setup = metric$setup,
    eval = function(X, x0, ctx) rowSums(abs(sweep(X, 2, x0, "-")) *
      matrix(ctx$inv_s, nrow(X), ncol(X), byrow = TRUE)))
  for (fun in list(
    function() cate_dnn(X = X, X0 = targets, pseudo_outcome = signal,
      distance = metric, ctx = ctx, s = 3),
    function() affine_ratio_dnn(X = X, X0 = targets, numerator = signal,
      denominator = rep(1, 7), distance = metric, ctx = ctx, s = 3))) {
    fit <- fun()
    expect_equal(predict(unserialize(serialize(fit, NULL)), X0 = targets),
      if (inherits(fit, "cate_result")) fit$tau else fit$theta)
  }
})
