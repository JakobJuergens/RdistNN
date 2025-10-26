# tests/testthat/test-cate-score.R
# -------------------------------------------------------------------
# CATE score test suite
# Purpose: All tests specific to the CATE score (factory, setup, eval,
#          compiled vs R fallback, algebraic properties).
#
# Layout (by section tag):
#   [smoke]   Symbols / environment
#   [engine]  Engine parity & fallbacks (compiled vs R)
#   [ties]    Determinism & tie-breaking
#   [factory] Factory-level validation (rho XOR rho_col)
#   [setup]   Setup/order, finite outputs, basic sanity
#   [cohere]  Single-subset vs batch coherence
#   [perm]    Kernel permutation invariance
#   [rho]     Non-finite rho handling
#   [symbols] Compiled symbol presence
#   [convex]  Convex-hull localization property
#   [inputs]  Early/late input validation & rho_col equivalence
#   [equiv]   dnnk with k=1 equals dnn1
#   [bandw]   Kernel auto-bandwidth median rule (even s)
#   [compiled]Compiled path: finite values & subset-size respect
#
# Conventions:
#   * XC refers to the X-column names used in demo data.
#   * We prefer explicit seeds near stochastic calls.
# -------------------------------------------------------------------

testthat::local_edition(3)
XC <- c("X1", "X2")

# [smoke] Ensure core entry points exist where needed -----------------

test_that("[smoke] u_statistic() exists for integration tests", {
  skip_if_not(exists("u_statistic", mode = "function"),
              "u_statistic() not found; integration tests will skip.")
  succeed()
})

# [engine] Engine parity & fallback ----------------------------------

test_that("[engine] compiled and R fallback produce identical fit", {
  skip_if_not(exists("u_statistic", mode = "function"),
              "u_statistic() not found")
  skip_if_not(has_cate_cpp(), "C++ CATE functions not available")

  demo <- cate_demo_data(n = 180, p = 2, seed = 77)
  df <- demo$data; x0 <- demo$x0; rho <- df$rho_oracle

  sc <- score_cate_rho(x0 = x0, x_cols = XC, rho = rho, localizer = "dnnk", k = 5)

  # Compiled path
  u_cpp <- u_statistic(data = df, score = sc, sampling = list(mode = "fixedN", N = 25000))
  fit_cpp <- fit(u_cpp)$value

  # Force R fallback (clear compiled ptr in ctx)
  ctx <- score_setup(sc, df); ctx$ptr_ctx <- NULL
  sc_r <- sc
  sc_r$setup <- function(data, ...) { ctx }
  u_r <- u_statistic(data = df, score = sc_r, sampling = list(mode = "fixedN", N = 25000))
  fit_r <- fit(u_r)$value

  expect_equal(fit_cpp, fit_r, tolerance = 1e-8)
})

# [ties] Determinism & tie-breaking ----------------------------------

test_that("[ties] exact ties break deterministically (lower index first)", {
  n <- 10
  df <- data.frame(
    X1 = c(0, 0, rep(10, n-2)),
    X2 = c(1, -1, rep(10, n-2)),
    Y = 0, D = 0
  )
  rho <- seq_len(n)   # reveals which indices are chosen
  x0  <- c(0, 0)

  sc  <- score_cate_rho(x0, XC, rho = rho, localizer = "dnnk", k = 2)
  ctx <- score_setup(sc, df)

  subs <- rbind(c(1,2,3,4,5,6), c(2,1,3,4,5,6)) # both contain the tied pair
  v_cpp <- score_eval_batch(sc, subs, df, ctx)   # compiled if available
  ctx$ptr_ctx <- NULL
  v_r   <- score_eval_batch(sc, subs, df, ctx)   # R fallback

  expect_equal(v_cpp[1], 1.5, tolerance = 1e-12) # (1 + 2)/2
  expect_equal(v_cpp, v_r, tolerance = 1e-12)
})

# [factory] Factory-level validation ---------------------------------

test_that("[factory] score_cate_rho enforces XOR: provide exactly one of rho or rho_col", {
  demo <- cate_demo_data(n = 50, p = 2, seed = 123)
  df <- demo$data
  x0 <- demo$x0
  rho <- df$rho_oracle

  expect_error(
    score_cate_rho(x0 = x0, x_cols = XC),
    regexp = "exactly one|mutually exclusive", ignore.case = TRUE
  )
  expect_error(
    score_cate_rho(x0 = x0, x_cols = XC, rho_col = "rho", rho = rho),
    regexp = "exactly one|mutually exclusive", ignore.case = TRUE
  )

  expect_true(is_score(score_cate_rho(x0 = x0, x_cols = XC, rho = rho)))
  expect_true(is_score(score_cate_rho(x0 = x0, x_cols = XC, rho_col = "rho")))
})

# [setup] Setup/order, finite outputs, sanity ------------------------

test_that("[setup] setup builds ctx; order is integer; batch outputs are finite", {
  demo <- cate_demo_data(n = 80, p = 2, seed = 1)
  df <- demo$data; x0 <- demo$x0; rho <- df$rho_oracle

  for (loc in c("dnn1","dnnk","kernel")) {
    sc <- score_cate_rho(x0 = x0, x_cols = XC, rho = rho, localizer = loc, k = 5)
    ctx <- score_setup(sc, df)
    s   <- score_order(sc, nrow(df), df, ctx)
    expect_type(s, "integer")
    expect_true(s >= 1L)

    set.seed(42)
    subs <- cate_demo_subsets(n = nrow(df), K = 25, s = 8)
    vals <- score_eval_batch(sc, subs, df, ctx)

    expect_length(vals, nrow(subs))
    expect_true(all(is.finite(vals)))
  }
})

# [cohere] Single-subset vs batch coherence --------------------------

test_that("[cohere] single-subset eval equals corresponding batch value", {
  demo <- cate_demo_data(n = 60, p = 2, seed = 2)
  df <- demo$data; x0 <- demo$x0; rho <- df$rho_oracle

  sc  <- score_cate_rho(x0 = x0, x_cols = XC, rho = rho, localizer = "dnnk", k = 3)
  ctx <- score_setup(sc, df)
  ctx$ptr_ctx <- NULL  # force R path

  set.seed(7)
  subs <- cate_demo_subsets(n = nrow(df), K = 10, s = 7)
  vals_b <- score_eval_batch(sc, subs, df, ctx)

  vals_s <- vapply(
    seq_len(nrow(subs)),
    function(i) score_eval(sc, subs[i, ], df, ctx),
    numeric(1)
  )
  expect_equal(vals_b, vals_s, tolerance = 1e-12)
})

# [perm] Kernel permutation invariance -------------------------------

test_that("[perm] kernel localizer is permutation-invariant within subsets", {
  demo <- cate_demo_data(n = 75, p = 2, seed = 3)
  df <- demo$data; x0 <- demo$x0; rho <- df$rho_oracle

  sc  <- score_cate_rho(x0, XC, rho = rho, localizer = "kernel", h = NULL)
  ctx <- score_setup(sc, df)
  ctx$ptr_ctx <- NULL  # force R path

  set.seed(11)
  subs <- cate_demo_subsets(n = nrow(df), K = 8, s = 6)
  for (i in seq_len(nrow(subs))) {
    idxs <- subs[i, ]
    perm <- sample(idxs)
    v1 <- score_eval(sc, idxs, df, ctx)
    v2 <- score_eval(sc, perm, df, ctx)
    expect_equal(v1, v2, tolerance = 1e-12)
  }
})

# [rho] Non-finite rho handling --------------------------------------

test_that("[rho] rho vector must be finite", {
  demo <- cate_demo_data(n = 20, p = 2, seed = 1)
  df <- demo$data
  bad_rho <- demo$data$rho_oracle
  bad_rho[3] <- NA_real_

  sc <- score_cate_rho(x0 = demo$x0, x_cols = XC, rho = bad_rho)
  ctx <- score_setup(sc, df); ctx$ptr_ctx <- NULL

  subs <- matrix(c(1,2,3,4,5,6,7,8), nrow = 1)
  expect_error(score_eval_batch(sc, subs, df, ctx), "non-finite", ignore.case = TRUE)
})

# [symbols] Compiled symbol presence ---------------------------------

test_that("[symbols] compiled symbols exported when C++ is available", {
  testthat::skip_if_not(has_cate_cpp(), "C++ CATE functions not available")
  expect_true(is.function(rustat_cate_ctx_create))
  expect_true(is.function(rustat_cate_eval_batch))
})

# [convex] Convex-hull localization property -------------------------

test_that("[convex] localized value lies within convex hull of subset rhos", {
  demo <- cate_demo_data(n = 60, p = 2, seed = 3)
  df <- demo$data; x0 <- demo$x0; rho <- df$rho_oracle

  for (loc in c("dnn1","dnnk","kernel")) {
    sc  <- score_cate_rho(x0, XC, rho = rho, localizer = loc, k = 3, h = NULL)
    ctx <- score_setup(sc, df); ctx$ptr_ctx <- NULL  # test R path

    set.seed(10)
    subs <- cate_demo_subsets(nrow(df), K = 25, s = 6)
    vals <- score_eval_batch(sc, subs, df, ctx)

    mins <- apply(subs, 1L, function(id) min(rho[id]))
    maxs <- apply(subs, 1L, function(id) max(rho[id]))
    expect_true(all(vals >= mins - 1e-12))
    expect_true(all(vals <= maxs + 1e-12))
  }
})

# [inputs] Early/late validation & rho_col equivalence ---------------

test_that("[inputs] factory/setup validate dimensions; rho_col equals rho", {
  demo <- cate_demo_data(n = 30, p = 3, seed = 1)
  df <- demo$data
  rho <- df$rho_oracle

  expect_error(score_cate_rho(x0 = rep(0, 3), x_cols = 1:3),
               "exactly one", ignore.case = TRUE)
  expect_error(score_cate_rho(x0 = rep(0, 3), x_cols = 1:3, rho = rho, rho_col = "rho"),
               "exactly one", ignore.case = TRUE)

  expect_error(score_cate_rho(x0 = 0, x_cols = 1:3, rho = rho) |> score_setup(df),
               regexp = "length|ncol", ignore.case = TRUE)
  expect_error(score_cate_rho(x0 = rep(0, 3), x_cols = 1:3, rho = rho[-1]) |> score_setup(df),
               regexp = "length", ignore.case = TRUE)

  df$rho <- rho
  sc_vec <- score_cate_rho(x0 = rep(0, 3), x_cols = 1:3, rho = rho)
  sc_col <- score_cate_rho(x0 = rep(0, 3), x_cols = 1:3, rho_col = "rho")

  ctx_vec <- score_setup(sc_vec, df); ctx_col <- score_setup(sc_col, df)
  subs <- cate_demo_subsets(nrow(df), K = 25, s = 6)
  ctx_vec$ptr_ctx <- NULL; ctx_col$ptr_ctx <- NULL
  expect_equal(score_eval_batch(sc_vec, subs, df, ctx_vec),
               score_eval_batch(sc_col, subs, df, ctx_col),
               tolerance = 1e-12)
})

# [equiv] dnnk with k=1 equals dnn1 ----------------------------------

test_that("[equiv] dnnk with k = 1 equals dnn1 exactly", {
  demo <- cate_demo_data(n = 80, p = 2, seed = 2)
  df <- demo$data; x0 <- demo$x0; rho <- df$rho_oracle

  sc1 <- score_cate_rho(x0, XC, rho = rho, localizer = "dnn1")
  sck <- score_cate_rho(x0, XC, rho = rho, localizer = "dnnk", k = 1)

  ctx1 <- score_setup(sc1, df); ctxk <- score_setup(sck, df)
  ctx1$ptr_ctx <- NULL; ctxk$ptr_ctx <- NULL

  subs <- cate_demo_subsets(nrow(df), K = 40, s = 7)
  expect_equal(score_eval_batch(sc1, subs, df, ctx1),
               score_eval_batch(sck, subs, df, ctxk),
               tolerance = 0)
})

# [bandw] Kernel auto-bandwidth median rule --------------------------

test_that("[bandw] kernel auto-bandwidth uses true median for even s", {
  demo <- cate_demo_data(n = 50, p = 2, seed = 5)
  df <- demo$data; x0 <- demo$x0; rho <- df$rho_oracle

  sc  <- score_cate_rho(x0, XC, rho = rho, localizer = "kernel", h = NULL)
  ctx <- score_setup(sc, df)

  set.seed(12)                               # even s to exercise median averaging
  subs <- cate_demo_subsets(n = nrow(df), K = 60, s = 8)

  v_cpp <- score_eval_batch(sc, subs, df, ctx)
  ctx$ptr_ctx <- NULL
  v_r   <- score_eval_batch(sc, subs, df, ctx)

  expect_equal(v_cpp, v_r, tolerance = 1e-8)
})

# [compiled] Compiled path behaviour ---------------------------------

test_that("[compiled] compiled batch matches R across localizers", {
  skip_if_not(has_cate_cpp(), "C++ CATE functions not available")

  demo <- cate_demo_data(n = 120, p = 2, design = "hetero", seed = 5)
  df <- demo$data; x0 <- demo$x0; rho <- df$rho_oracle

  for (loc in c("dnn1","dnnk","kernel")) {
    sc  <- score_cate_rho(x0 = x0, x_cols = XC, rho = rho, localizer = loc, k = 5)
    ctx <- score_setup(sc, df)

    set.seed(99)
    subs <- cate_demo_subsets(n = nrow(df), K = 200, s = 8)

    vals_cpp <- score_eval_batch(sc, subs, df, ctx)
    ctx_r <- ctx; ctx_r$ptr_ctx <- NULL
    vals_r <- score_eval_batch(sc, subs, df, ctx_r)

    expect_equal(vals_cpp, vals_r, tolerance = 1e-8, ignore_attr = TRUE)
  }
})

test_that("[compiled] compiled path returns finite values and respects subset size", {
  skip_if_not(has_cate_cpp(), "C++ CATE functions not available")

  demo <- cate_demo_data(n = 90, p = 2, seed = 10)
  df <- demo$data; x0 <- demo$x0; rho <- df$rho_oracle

  sc  <- score_cate_rho(x0 = x0, x_cols = XC, rho = rho, localizer = "dnnk", k = 4)
  ctx <- score_setup(sc, df)

  for (s in c(3L, 5L, 9L)) {
    set.seed(123 + s)
    subs <- cate_demo_subsets(n = nrow(df), K = 50, s = s)
    vals <- score_eval_batch(sc, subs, df, ctx)
    expect_true(all(is.finite(vals)))
    expect_length(vals, nrow(subs))
  }
})
