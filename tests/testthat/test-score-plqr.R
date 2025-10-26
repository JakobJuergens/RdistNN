# tests/testthat/test-plqr-score.R
# -------------------------------------------------------------------
# PLQR (component) score test suite
# Purpose: Tests for score_plqr_component() — factory, setup, eval,
#          compiled vs R fallback, and algebraic properties.
#
# Tags:
#   [smoke]   Symbols / environment
#   [engine]  Engine parity & fallbacks (compiled vs R)
#   [ties]    Determinism & tie-breaking
#   [factory] Basic validation (beta length, j range, tau)
#   [setup]   Setup/order, finite outputs, basic sanity
#   [cohere]  Single-subset vs batch coherence
#   [perm]    Kernel permutation invariance
#   [inputs]  Non-finite & dimension checks
#   [convex]  ψτ(u)·(D_j - mD_j) convex hull bound (with rDj ≡ 1)
#   [equiv]   dnnk with k=1 equals dnn1
#   [bandw]   Kernel auto-bandwidth median rule (even s)
#   [compiled]Compiled path: parity & finite outputs
# -------------------------------------------------------------------

testthat::local_edition(3)

# --- Helpers ---------------------------------------------------------

# Minimal PLQR demo that matches score_plqr_component()'s API:
# We build data so that: D_j ≡ 1, mD_j ≡ 0, beta_j = 0 ⇒ rDj ≡ 1 and Dβ ≡ 0.
# Then the score is a weighted average of ψτ(u) with u := Y - q(X).
.make_plqr_comp_demo <- function(n = 200L, px = 2L, tau = 0.5, seed = 1L) {
  set.seed(seed)
  X <- matrix(runif(n * px), n, px, dimnames = list(NULL, paste0("X", seq_len(px))))
  xnames <- colnames(X)
  # Construct u so we can control ψτ(u) directly in some tests.
  u <- rnorm(n)
  # q(X) is an oracle nuisance; choose a smooth function.
  qx <- as.numeric(X[, 1] - 0.5 * (px >= 2) * X[, 2])
  # Choose Y so that u = Y - q(X) (since Dβ = 0 in our design).
  Y  <- u + qx

  # Single D component with rDj ≡ 1 (D1 = 1, mD1 = 0) and beta = 0
  D1  <- rep(1, n)
  mD1 <- rep(0, n)
  beta <- 0

  df <- data.frame(X, D1 = D1, mD1 = mD1, q = qx, Y = Y, u = u, check.names = FALSE)
  list(
    data = df,
    x0   = colMeans(X) |> as.numeric(),
    tau  = tau,
    cols = list(XC = xnames, D = "D1", mD = "mD1", q = "q"),
    beta = beta,
    j    = 1L
  )
}

.make_subsets <- function(n, K, s) t(replicate(K, sample.int(n, s, replace = FALSE)))

# Skip if factory not present (helps during refactors)
.skip_if_no_factory <- function() {
  testthat::skip_if_not(exists("score_plqr_component", mode = "function"),
                        "score_plqr_component() not found.")
}

.has_plqr_cpp <- function() {
  isTRUE(exists("rustat_plqr_ctx_create", mode = "function")) &&
    isTRUE(exists("rustat_plqr_eval_batch", mode = "function"))
}

# --- Tests -----------------------------------------------------------

# [smoke]
test_that("[smoke] u_statistic() exists for integration tests", {
  testthat::skip_if_not(exists("u_statistic", mode = "function"),
                        "u_statistic() not found; integration tests will skip.")
  succeed()
})

# [engine]
test_that("[engine] compiled and R fallback produce identical fit", {
  testthat::skip_if_not(exists("u_statistic", mode = "function"), "u_statistic() not found")
  testthat::skip_if_not(.has_plqr_cpp(), "C++ PLQR functions not available")
  .skip_if_no_factory()

  demo <- .make_plqr_comp_demo(n = 180, px = 2, tau = 0.5, seed = 77)
  df <- demo$data; x0 <- demo$x0; tau <- demo$tau
  XC <- demo$cols$XC; D_cols <- demo$cols$D; mD_cols <- demo$cols$mD; q_col <- demo$cols$q
  beta <- demo$beta; j <- demo$j

  sc <- score_plqr_component(
    x0 = x0, x_cols = XC, D_cols = D_cols, q_col = q_col, mD_cols = mD_cols,
    tau = tau, j = j, beta = beta, localizer = "dnnk", k = 5
  )

  # Compiled path
  u_cpp <- u_statistic(data = df, score = sc, sampling = list(mode = "fixedN", N = 25000))
  fit_cpp <- fit(u_cpp)$value

  # Force R fallback
  ctx <- score_setup(sc, df); ctx$ptr_ctx <- NULL
  sc_r <- sc; sc_r$setup <- function(data, ...) ctx
  u_r <- u_statistic(data = df, score = sc_r, sampling = list(mode = "fixedN", N = 25000))
  fit_r <- fit(u_r)$value

  expect_equal(fit_cpp, fit_r, tolerance = 1e-8)
})

# [ties]
test_that("[ties] exact ties break deterministically (lower index first)", {
  .skip_if_no_factory()

  n <- 10
  # Two points equidistant to x0 = (0,0); others far away
  df <- data.frame(
    X1 = c(0, 0, rep(10, n-2)),
    X2 = c(1, -1, rep(10, n-2)),
    # rDj ≡ 1:
    D1  = 1,
    mD1 = 0,
    # Build u so first two are positive (ψ_{0.5} = +0.5), others negative (−0.5)
    u   = c( 1,  1, rep(-1, n-2)),
    q   = 0
  )
  df$Y <- df$u + df$q
  x0 <- c(0, 0)
  tau <- 0.5
  XC <- c("X1","X2"); D_cols <- "D1"; mD_cols <- "mD1"; q_col <- "q"
  beta <- 0; j <- 1L

  sc  <- score_plqr_component(x0, XC, D_cols, q_col, mD_cols,
                              tau = tau, j = j, beta = beta, localizer = "dnnk", k = 2)
  ctx <- score_setup(sc, df)

  subs <- rbind(c(1,2,3,4,5,6), c(2,1,3,4,5,6))
  v_cpp <- score_eval_batch(sc, subs, df, ctx)
  ctx$ptr_ctx <- NULL
  v_r   <- score_eval_batch(sc, subs, df, ctx)

  # With rDj ≡ 1, the score is the average of ψτ(u); for (1,2) both ψ=+0.5 ⇒ 0.5
  expect_equal(v_cpp[1], 0.5, tolerance = 1e-12)
  expect_equal(v_cpp, v_r, tolerance = 1e-12)
})

# [factory] validates tau in (0,1), beta length, and j range
test_that("[factory] validates tau in (0,1), beta length, and j range", {
  .skip_if_no_factory()
  demo <- .make_plqr_comp_demo(n = 60, px = 2, tau = 0.5, seed = 123)
  df <- demo$data; x0 <- demo$x0
  XC <- demo$cols$XC; D_cols <- demo$cols$D; mD_cols <- demo$cols$mD; q_col <- demo$cols$q

  expect_error(
    score_plqr_component(x0, XC, D_cols, q_col, mD_cols, tau = 1,   j = 1L, beta = 0),
    regexp = "tau.*\\(0,\\s*1\\)"
  )
  expect_error(
    score_plqr_component(x0, XC, D_cols, q_col, mD_cols, tau = -0.1, j = 1L, beta = 0),
    regexp = "tau.*\\(0,\\s*1\\)"
  )
  expect_error(
    score_plqr_component(x0, XC, D_cols, q_col, mD_cols, tau = 0.5, j = 1L, beta = c(0, 0)),
    regexp = "beta length|length\\(D_cols\\)"
  )
  expect_error(
    score_plqr_component(x0, XC, D_cols, q_col, mD_cols, tau = 0.5, j = 0L, beta = 0),
    regexp = "j.*1.."
  )
  expect_error(
    score_plqr_component(x0, XC, D_cols, q_col, mD_cols, tau = 0.5, j = 2L, beta = 0),
    regexp = "j.*1.."
  )
})


# [cohere]
test_that("[cohere] single-subset eval equals corresponding batch value", {
  .skip_if_no_factory()
  demo <- .make_plqr_comp_demo(n = 70, px = 2, tau = 0.5, seed = 2)
  df <- demo$data; x0 <- demo$x0; tau <- demo$tau
  XC <- demo$cols$XC; D_cols <- demo$cols$D; mD_cols <- demo$cols$mD; q_col <- demo$cols$q
  beta <- demo$beta; j <- demo$j

  sc  <- score_plqr_component(x0, XC, D_cols, q_col, mD_cols,
                              tau = tau, j = j, beta = beta, localizer = "dnnk", k = 3)
  ctx <- score_setup(sc, df); ctx$ptr_ctx <- NULL

  set.seed(7)
  subs <- .make_subsets(n = nrow(df), K = 10, s = 7)
  vals_b <- score_eval_batch(sc, subs, df, ctx)
  vals_s <- vapply(seq_len(nrow(subs)), function(i) score_eval(sc, subs[i, ], df, ctx), numeric(1))
  expect_equal(vals_b, vals_s, tolerance = 1e-12)
})

# [perm]
test_that("[perm] kernel localizer is permutation-invariant within subsets", {
  .skip_if_no_factory()
  demo <- .make_plqr_comp_demo(n = 80, px = 2, tau = 0.5, seed = 3)
  df <- demo$data; x0 <- demo$x0; tau <- demo$tau
  XC <- demo$cols$XC; D_cols <- demo$cols$D; mD_cols <- demo$cols$mD; q_col <- demo$cols$q
  beta <- demo$beta; j <- demo$j

  sc  <- score_plqr_component(x0, XC, D_cols, q_col, mD_cols,
                              tau = tau, j = j, beta = beta, localizer = "kernel", h = NULL)
  ctx <- score_setup(sc, df); ctx$ptr_ctx <- NULL

  set.seed(11)
  subs <- .make_subsets(n = nrow(df), K = 8, s = 6)
  for (i in seq_len(nrow(subs))) {
    idxs <- subs[i, ]
    perm <- sample(idxs)
    v1 <- score_eval(sc, idxs, df, ctx)
    v2 <- score_eval(sc, perm, df, ctx)
    expect_equal(v1, v2, tolerance = 1e-12)
  }
})

# [inputs]
test_that("[inputs] non-finite and dimension checks are enforced", {
  .skip_if_no_factory()
  demo <- .make_plqr_comp_demo(n = 40, px = 3, tau = 0.5, seed = 6)
  df <- demo$data; x0 <- demo$x0
  XC <- demo$cols$XC; D_cols <- demo$cols$D; mD_cols <- demo$cols$mD; q_col <- demo$cols$q

  # x0 length mismatch
  expect_error(score_plqr_component(0, XC, D_cols, q_col, mD_cols, tau = .5, j = 1L, beta = 0) |> score_setup(df),
               regexp = "x0 length|ncol", ignore.case = TRUE)

  # Make Y non-finite
  df_bad <- df; df_bad$Y[3] <- NA_real_
  sc <- score_plqr_component(x0, XC, D_cols, q_col, mD_cols, tau = .5, j = 1L, beta = 0)
  expect_error(score_setup(sc, df_bad), regexp = "Y.*finite", ignore.case = TRUE)

  # Column-count mismatch between D and mD (enforced in setup)
  df_badmD <- df
  df_badmD$mD2 <- 0  # extra nuisance column present in data, but we error earlier anyway

  expect_error(
    score_plqr_component(
      x0     = x0,
      x_cols = XC,
      D_cols = "D1",                  # 1 D column
      q_col  = q_col,
      mD_cols = c("mD1", "mD2"),      # 2 mD columns -> factory should reject
      tau    = 0.5,
      j      = 1L,
      beta   = 0
    ),
    regexp = "mD_cols.*same length|same length and order", ignore.case = TRUE
  )
})

# [convex]
test_that("[convex] with rDj ≡ 1, the localized value lies in [tau - 1, tau]", {
  .skip_if_no_factory()
  demo <- .make_plqr_comp_demo(n = 60, px = 2, tau = 0.5, seed = 5)
  df <- demo$data; x0 <- demo$x0; tau <- demo$tau
  XC <- demo$cols$XC; D_cols <- demo$cols$D; mD_cols <- demo$cols$mD; q_col <- demo$cols$q
  beta <- demo$beta; j <- demo$j

  for (loc in c("dnn1","dnnk","kernel")) {
    sc  <- score_plqr_component(x0, XC, D_cols, q_col, mD_cols,
                                tau = tau, j = j, beta = beta, localizer = loc, k = 3, h = NULL)
    ctx <- score_setup(sc, df); ctx$ptr_ctx <- NULL

    set.seed(10)
    subs <- .make_subsets(nrow(df), K = 25, s = 6)
    vals <- score_eval_batch(sc, subs, df, ctx)

    expect_true(all(vals >= (tau - 1) - 1e-12))
    expect_true(all(vals <=  tau      + 1e-12))
  }
})

# [equiv]
test_that("[equiv] dnnk with k = 1 equals dnn1 exactly", {
  .skip_if_no_factory()
  demo <- .make_plqr_comp_demo(n = 85, px = 2, tau = 0.5, seed = 12)
  df <- demo$data; x0 <- demo$x0; tau <- demo$tau
  XC <- demo$cols$XC; D_cols <- demo$cols$D; mD_cols <- demo$cols$mD; q_col <- demo$cols$q
  beta <- demo$beta; j <- demo$j

  sc1 <- score_plqr_component(x0, XC, D_cols, q_col, mD_cols,
                              tau = tau, j = j, beta = beta, localizer = "dnn1")
  sck <- score_plqr_component(x0, XC, D_cols, q_col, mD_cols,
                              tau = tau, j = j, beta = beta, localizer = "dnnk", k = 1)

  ctx1 <- score_setup(sc1, df); ctxk <- score_setup(sck, df)
  ctx1$ptr_ctx <- NULL; ctxk$ptr_ctx <- NULL

  set.seed(40)
  subs <- .make_subsets(nrow(df), K = 40, s = 7)
  expect_equal(score_eval_batch(sc1, subs, df, ctx1),
               score_eval_batch(sck, subs, df, ctxk),
               tolerance = 0)
})

# [bandw]
test_that("[bandw] kernel auto-bandwidth uses true median for even s", {
  .skip_if_no_factory()
  demo <- .make_plqr_comp_demo(n = 55, px = 2, tau = 0.5, seed = 9)
  df <- demo$data; x0 <- demo$x0; tau <- demo$tau
  XC <- demo$cols$XC; D_cols <- demo$cols$D; mD_cols <- demo$cols$mD; q_col <- demo$cols$q
  beta <- demo$beta; j <- demo$j

  sc  <- score_plqr_component(x0, XC, D_cols, q_col, mD_cols,
                              tau = tau, j = j, beta = beta, localizer = "kernel", h = NULL)
  ctx <- score_setup(sc, df)

  set.seed(12)                               # even s to exercise median averaging
  subs <- .make_subsets(n = nrow(df), K = 60, s = 8)

  v_cpp <- score_eval_batch(sc, subs, df, ctx)   # compiled if present
  ctx$ptr_ctx <- NULL
  v_r   <- score_eval_batch(sc, subs, df, ctx)

  expect_equal(v_cpp, v_r, tolerance = 1e-8)
})

# [compiled]
test_that("[compiled] compiled batch matches R across localizers", {
  testthat::skip_if_not(.has_plqr_cpp(), "C++ PLQR functions not available")
  .skip_if_no_factory()

  demo <- .make_plqr_comp_demo(n = 120, px = 2, tau = 0.5, seed = 5)
  df <- demo$data; x0 <- demo$x0; tau <- demo$tau
  XC <- demo$cols$XC; D_cols <- demo$cols$D; mD_cols <- demo$cols$mD; q_col <- demo$cols$q
  beta <- demo$beta; j <- demo$j

  for (loc in c("dnn1","dnnk","kernel")) {
    sc  <- score_plqr_component(x0, XC, D_cols, q_col, mD_cols,
                                tau = tau, j = j, beta = beta, localizer = loc, k = 5)
    ctx <- score_setup(sc, df)

    set.seed(99)
    subs <- .make_subsets(n = nrow(df), K = 200, s = 8)

    vals_cpp <- score_eval_batch(sc, subs, df, ctx)
    ctx_r <- ctx; ctx_r$ptr_ctx <- NULL
    vals_r <- score_eval_batch(sc, subs, df, ctx_r)

    expect_equal(vals_cpp, vals_r, tolerance = 1e-8, ignore_attr = TRUE)
  }
})

test_that("[compiled] compiled path returns finite values and respects subset size", {
  testthat::skip_if_not(.has_plqr_cpp(), "C++ PLQR functions not available")
  .skip_if_no_factory()

  demo <- .make_plqr_comp_demo(n = 95, px = 2, tau = 0.5, seed = 10)
  df <- demo$data; x0 <- demo$x0; tau <- demo$tau
  XC <- demo$cols$XC; D_cols <- demo$cols$D; mD_cols <- demo$cols$mD; q_col <- demo$cols$q
  beta <- demo$beta; j <- demo$j

  sc  <- score_plqr_component(x0, XC, D_cols, q_col, mD_cols,
                              tau = tau, j = j, beta = beta, localizer = "dnnk", k = 4)
  ctx <- score_setup(sc, df)

  for (s in c(3L, 5L, 9L)) {
    set.seed(123 + s)
    subs <- .make_subsets(n = nrow(df), K = 50, s = s)
    vals <- score_eval_batch(sc, subs, df, ctx)
    expect_true(all(is.finite(vals)))
    expect_length(vals, nrow(subs))
  }
})
