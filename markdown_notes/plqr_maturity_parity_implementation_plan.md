# PLQR Maturity Parity Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Bring `plqr_dnn()` and its supporting methods/tests/docs in the `RdistNN` package up to the same maturity level as `cate_dnn()`, closing concrete, verified gaps in test coverage, internal consistency, and theory-to-code documentation — without adding a paper-simulation harness (explicitly out of scope, see below) and without modifying `cate_dnn()`'s own behavior.

**Architecture:** No new estimator functionality is added. All work is (a) new regression/coverage tests that lock in already-correct existing behavior, (b) two small, targeted code fixes for genuine internal inconsistencies, (c) a new numbered design note plus roxygen/vignette/NEWS updates that make the PLQR variance-estimation design decision traceable to its theoretical justification.

**Tech Stack:** R (testthat edition 3, roxygen2, pkgdown), Rcpp (read-only in this plan), existing RdistNN internal helpers (`.dnn_col_mat`, `.dnn_col_vec`, `.abort_bad_arg`, `.abort_bad_state`, `.validate_*`).

## Global Constraints

- Every step that touches `R/*.R` runs `Rscript -e 'devtools::document()'` afterward if roxygen changed, and the full suite (`Rscript -e 'devtools::test()'`) must show 0 failures before a task is considered done.
- Do not hand-edit `NAMESPACE`, `R/RcppExports.R`, or `src/RcppExports.cpp` — regenerate via `devtools::document()` only.
- Do not modify `cate_dnn()`, `np_dnn()`, or any of their existing test files' assertions — this plan only reads them as reference patterns to port.
- Follow `markdown_notes/design/01`–`05`, `08`, `09`, `13` (already-established package conventions): no `stopifnot()` in user-facing paths, `call. = FALSE` on user-facing diagnostics, validation before heavy computation, argument names in backticks in error messages.
- **Out of scope for this plan** (found during research, deliberately not acted on): a PLQR paper-simulation harness (`R/sim-plqr-*.R` + `paper-plqr-simulation.Rmd`) — the user explicitly chose to skip this, since the Cond_Z_Estimation paper has no PLQR numerical study to mirror; `cate_dnn()`'s orphaned compiled score path (`rustat_cate_ctx_create` unused); `vcov.cate_result()`; aligning `score_cate_rho()`'s deferred-validation timing to `score_plqr_component()`'s eager validation (would break existing CATE tests and isn't required for PLQR parity).

---

## Task 1: Consolidate duplicated PLQR test-data helpers into a shared helper file

Two test files each define their own near-identical local data-generator (`test-user-plqr-dnn.R`'s `.make_plqr_data()` and `test-user-plqr-inference.R`'s `.make_plqr_inference_data()` + `.manual_plqr_sandwich()`). CATE avoids this duplication via `tests/testthat/helper-user-cate-dnn.R`. Later tasks in this plan (2, 3, 6) add new PLQR test files that need the same fixtures, so centralize first.

**Files:**
- Create: `tests/testthat/helper-user-plqr-dnn.R`
- Modify: `tests/testthat/test-user-plqr-dnn.R:12-28` (remove local `.make_plqr_data()`)
- Modify: `tests/testthat/test-user-plqr-inference.R:10-40` (remove local `.make_plqr_inference_data()` and `.manual_plqr_sandwich()`)

**Interfaces:**
- Produces: `.make_plqr_data(n = 80L, p = 2L, p_d = 2L, m = 4L, tau = 0.5)` -> `list(df, X, D, mD, q, Y, X0, beta, tau, x_cols, D_cols, mD_cols)`; `.make_plqr_inference_data(n = 90L, p = 2L, p_d = 2L, m = 2L, tau = 0.5)` -> same shape plus `density`; `.manual_plqr_sandwich(d, w)` -> list with `omega`/`se` used by inference tests. Both are used, unmodified in behavior, by Tasks 2, 3, 6 below.

- [ ] **Step 1: Create the shared helper file with both generators moved verbatim**

```r
# tests/testthat/helper-user-plqr-dnn.R
# Shared helpers for plqr_dnn() user-level tests.

.make_plqr_data <- function(n = 80L, p = 2L, p_d = 2L, m = 4L, tau = 0.5) {
  demo <- plqr_demo_data(n = n, p = p, p_d = p_d, tau = tau, seed = 1)
  df <- demo$data
  X <- as.matrix(df[paste0("X", seq_len(p))])
  D <- as.matrix(df[paste0("D", seq_len(p_d))])
  mD <- as.matrix(df[paste0("mD", seq_len(p_d))])
  q <- as.numeric(df$q_oracle)
  Y <- as.numeric(df$Y)
  X0 <- matrix(rnorm(m * p), m, p)
  list(
    df = df, X = X, D = D, mD = mD, q = q, Y = Y, X0 = X0,
    beta = demo$truth$beta_true, tau = tau,
    x_cols = paste0("X", seq_len(p)),
    D_cols = paste0("D", seq_len(p_d)),
    mD_cols = paste0("mD", seq_len(p_d))
  )
}

.make_plqr_inference_data <- function(n = 90L, p = 2L, p_d = 2L, m = 2L, tau = 0.5) {
  demo <- plqr_demo_data(n = n, p = p, p_d = p_d, tau = tau, seed = 818)
  df <- demo$data
  X <- as.matrix(df[paste0("X", seq_len(p))])
  D <- as.matrix(df[paste0("D", seq_len(p_d))])
  mD <- as.matrix(df[paste0("mD", seq_len(p_d))])
  q <- as.numeric(df$q_oracle)
  Y <- as.numeric(df$Y)
  X0 <- matrix(rnorm(m * p), m, p)
  density <- 0.75 + 0.1 * seq_len(n) / n
  df$q_density <- density
  list(
    df = df, X = X, D = D, mD = mD, q = q, Y = Y, X0 = X0,
    density = density, beta = demo$truth$beta_true, tau = tau,
    x_cols = paste0("X", seq_len(p)),
    D_cols = paste0("D", seq_len(p_d)),
    mD_cols = paste0("mD", seq_len(p_d))
  )
}

.manual_plqr_sandwich <- function(d, w) {
  d_resid <- d$D - d$mD
  S <- d$tau * (1 - d$tau) * crossprod(d_resid, d_resid * w)
  J <- crossprod(d_resid, d$D * as.numeric(d$density * w))
  J_inv <- solve(J)
  omega <- J_inv %*% S %*% t(J_inv)
  vcov <- omega / (1 / sum(w^2))
  list(S = S, J = J, omega = omega, vcov = vcov)
}
```

**Correction recorded during execution (2026-07-08):** the plan's first draft of
`.manual_plqr_sandwich()` above was a paraphrase and did not match the original
function's return fields (`S`/`J`/`omega`/`vcov`, not `se`) or its inline `vcov`
computation. The version shown here has been corrected to match
`tests/testthat/test-user-plqr-inference.R`'s original byte-for-byte, which is
what `test-user-plqr-inference.R:60-65` actually consumes (`manual$S`,
`manual$J`, `manual$omega`, `manual$vcov`). Always verify plan-embedded code
against the actual current file before relying on it verbatim.

- [ ] **Step 2: Remove the now-duplicated definitions from the two original files**

In `tests/testthat/test-user-plqr-dnn.R`, delete lines 12-28 (the `.make_plqr_data <- function(...) {...}` block) — the `skip_if_not(...)` block at lines 1-8 and `set.seed(321)` at line 10 stay.

In `tests/testthat/test-user-plqr-inference.R`, delete the `.make_plqr_inference_data <- function(...) {...}` and `.manual_plqr_sandwich <- function(...) {...}` blocks (originally lines 10-40) — any `skip_if_not(...)` / `set.seed(...)` lines around them stay.

- [ ] **Step 3: Run the full test suite to confirm no regression**

Run: `Rscript -e 'devtools::test(filter = "user-plqr")'`
Expected: same pass count as before this task (0 failures) — this task only relocates code, it does not change behavior.

- [ ] **Step 4: Commit**

```bash
git add tests/testthat/helper-user-plqr-dnn.R tests/testthat/test-user-plqr-dnn.R tests/testthat/test-user-plqr-inference.R
git commit -m "test: centralize PLQR test data helpers into a shared helper file"
```

---

## Task 2: Add custom-distance / custom-rank-weight end-to-end test for PLQR

CATE has `tests/testthat/test-user-cate-dnn-custom.R` proving a user-registered custom `distance` and a user-supplied custom `rank_weight` spec both flow end-to-end through the estimator and agree with a manually computed weighted moment. PLQR has no equivalent — its only "custom distance" coverage uses a canned `distance_l2()` spec object, not a newly `register_distance()`-ed function. Port the CATE pattern onto `plqr_dnn()`.

**Files:**
- Create: `tests/testthat/test-user-plqr-dnn-custom.R`

**Interfaces:**
- Consumes: `.make_plqr_data()` from Task 1's `helper-user-plqr-dnn.R`; existing exported `register_distance()`, `get_distance()`, `weights_dnn()`, `weights_matrix()`, `as_rank_weight()`, `rank_weight_eval_batch()`, `distance_eval_batch()`, `ord_to_ranks()`; existing `plqr_dnn()` matrix interface (`X, D, q, mD, Y, beta, tau, X0, distance, k, s, method, se, return`).

- [ ] **Step 1: Write the test**

```r
# tests/testthat/test-user-plqr-dnn-custom.R

skip_if_not(
  all(vapply(
    c("plqr_dnn", "weights_dnn", "as_rank_weight"),
    exists, logical(1),
    mode = "function"
  )),
  "plqr_dnn() or required helpers not found; skipping file"
)

test_that("custom distance works in plqr_dnn; custom rank_weight works via weights_matrix", {
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

  set.seed(124)
  d <- .make_plqr_data(n = 40L, p = 3L, p_d = 2L, m = 4L, tau = 0.5)

  nm <- paste0("linf_plqr_", as.integer(runif(1, 1e6, 2e6)))
  register_distance(nm, function(X, x0, ctx) {
    apply(abs(X - matrix(x0, nrow(X), ncol(X), TRUE)), 1L, max)
  })
  spec <- get_distance(nm)

  wt <- weights_dnn(
    X = d$X, X0 = d$X0, distance = nm,
    k = 3L, s = 6L, method = "exact",
    return = "W"
  )
  g_manual <- (d$tau - as.numeric(d$Y <= as.vector(d$D %*% d$beta) + d$q)) * (d$D - d$mD)
  moment_manual <- crossprod(g_manual, wt$W)

  out1 <- plqr_dnn(
    X = d$X, D = d$D, q = d$q, mD = d$mD, Y = d$Y,
    beta = d$beta, tau = d$tau, X0 = d$X0,
    distance = nm, k = 3L, s = 6L, method = "exact",
    se = FALSE, return = "W"
  )
  expect_equal(unname(out1$moment), unname(moment_manual), tolerance = 1e-12)

  DX <- distance_eval_batch(spec, d$X, d$X0)
  ORD <- apply(DX, 2L, function(col) order(col, seq_along(col)))
  RKS <- apply(ORD, 2L, ord_to_ranks)

  rw_power <- function(alpha) {
    as_rank_weight(function(r, ctx = NULL) {
      w <- 1 / (as.numeric(r)^alpha)
      w / sum(w)
    }, name = paste0("inv_rank^", alpha))
  }
  rw1 <- rw_power(1)
  rw2 <- rw_power(2)

  W1 <- rank_weight_eval_batch(rw1, RKS)
  W2 <- rank_weight_eval_batch(rw2, RKS)
  moment_rw1 <- crossprod(g_manual, W1)
  moment_rw2 <- crossprod(g_manual, W2)

  expect_true(all(abs(colSums(W1) - 1) < 1e-12))
  expect_true(all(abs(colSums(W2) - 1) < 1e-12))
  expect_false(isTRUE(all.equal(unname(moment_rw1), unname(moment_rw2))))
})
```

- [ ] **Step 2: Run the test**

Run: `Rscript -e 'devtools::test(filter = "user-plqr-dnn-custom")'`
Expected: PASS, 0 failures. This locks in existing correct behavior (`plqr_dnn()` already supports registered custom distances via its `distance` argument, and `rank_weight_eval_batch()` already supports custom rank-weight specs); it should pass without any production-code change. If it fails, that is a genuine, previously-undetected gap — stop and report it rather than adjusting the test to match broken behavior.

- [ ] **Step 3: Commit**

```bash
git add tests/testthat/test-user-plqr-dnn-custom.R
git commit -m "test: add custom-distance/custom-rank-weight coverage for plqr_dnn()"
```

---

## Task 3: Add tied-covariate / rank-boundary test for PLQR

CATE has `tests/testthat/test-user-cate-dnn-ties.R` verifying deterministic tie-breaking (distance ties broken by index) and stable output under the ties, including the `k = s` degenerate-weight warning path. PLQR has no analog. Port it.

**Files:**
- Create: `tests/testthat/test-user-plqr-dnn-ties.R`

**Interfaces:**
- Consumes: `.make_plqr_data()`-style matrix construction (built inline here, matching the CATE ties test's inline style rather than the shared helper, since it needs a specific tied-`X` design); existing `get_distance()`, `distance_eval_batch()`, `plqr_dnn()`.

- [ ] **Step 1: Write the test**

```r
# tests/testthat/test-user-plqr-dnn-ties.R

skip_if_not(
  all(vapply(
    c("plqr_dnn", "get_distance", "distance_eval_batch"),
    exists, logical(1),
    mode = "function"
  )),
  "required APIs not available"
)

test_that("deterministic ties in distances lead to deterministic ordering and stable moment", {
  set.seed(43)
  n <- 10L
  p <- 2L
  p_d <- 1L
  m <- 3L
  X <- cbind(seq_len(n), rep(0, n))
  X0 <- matrix(
    c(
      5, 0,
      5, 1,
      5, -1
    ),
    m, p,
    byrow = TRUE
  )

  D <- matrix(rnorm(n * p_d), n, p_d)
  mD <- matrix(rnorm(n * p_d), n, p_d)
  q <- rnorm(n)
  Y <- rnorm(n)
  beta <- rep(0.3, p_d)

  # k = s is intentional here (uniform weights make tie ordering the only
  # source of variation), so the degenerate-localization warning is expected.
  out <- suppressWarnings(plqr_dnn(
    X = X, D = D, q = q, mD = mD, Y = Y,
    beta = beta, tau = 0.5, X0 = X0,
    distance = "l2", k = 3L, s = 3L, method = "exact",
    se = FALSE, return = c("W", "ord")
  ))
  expect_true("ord" %in% names(out))

  l2 <- get_distance("l2")
  DX <- distance_eval_batch(l2, X, X0)
  exp_ord1 <- order(DX[, 1], seq_along(DX[, 1]))
  expect_equal(out$ord[, 1], exp_ord1)

  out2 <- suppressWarnings(plqr_dnn(
    X = X, D = D, q = q, mD = mD, Y = Y,
    beta = beta, tau = 0.5, X0 = X0,
    distance = "l2", k = 3L, s = 3L, method = "exact",
    se = FALSE
  ))
  expect_equal(out$moment, out2$moment, tolerance = 1e-12)
})
```

- [ ] **Step 2: Run the test**

Run: `Rscript -e 'devtools::test(filter = "user-plqr-dnn-ties")'`
Expected: PASS, 0 failures (regression-lock-in of existing tie-breaking behavior, shared with `weights_dnn()`/`.resolve_ranking()` internals already exercised by the CATE ties test).

- [ ] **Step 3: Commit**

```bash
git add tests/testthat/test-user-plqr-dnn-ties.R
git commit -m "test: add tied-covariate rank-boundary coverage for plqr_dnn()"
```

---

## Task 4: Add missing trivial validation-branch tests to `plqr_dnn()`'s test suite

The following validation branches already exist and work correctly in `plqr_dnn()` (verified against `R/user-dnn-plqr.R`) but have no regression test, unlike their CATE counterparts. Add them as new `test_that()` blocks appended to `tests/testthat/test-user-plqr-dnn.R`.

**Files:**
- Modify: `tests/testthat/test-user-plqr-dnn.R` (append at end of file)

**Interfaces:**
- Consumes: `.make_plqr_data()` from Task 1.

- [ ] **Step 1: Append the following test blocks**

```r
test_that("X0 ncol mismatch vs X is rejected", {
  d <- .make_plqr_data(n = 30L, p = 2L, p_d = 2L, m = 2L, tau = 0.5)
  expect_error(
    plqr_dnn(
      X = d$X, D = d$D, q = d$q, mD = d$mD, Y = d$Y,
      beta = d$beta, tau = d$tau,
      X0 = matrix(0, 2, ncol(d$X) + 1),
      distance = "l2", k = 2L, s = 3L, se = FALSE
    ),
    regexp = "ncol\\(X0\\) .* ncol\\(X\\)", ignore.case = TRUE
  )
})

test_that("se = FALSE bypasses delete-1 jackknife feasibility check for plqr_dnn", {
  d <- .make_plqr_data(n = 10L, p = 2L, p_d = 1L, m = 2L, tau = 0.5)
  expect_no_error(
    plqr_dnn(
      X = d$X, D = d$D, q = d$q, mD = d$mD, Y = d$Y,
      beta = d$beta, tau = d$tau, X0 = d$X0,
      distance = "l2", k = 2L, s = nrow(d$X), method = "exact",
      se = FALSE, jackknife = "loo"
    )
  )
})

test_that("invalid jackknife argument errors clearly for plqr_dnn", {
  d <- .make_plqr_data(n = 20L, p = 2L, p_d = 1L, m = 3L, tau = 0.5)
  expect_error(
    plqr_dnn(
      X = d$X, D = d$D, q = d$q, mD = d$mD, Y = d$Y,
      beta = d$beta, tau = d$tau, X0 = d$X0,
      distance = "l2", k = 2L, s = 3L,
      se = TRUE, jackknife = "bad_option"
    ),
    regexp = "one of"
  )
})

test_that("constant PLQR moment component implies zero jackknife SE", {
  d <- .make_plqr_data(n = 18L, p = 2L, p_d = 1L, m = 2L, tau = 0.5)
  # Force a degenerate case where D == mD everywhere, so g is identically 0
  # regardless of Y/q/beta, and the fast-jackknife SE on `moment` must be 0.
  d$mD <- d$D
  out <- plqr_dnn(
    X = d$X, D = d$D, q = d$q, mD = d$mD, Y = d$Y,
    beta = d$beta, tau = d$tau, X0 = d$X0,
    distance = "l2", k = 2L, s = 4L, method = "exact",
    se = TRUE, jackknife = "fast"
  )
  expect_true(all(abs(out$moment) < 1e-12))
  expect_true(all(abs(out$se) < 1e-12))
})

test_that("matrix interface reports each missing D/q/mD/Y argument individually", {
  d <- .make_plqr_data(n = 15L, p = 2L, p_d = 1L, m = 2L, tau = 0.5)

  expect_error(
    plqr_dnn(X = d$X, q = d$q, mD = d$mD, Y = d$Y, beta = d$beta, tau = d$tau, X0 = d$X0, se = FALSE),
    regexp = "missing: D", ignore.case = TRUE
  )
  expect_error(
    plqr_dnn(X = d$X, D = d$D, mD = d$mD, Y = d$Y, beta = d$beta, tau = d$tau, X0 = d$X0, se = FALSE),
    regexp = "missing: q", ignore.case = TRUE
  )
  expect_error(
    plqr_dnn(X = d$X, D = d$D, q = d$q, Y = d$Y, beta = d$beta, tau = d$tau, X0 = d$X0, se = FALSE),
    regexp = "missing: mD", ignore.case = TRUE
  )
  expect_error(
    plqr_dnn(X = d$X, D = d$D, q = d$q, mD = d$mD, beta = d$beta, tau = d$tau, X0 = d$X0, se = FALSE),
    regexp = "missing: Y", ignore.case = TRUE
  )
})

test_that("data-frame API distinguishes NA vs non-finite x_cols index diagnostics for plqr_dnn", {
  d <- .make_plqr_data(n = 25L, p = 2L, p_d = 1L, m = 2L, tau = 0.5)

  expect_error(
    plqr_dnn(
      data = d$df, x_cols = c(1L, NA_integer_), D_cols = d$D_cols,
      q_col = "q_oracle", mD_cols = d$mD_cols, y_col = "Y",
      beta = d$beta, tau = d$tau, X0 = d$X0, k = 2L, s = 3L, se = FALSE
    ),
    regexp = "x_cols.*NA", ignore.case = TRUE
  )

  expect_error(
    plqr_dnn(
      data = d$df, x_cols = c(1L, Inf), D_cols = d$D_cols,
      q_col = "q_oracle", mD_cols = d$mD_cols, y_col = "Y",
      beta = d$beta, tau = d$tau, X0 = d$X0, k = 2L, s = 3L, se = FALSE
    ),
    regexp = "x_cols.*non-finite", ignore.case = TRUE
  )
})

test_that("data-frame API distinguishes NA vs non-finite covariate diagnostics for plqr_dnn", {
  d <- .make_plqr_data(n = 25L, p = 2L, p_d = 1L, m = 2L, tau = 0.5)
  base <- function(df) {
    tryCatch(
      plqr_dnn(
        data = df, x_cols = d$x_cols, D_cols = d$D_cols,
        q_col = "q_oracle", mD_cols = d$mD_cols, y_col = "Y",
        beta = d$beta, tau = d$tau, X0 = d$X0, k = 2L, s = 3L, se = FALSE
      ),
      error = conditionMessage
    )
  }

  df_na <- d$df
  df_na[[d$x_cols[[1L]]]][2L] <- NA_real_
  df_inf <- d$df
  df_inf[[d$x_cols[[1L]]]][2L] <- Inf

  msg_na <- base(df_na)
  msg_inf <- base(df_inf)

  expect_match(msg_na, "X.*NA", ignore.case = TRUE)
  expect_match(msg_inf, "X.*non-finite", ignore.case = TRUE)
  expect_false(identical(msg_na, msg_inf))
})

test_that("matrix API rejects non-finite X for plqr_dnn", {
  d <- .make_plqr_data(n = 20L, p = 2L, p_d = 1L, m = 2L, tau = 0.5)
  X_bad <- d$X
  X_bad[3L, 1L] <- Inf
  expect_error(
    plqr_dnn(
      X = X_bad, D = d$D, q = d$q, mD = d$mD, Y = d$Y,
      beta = d$beta, tau = d$tau, X0 = d$X0, k = 2L, s = 3L, se = FALSE
    ),
    regexp = "X.*non-finite", ignore.case = TRUE
  )
})
```

- [ ] **Step 2: Run the appended tests**

Run: `Rscript -e 'devtools::test(filter = "user-plqr-dnn")'`
Expected: PASS, 0 failures for the new blocks. If `"matrix interface reports each missing D/q/mD/Y argument individually"` fails because the actual error text differs from `"missing: D"` etc., read the exact message produced by `R/user-dnn-plqr.R`'s `missing_args` block (around line 362-371) and adjust the `regexp=` to match the real wording — do not change the production message just to satisfy the test.

- [ ] **Step 3: Commit**

```bash
git add tests/testthat/test-user-plqr-dnn.R
git commit -m "test: close trivial validation-branch coverage gaps in plqr_dnn()"
```

---

## Task 5: Add missing `X`-finiteness setup-guard test to the PLQR score layer

`score_plqr_component()`'s `setup()` validates `X` for NA/non-finite (`R/score-plqr.R:110-113`), but the existing "finiteness guards fire at setup" test in `tests/testthat/test-score-plqr.R` exercises `Y`, `q`, `D`, `mD` only — never `X`.

**Files:**
- Modify: `tests/testthat/test-score-plqr.R` (inside the existing `test_that("[inputs/setup] finiteness guards fire at setup", ...)` block)

**Interfaces:**
- Consumes: existing test fixture `.get_plqr_fx(n, p, p_d)` -> `fx$ok_clean` with fields `x0, x_cols, d_cols, cols$q, cols$mD, cols$Y, beta_true, tau, data`.

- [ ] **Step 1: Add an `X` case to the existing test**

Add this block inside the existing `test_that("[inputs/setup] finiteness guards fire at setup", ...)` test in `tests/testthat/test-score-plqr.R`, alongside the existing `Y`/`q`/`D`/`mD` cases:

```r
  # X non-finite
  df_bad <- df
  df_bad[[scn$x_cols[1]]][4] <- Inf
  expect_error(score_setup(sc, df_bad), "X cannot contain non-finite values", ignore.case = FALSE)

  # X NA
  df_bad <- df
  df_bad[[scn$x_cols[1]]][6] <- NA_real_
  expect_error(score_setup(sc, df_bad), "X cannot contain `NA` values\\.", ignore.case = FALSE)
```

- [ ] **Step 2: Run the test**

Run: `Rscript -e 'devtools::test(filter = "score-plqr")'`
Expected: PASS, 0 failures.

- [ ] **Step 3: Commit**

```bash
git add tests/testthat/test-score-plqr.R
git commit -m "test: cover X finiteness guard in score_plqr_component() setup"
```

---

## Task 6: Add `.model`-state-corruption tests for `predict.plqr_result()`

`.plqr_validate_prediction_model()` (`R/user-dnn-plqr-methods.R:14-145`) already validates `object$.model`'s `x_mat`, `d_mat`/`mD_mat` shape agreement, `beta` length, `tau`, `k <= s`, `s <= n`, `method`, and `normalize` — field by field, with specific messages — but none of it is under test. CATE's analogous corruption suite is `tests/testthat/test-user-cate-methods.R:209-258`. Port it onto `plqr_result`.

**Files:**
- Modify: `tests/testthat/test-user-plqr-inference.R` (append a new `test_that()` block)

**Interfaces:**
- Consumes: `.make_plqr_inference_data()` from Task 1.

- [ ] **Step 1: Append the test**

```r
test_that("predict.plqr_result validates fitted model state before computing", {
  d <- .make_plqr_inference_data(n = 60L, m = 2L)
  fit <- plqr_dnn(
    X = d$X, D = d$D, q = d$q, mD = d$mD, Y = d$Y,
    beta = d$beta, tau = d$tau, X0 = d$X0,
    method = "exact", k = 2L, s = 6L, se = FALSE
  )

  bad_x <- fit
  bad_x$.model$x_mat[1L, 1L] <- Inf
  expect_error(
    predict(bad_x, X0 = d$X0),
    "object$.model$x_mat cannot contain non-finite values",
    fixed = TRUE
  )

  bad_d <- fit
  bad_d$.model$d_mat <- bad_d$.model$d_mat[-1L, , drop = FALSE]
  expect_error(
    predict(bad_d, X0 = d$X0),
    "object$.model$d_mat/object$.model$mD_mat",
    fixed = TRUE
  )

  bad_mD_cols <- fit
  bad_mD_cols$.model$mD_mat <- cbind(bad_mD_cols$.model$mD_mat, 1)
  expect_error(
    predict(bad_mD_cols, X0 = d$X0),
    "object$.model$mD_mat",
    fixed = TRUE
  )

  bad_beta <- fit
  bad_beta$.model$beta <- bad_beta$.model$beta[-1L]
  expect_error(
    predict(bad_beta, X0 = d$X0),
    "object$.model$beta",
    fixed = TRUE
  )

  bad_k <- fit
  bad_k$.model$k <- bad_k$.model$s + 1L
  expect_error(
    predict(bad_k, X0 = d$X0),
    "object$.model$k must satisfy k <= s",
    fixed = TRUE
  )

  bad_s <- fit
  bad_s$.model$s <- nrow(bad_s$.model$x_mat) + 1L
  expect_error(
    predict(bad_s, X0 = d$X0),
    "object$.model$s must satisfy s <= n",
    fixed = TRUE
  )

  bad_method <- fit
  bad_method$.model$method <- "bogus"
  expect_error(
    predict(bad_method, X0 = d$X0),
    "object$.model$method must be one of",
    fixed = TRUE
  )

  bad_normalize <- fit
  bad_normalize$.model$normalize <- NA
  expect_error(
    predict(bad_normalize, X0 = d$X0),
    "object$.model$normalize must be TRUE or FALSE",
    fixed = TRUE
  )
})
```

- [ ] **Step 2: Run the test**

Run: `Rscript -e 'devtools::test(filter = "user-plqr-inference")'`
Expected: PASS, 0 failures — this exercises existing, already-implemented validation code. If any sub-case fails because the actual message differs from what's quoted above, re-read the corresponding branch in `R/user-dnn-plqr-methods.R:14-145` and correct the `regexp`/`fixed` string in the test to the actual message — do not weaken the assertion to something non-specific.

- [ ] **Step 3: Commit**

```bash
git add tests/testthat/test-user-plqr-inference.R
git commit -m "test: cover .model state validation in predict.plqr_result()"
```

---

## Task 7: Fix `y_col` default inconsistency in `plqr_dnn()`

`plqr_dnn()` declares `y_col = "Y"` (a non-`NULL` default) while every other required data-frame-interface argument (`x_cols`, `D_cols`, `q_col`, `mD_cols`) defaults to `NULL`. This means the "must be provided with `data=` interface" check (`R/user-dnn-plqr.R:319-330`) can never actually flag a missing `y_col` — an absent outcome column instead falls through to a less specific "column not found" error from `.dnn_col_vec()`. This is a real internal inconsistency (not a design fork): every other required column argument in this wrapper, and every required column argument in `cate_dnn()`, uses `NULL` as the "not supplied" sentinel. Fix it for consistency.

**Files:**
- Modify: `R/user-dnn-plqr.R:266` (default value) and `R/user-dnn-plqr.R:239-257` (the `@examples` block, which must still pass `y_col = "Y"` explicitly)
- Test: `tests/testthat/test-user-plqr-dnn.R`

**Interfaces:**
- Produces: `plqr_dnn(..., y_col = NULL, ...)` as the new default signature (was `y_col = "Y"`). Callers that already pass `y_col = "Y"` explicitly (all existing examples and tests do, since the roxygen example already passes it) are unaffected.

- [ ] **Step 1: Write the failing test**

Append to `tests/testthat/test-user-plqr-dnn.R`:

```r
test_that("data-frame API reports y_col as missing when omitted", {
  d <- .make_plqr_data(n = 20L, p = 2L, p_d = 1L, m = 2L, tau = 0.5)
  expect_error(
    plqr_dnn(
      data = d$df, x_cols = d$x_cols, D_cols = d$D_cols,
      q_col = "q_oracle", mD_cols = d$mD_cols,
      beta = d$beta, tau = d$tau, X0 = d$X0, k = 2L, s = 3L, se = FALSE
    ),
    regexp = "missing: y_col", ignore.case = TRUE
  )
})
```

- [ ] **Step 2: Run test to verify it fails**

Run: `Rscript -e 'devtools::test(filter = "user-plqr-dnn")'`
Expected: FAIL — with the current `y_col = "Y"` default, `plqr_dnn()` instead looks for a literal `"Y"` column in `d$df` (which happens to exist in the demo data, so it currently succeeds instead of erroring — the test should fail with "expected error" not being thrown, or produce a different error than expected).

- [ ] **Step 3: Change the default**

In `R/user-dnn-plqr.R`, change line 266 from:

```r
  y_col = "Y",
```

to:

```r
  y_col = NULL,
```

- [ ] **Step 4: Run test to verify it passes**

Run: `Rscript -e 'devtools::test(filter = "user-plqr-dnn")'`
Expected: PASS, 0 failures, including all pre-existing PLQR tests (they all pass `y_col = "Y"` explicitly already, so the default change does not affect them).

- [ ] **Step 5: Regenerate docs and run the full suite**

Run: `Rscript -e 'devtools::document()'`
Run: `Rscript -e 'devtools::test()'`
Expected: 0 failures package-wide.

- [ ] **Step 6: Commit**

```bash
git add R/user-dnn-plqr.R man/plqr_dnn.Rd tests/testthat/test-user-plqr-dnn.R
git commit -m "fix: default y_col to NULL in plqr_dnn() for validation consistency"
```

---

## Task 8: Add `coef.cate_result()`

`plqr_result` has `coef.plqr_result()` (`R/user-dnn-plqr-methods.R:296`) but `cate_result` has no `coef()` method, unlike every other S3 accessor which exists symmetrically on both classes. Add it — this is a one-line, unambiguous accessor (the estimate is `object$tau`), matching design note 13 §8.3's "class consistency" requirement.

**Files:**
- Modify: `R/user-dnn-cate-methods.R` (add new exported function)
- Test: `tests/testthat/test-user-cate-methods.R`

**Interfaces:**
- Produces: `coef.cate_result(object, ...)` -> `object$tau` (named numeric vector, one entry per target, matching the existing `print.cate_result()`/`summary.cate_result()` treatment of `tau`).

- [ ] **Step 1: Write the failing test**

Append to `tests/testthat/test-user-cate-methods.R`:

```r
test_that("coef.cate_result returns the tau vector", {
  obj <- .make_cate_fit(se = FALSE)
  expect_equal(coef(obj$fit), obj$fit$tau)
})
```

- [ ] **Step 2: Run test to verify it fails**

Run: `Rscript -e 'devtools::test(filter = "user-cate-methods")'`
Expected: FAIL with `could not find function "coef.cate_result"` / `no applicable method for 'coef' applied to an object of class "cate_result"`.

- [ ] **Step 3: Implement the method**

Add to `R/user-dnn-cate-methods.R` (near the other simple accessors, e.g. after `print.cate_result` or near `confint.cate_result`):

```r
#' Extract CATE point estimates
#'
#' @param object A `cate_result` object from [cate_dnn()].
#' @param ... Unused; present for S3 consistency.
#' @return Named numeric vector of point estimates (`object$tau`), one entry
#'   per target point.
#' @export
coef.cate_result <- function(object, ...) {
  object$tau
}
```

- [ ] **Step 4: Run test to verify it passes**

Run: `Rscript -e 'devtools::document()'`
Run: `Rscript -e 'devtools::test(filter = "user-cate-methods")'`
Expected: PASS, 0 failures. Confirm `man/coef.cate_result.Rd` was generated and `NAMESPACE` gained `S3method(coef,cate_result)`.

- [ ] **Step 5: Check pkgdown reference coverage**

Run: `Rscript tools/ci/check_pkgdown_reference_coverage.R`
Expected: passes, or reports that `coef.cate_result` needs a `_pkgdown.yml` reference entry — if so, add it alongside the existing `coef.plqr_result` entry in `_pkgdown.yml`.

- [ ] **Step 6: Commit**

```bash
git add R/user-dnn-cate-methods.R man/coef.cate_result.Rd NAMESPACE _pkgdown.yml tests/testthat/test-user-cate-methods.R
git commit -m "feat: add coef.cate_result() for S3 accessor parity with plqr_result"
```

---

## Task 9: Remove dead `has_cate_cpp()` test helper

`tests/testthat/helper-skip.R` defines `has_cate_cpp()`, which checks for `rustat_cate_eval_batch` — a function that does not exist anywhere in the package (only `rustat_cate_ctx_create`, `rustat_cate_g_all`, `rustat_cate_g_indexed` exist). It has zero call sites anywhere in `tests/`. This is dead, stale code left over from an earlier compiled-path design; remove it per "if you are certain something is unused, delete it completely."

**Files:**
- Modify: `tests/testthat/helper-skip.R`

- [ ] **Step 1: Confirm zero call sites**

Run: `grep -rn "has_cate_cpp" tests/`
Expected: only the definition itself in `helper-skip.R`, no call sites elsewhere.

- [ ] **Step 2: Remove the dead function**

Open `tests/testthat/helper-skip.R` and delete the `has_cate_cpp <- function(...) {...}` definition entirely (confirm via Read before editing — do not guess line numbers).

- [ ] **Step 3: Run the full suite**

Run: `Rscript -e 'devtools::test()'`
Expected: 0 failures (removing dead, uncalled code cannot break anything).

- [ ] **Step 4: Commit**

```bash
git add tests/testthat/helper-skip.R
git commit -m "chore: remove dead has_cate_cpp() test helper referencing a nonexistent symbol"
```

---

## Task 10: Add design note 15 — PLQR variance-estimation rationale and nuisance notation mapping

This is the load-bearing documentation deliverable. It records, in the repository's own numbered-design-note convention, **why** `plqr_dnn()`'s `confint()`/`predict()` require the direct `J, S` sandwich (`density`/`density_col`) rather than accepting the ordinary jackknife SE the way `cate_dnn()`/`np_dnn()` do — this was verified during planning to be a theory-driven design choice, not a missing feature:

- The Cond_Z_Estimation paper's active "coupled-score" theory (`Paper/Parts/coupled_score/05_Applications/05b_Coupled_Nonsmooth_Applications.tex`, Corollary `cor:coupled_plqr_worked_point_inference`, condition **PLQR-E**) states explicitly that PLQR uses a direct plug-in `J, S` sandwich rather than the generic jackknife score covariance, and gives the exact formulas `S(bx) = tau(1-tau) * sum_i w_i tilde_D_i tilde_D_i'`, `J(bx) = sum_i w_i q_i tilde_D_i D_i'`, `Omega = J^{-1} S J^{-T}` — which is *exactly* what `.plqr_direct_sandwich()` in `R/user-dnn-plqr.R:56-114` already computes (verified line-by-line during planning: `S` at line 71-72, `J` at line 74, `omega`/`vcov` at lines 90-93).
- JK_GenUStat's jackknife ratio-consistency theorems (`thm:JK_Consistency_compl`/`thm:JK_Consistency_incompl`) are proven for generalized U-statistics, and its own conclusion (`Main_End.tex:9-13`) explicitly states conditional Z-estimation (score-based kernels, which is what the PLQR quantile check-function moment is) is **not resolved** by that paper. So there is no proof anywhere in this repository's source papers that a delete-1 jackknife SE on the PLQR moment translates into a valid coefficient-level CI — which is exactly why the package restricts `confint`/`predict` CI construction to the `density`-gated direct-sandwich route.
- `cate_dnn()`'s jackknife SE, by contrast, *is* covered: its pseudo-outcome score is a function of a single observation (matching the "selector times a function of a single observation" structure JK_GenUStat's TDNN proof (`Supp_DNN_HajekDominance.tex:104-108`) verifies Hájek dominance for), so `thm:JK_Consistency_TDNN` applies to it.

Also record the nuisance-notation mapping between RdistNN's argument names and the paper's symbols, since RdistNN's `q_col`/`q` argument (documented as "the `q_tau(X)` nuisance") is actually playing the role of the paper's `g_0(U)` (the nonparametric partially-linear shape, `U = (X, Z)`), not a "conditional quantile function of X alone" — the paper's own `q_{0,bx}(D,Z)` (a density, entering only the Jacobian) is what RdistNN calls `density`/`density_col`.

**Files:**
- Create: `markdown_notes/design/15_plqr_variance_estimation_rationale.md`
- Modify: `markdown_notes/design/README.md` (add note 15 to the read order and structure)

- [ ] **Step 1: Write the design note**

```markdown
# PLQR Variance-Estimation Rationale and Nuisance Notation

This note records why `plqr_dnn()`'s coefficient-level inference (`confint()`,
`predict(..., se.fit = TRUE)`, `summary()`) is gated on the direct coefficient
sandwich (`density`/`density_col`) rather than accepting the ordinary
delete-1 jackknife SE the way `np_dnn()` and `cate_dnn()` do, and maps
RdistNN's PLQR argument names onto the theory papers' notation. This is a
permanent, theory-driven design decision, not a temporary limitation.

## 1. Why PLQR requires the direct sandwich

`np_dnn()` and `cate_dnn()`'s delete-1 jackknife standard errors are backed
by `JK_GenUStat`'s jackknife ratio-consistency theorems
(`thm:JK_Consistency_compl`, `thm:JK_Consistency_incompl`,
`Main_GenUStat.tex`), specialized to the DNN/TDNN kernel in
`thm:JK_Consistency_TDNN` (`Main_TDNN.tex`). That specialization works
because the underlying score (`Y` for `np_dnn()`, the AIPW/Robinson
pseudo-outcome for `cate_dnn()`) is a function of a single observation —
the same "selector times a function of a single observation" shape the
DNN/TDNN Hájek-dominance proof verifies
(`Supp_DNN_HajekDominance.tex:104-108`).

`JK_GenUStat`'s own conclusion is explicit that this does **not** extend to
conditional Z-estimation for free:

> "A related direction is conditional Z-estimation, where DNN localizer
> weights enter a score-based rather than regression-based kernel; the
> Hájek-dominance machinery applies once the score-kernel first-projection
> variance is controlled ... The present paper does not resolve those
> broader questions." (`Main_End.tex:9-13`)

The PLQR moment (a quantile check-function score) is exactly this
unresolved case. There is no proof anywhere in this repository's source
papers that a delete-1 jackknife SE on the PLQR *moment* has a valid
translation into a coefficient-level (`theta`) confidence interval.

`Cond_Z_Estimation`'s active "coupled-score" theory
(`Paper/Parts/coupled_score/05_Applications/05b_Coupled_Nonsmooth_Applications.tex`)
instead derives a complete, self-contained asymptotic theory for PLQR that
does *not* route through jackknife at all: Corollary
`cor:coupled_plqr_worked_point_inference` gives
`sqrt(neff(x)) (hat theta(bx) - theta_0(bx)) -> N(0, Omega_PLQR(bx))`,
`Omega_PLQR(bx) = J(bx)^{-1} S(bx) J(bx)^{-T}`, with condition **PLQR-E**
specifying the direct plug-in estimators

```
S(bx) = tau(1-tau) * sum_i w_i * tilde_D_i tilde_D_i'
J(bx) = sum_i w_i * q_i * tilde_D_i * D_i'        (q_i estimates the local conditional density)
```

`.plqr_direct_sandwich()` (`R/user-dnn-plqr.R:56-114`) implements exactly
this: `S` at lines 71-72, `J` at line 74, `Omega = J^{-1} S J^{-T}` at lines
90-91, `vcov = Omega / neff` at line 92. This was verified line-by-line
against the paper's formulas during the PLQR-maturity planning pass
(2026-07-08) and matches exactly.

**Consequence for the package's public API:** `confint.plqr_result()` and
`predict.plqr_result(..., se.fit = TRUE)` only produce coefficient SEs/CIs
when `density`/`density_col` was supplied at fit time (the direct-sandwich
route), because that is the only route this repository's theory has proven
valid. The moment-level jackknife SE (`fit$se` when `density` was not
supplied) remains available for backward-compatible moment diagnostics, but
must **not** be treated as a valid basis for a coefficient-level CI on
`theta` — no theorem in this repository's paper corpus covers that
translation.

## 2. Nuisance notation mapping

RdistNN's PLQR argument names predate the paper's current "coupled-score"
notation and use a different (looser) naming convention. This table is the
authoritative mapping for anyone connecting code to theory:

| RdistNN name         | Paper symbol                  | Meaning                                                                 | Domain (per paper)          |
|-----------------------|-------------------------------|--------------------------------------------------------------------------|------------------------------|
| `q_col` / `q`          | `g_0(U)`                       | Nonparametric partially-linear shape nuisance; enters the threshold as `D'theta + g_0(U)`. It is **not itself a quantile function** — the sum `D'theta_0 + g_0` is the conditional quantile. The name `q_tau(X)` in RdistNN's docs/examples is a legacy label kept for backward compatibility; treat it as `g_0`, not as "the quantile." | `U = (X, Z)`, the full adjustment state — not `X` alone. |
| `mD_cols` / `mD`       | `m_0(U) = E[D \| U]`            | Conditional mean used to residualize the design: `tilde_D = D - m_0(U)`. | `U = (X, Z)` |
| `density` / `density_col` | `q_{0,bx}(D, Z)`            | Conditional density of `Y` at the fitted quantile value; enters **only** the Jacobian `J(bx)`, never the moment or `S(bx)`. | function of `(D, Z)` given `X = bx` |

Both `q_col`/`q` and `mD_cols`/`mD` are supplied externally by the user
(no nuisance is estimated internally by `plqr_dnn()`), so nothing in the
implementation actually restricts them to depend on `X` alone — the
restriction was only in the documentation's phrasing. Task 11 of the
2026-07-08 PLQR-maturity plan corrects the roxygen wording to match this
table; this note is the source of truth going forward.

## 3. Scope note

This note does not attempt to prove Hájek dominance or jackknife
consistency for the PLQR score kernel — per `Main_End.tex`, that remains
open in the source papers. It only records why the package's existing
`confint`/`predict` gating on `density`/`density_col` is the theoretically
correct behavior, not a placeholder limitation, so a future contributor
does not "fix" the asymmetry with `cate_dnn()` by adding an unproven
jackknife-based coefficient CI.
```

- [ ] **Step 2: Add the note to the design-note README's read order**

In `markdown_notes/design/README.md`, add `15_plqr_variance_estimation_rationale.md` under the "Execution" section (after `14_integrated_tuning_grid_reuse.md`) in both the numbered "Read order" list and the "Structure" list, matching the existing format for note 14.

- [ ] **Step 3: Commit**

```bash
git add markdown_notes/design/15_plqr_variance_estimation_rationale.md markdown_notes/design/README.md
git commit -m "docs: add design note on PLQR variance-estimation rationale and nuisance notation"
```

---

## Task 11: Align `plqr_dnn()` roxygen with the paper's nuisance notation

Update the moment-formula and `@param` documentation in `R/user-dnn-plqr.R` to match the mapping recorded in Task 10's design note, and cross-reference it.

**Files:**
- Modify: `R/user-dnn-plqr.R:116-236` (roxygen block for `plqr_dnn()`)

- [ ] **Step 1: Update the `@details` nuisance paragraph**

In `R/user-dnn-plqr.R`, replace the current paragraph at lines 138-140:

```r
#' No nuisance models are fit internally; this wrapper assumes `q_tau(X)` and
#' `m_D(X)` are supplied. Density values for direct coefficient inference are
#' also supplied by the user rather than estimated internally.
```

with:

```r
#' No nuisance models are fit internally; this wrapper assumes `q_tau(X)` and
#' `m_D(X)` are supplied. Despite the `q_tau(X)` name (kept for backward
#' compatibility), this nuisance plays the role of the partially-linear
#' quantile model's nonparametric shape function \eqn{g_0(U)} in the
#' conditional Z-estimation theory this estimator implements: the threshold
#' is \eqn{D'\beta + g_0(U)}, and \eqn{g_0} alone is not itself a quantile.
#' Both `q_tau(X)` and `m_D(X)` may depend on covariates beyond `x_cols`
#' (the paper's \eqn{U = (X, Z)}); the package does not restrict their
#' domain, since they are supplied pre-computed. Density values for direct
#' coefficient inference are also supplied by the user rather than estimated
#' internally; they play the role of the paper's local conditional density
#' \eqn{q_{0,bx}(D,Z)}, entering only the coefficient Jacobian. See
#' `markdown_notes/design/15_plqr_variance_estimation_rationale.md` for the
#' full nuisance-notation mapping and for why coefficient-level inference
#' requires this direct-sandwich route rather than the jackknife SE used by
#' [cate_dnn()]/[np_dnn()].
```

- [ ] **Step 2: Regenerate docs**

Run: `Rscript -e 'devtools::document()'`
Expected: `man/plqr_dnn.Rd` updates; no errors.

- [ ] **Step 3: Run the full test suite**

Run: `Rscript -e 'devtools::test()'`
Expected: 0 failures — this is a documentation-only change.

- [ ] **Step 4: Commit**

```bash
git add R/user-dnn-plqr.R man/plqr_dnn.Rd
git commit -m "docs: align plqr_dnn() nuisance documentation with paper notation"
```

---

## Task 12: Document the direct-sandwich requirement on the inference methods themselves

`confint.plqr_result()`, `predict.plqr_result()`, and `summary.plqr_result()` currently say inference is "not available" / tell the user to "Refit with `density` or `density_col`" without saying *why* jackknife SE isn't an accepted fallback. Add the rationale to their roxygen and to the runtime error/message text, citing design note 15. Preserve existing test substrings exactly (`"density"` and `"not available"` are matched by existing tests with `ignore.case = TRUE`).

**Files:**
- Modify: `R/user-dnn-plqr-methods.R` (roxygen for `confint.plqr_result`, `predict.plqr_result`, `summary.plqr_result`; the `.abort_bad_state()` message in `predict.plqr_result()` at lines 506-512)

**Interfaces:**
- Consumes: existing `tests/testthat/test-user-plqr-inference.R` assertions matching `"density"` (ignore.case) and `"not available"` (ignore.case) — must keep matching after this change.

- [ ] **Step 1: Extend the `predict.plqr_result()` runtime message**

In `R/user-dnn-plqr-methods.R:506-512`, change:

```r
    .abort_bad_state(paste0(
      "Direct PLQR coefficient standard errors are not available for prediction ",
      "because this object does not retain density values. Refit with `density` ",
      "or `density_col`."
    ))
```

to:

```r
    .abort_bad_state(paste0(
      "Direct PLQR coefficient standard errors are not available for prediction ",
      "because this object does not retain density values. The jackknife SE on ",
      "the moment is not a validated basis for a coefficient-level interval here ",
      "(see markdown_notes/design/15_plqr_variance_estimation_rationale.md); ",
      "refit with `density` or `density_col` for coefficient-level inference."
    ))
```

- [ ] **Step 2: Add a roxygen `@details` note to `confint.plqr_result()`, `predict.plqr_result()`, and `summary.plqr_result()`**

For each of the three functions in `R/user-dnn-plqr-methods.R`, add (or extend an existing `@details` section with) this sentence:

```r
#' @details
#' Coefficient-level standard errors and intervals require the direct PLQR
#' sandwich (fit with `density` or `density_col`); the moment-level jackknife
#' SE available on fits without a density input has no proven translation
#' into a coefficient-level interval in this package's source theory — see
#' `markdown_notes/design/15_plqr_variance_estimation_rationale.md`.
```

- [ ] **Step 3: Run the existing inference test file to confirm no regression**

Run: `Rscript -e 'devtools::test(filter = "user-plqr-inference")'`
Expected: PASS, 0 failures — the existing `expect_error(..., "density", ignore.case = TRUE)` and `expect_error(..., "not available", ignore.case = TRUE)` assertions still match, since both words remain in the extended message.

- [ ] **Step 4: Regenerate docs**

Run: `Rscript -e 'devtools::document()'`

- [ ] **Step 5: Commit**

```bash
git add R/user-dnn-plqr-methods.R man/confint.plqr_result.Rd man/predict.plqr_result.Rd man/summary.plqr_result.Rd
git commit -m "docs: explain why PLQR coefficient inference requires the direct sandwich"
```

---

## Task 13: Extend the `conditional-z-estimation.Rmd` vignette with a PLQR theory section

The vignette currently only lists `plqr_dnn()` in a one-line bullet (`vignettes/conditional-z-estimation.Rmd:24`) with no grounding in the nuisance definitions or the variance-estimation design. Add a short, concrete section.

**Files:**
- Modify: `vignettes/conditional-z-estimation.Rmd`

- [ ] **Step 1: Read the current vignette structure**

Before editing, read the full current file to find the right insertion point (likely after the existing PLQR bullet/section, before or alongside the existing runnable example content) and to match its existing heading level and prose style.

- [ ] **Step 2: Add a new subsection**

Insert a subsection (using the vignette's existing heading style) titled something like "PLQR nuisances and inference" with this content, adapted to fit the surrounding prose:

```markdown
### PLQR nuisances and inference

`plqr_dnn()` localizes the moment

$$
\mathbb{E}\!\left[
  \left(\tau - \mathbf{1}\{Y \le D^\top \beta + q_\tau(X)\}\right)
  \left(D - m_D(X)\right)
  \;\middle|\; X \approx x_0
\right]
$$

Despite its name, `q_tau(X)` is the partially-linear model's nonparametric
shape nuisance (paper notation $g_0(U)$), not itself a quantile — the sum
$D'\theta_0 + g_0$ is the conditional quantile. `m_D(X)` residualizes the
design ($m_0(U) = E[D \mid U]$). Both are supplied by the user; no nuisance
model is fit internally.

Coefficient-level standard errors and confidence intervals require
supplying `density`/`density_col` (an estimate of the local conditional
density of `Y` at the fitted quantile). This is not an optional
convenience: the package's source theory (see
`markdown_notes/design/15_plqr_variance_estimation_rationale.md`) proves
coefficient-level inference only for this direct sandwich route, not for a
jackknife SE on the moment. Fits without a density input can still report
`se` (a jackknife SE on the *moment*, useful as a diagnostic) but
`confint()`/`predict(se.fit = TRUE)` will error, by design, until refit
with a density input.
```

- [ ] **Step 3: Confirm the vignette still builds**

Run: `Rscript -e 'devtools::build_vignettes()'` (or the repo's usual vignette-check command if different — check for one in `tools/ci/` first)
Expected: no errors.

- [ ] **Step 4: Commit**

```bash
git add vignettes/conditional-z-estimation.Rmd
git commit -m "docs(vignette): explain PLQR nuisance notation and the direct-sandwich requirement"
```

---

## Task 14: Add missing NEWS.md entries

Three merged, user-visible PLQR features have no NEWS.md entry: direct sandwich inference, the prediction method, and tuning-grid mode (for both `cate_dnn()`/`np_dnn()`/`plqr_dnn()`). Add them, plus a line recording the new design note and the simulation-harness scope decision.

**Files:**
- Modify: `NEWS.md` (add a new top section, since `## RdistNN (development)` already exists and these are unreleased additions)

- [ ] **Step 1: Read the current top of NEWS.md**

Read `NEWS.md` lines 1-10 to find exactly where to insert (immediately under the `## RdistNN (development)` heading, above the existing `### Plotting layer (Phase 3)` section, as a new subsection).

- [ ] **Step 2: Insert a new subsection**

```markdown
### PLQR inference, prediction, and tuning-grid mode

- Added the direct PLQR coefficient sandwich (`density`/`density_col`):
  `plqr_dnn()` now reports coefficient-level `se`, `vcov`, `omega`,
  `jacobian`, `score_covariance`, and `neff` when a local conditional
  density is supplied, matching the direct `J, S` sandwich this
  estimator's source theory specifies (see
  `markdown_notes/design/15_plqr_variance_estimation_rationale.md`).
  `confint.plqr_result()` and `summary.plqr_result()` report coefficient
  inference only through this route.
- Added `predict.plqr_result()`, including coefficient-level `se.fit`/
  `interval = "confidence"` support gated on the same direct-sandwich
  requirement.
- Added `tuning_mode = c("default", "grid")` to `cate_dnn()`, `np_dnn()`,
  and `plqr_dnn()`: grid mode evaluates a tuning object's feasible kernel-
  order candidates and attaches a tidy `tuning_grid` payload
  (`candidates`, `estimates`, `diagnostics`, `meta`), reusing cached ranks
  and rank-weight tables rather than refitting per candidate. See
  `markdown_notes/design/14_integrated_tuning_grid_reuse.md`.
- Recorded, but deliberately did not build, a PLQR analogue of the
  `sim-cate-*.R` paper-simulation harness: the Cond_Z_Estimation paper has
  no PLQR numerical study to mirror, so a simulation harness would require
  inventing a DGP and comparator methods from scratch. Left as an explicit,
  flagged future-work item rather than undertaken speculatively.
```

- [ ] **Step 3: Commit**

```bash
git add NEWS.md
git commit -m "docs: add NEWS.md entries for PLQR direct sandwich, prediction, and grid mode"
```

---

## Task 15: Final full verification pass

Mirror the verification protocol already used for recent PLQR-related commits in this repository (see `git log -1` message for commit `71c1ef0`), so this plan's changes are held to the same bar.

**Files:** none (verification only)

- [ ] **Step 1: Full test suite**

Run: `Rscript -e 'devtools::test()'`
Expected: 0 failures, 0 warnings (skips are acceptable if pre-existing and unrelated to this plan).

- [ ] **Step 2: Lint**

Run: `Rscript tools/ci/lint.R TRUE`
Expected: passes.

- [ ] **Step 3: Generated-file drift check**

Run: `Rscript tools/ci/check_generated_file_drift.R`
Expected: passes (confirms `NAMESPACE`/`man/*.Rd` were regenerated and committed for Tasks 7, 8, 11, 12).

- [ ] **Step 4: pkgdown reference coverage**

Run: `Rscript tools/ci/check_pkgdown_reference_coverage.R`
Expected: passes.

- [ ] **Step 5: Declared-test-dependencies check**

Run: `Rscript tools/ci/check_tests_declared_deps.R`
Expected: passes.

- [ ] **Step 6: Full package check**

Run: `Rscript -e 'res <- rcmdcheck::rcmdcheck(args = c("--no-manual", "--as-cran"), build_args = c("--no-manual", "--compact-vignettes=gs+qpdf"), error_on = "never", check_dir = tempdir()); print(res)'`
Expected: 0 errors, 0 warnings (notes are acceptable if consistent with the pre-existing baseline noted in recent commit `71c1ef0`: "0 errors, 0 warnings, 5 notes").

- [ ] **Step 7: Final commit and summary**

If all checks pass and no fixups were needed beyond what earlier tasks already committed, this step is a no-op. If any check required a fix, commit it:

```bash
git add -A
git commit -m "fix: address CI guardrail findings from PLQR maturity parity pass"
```

No `push` is performed as part of this plan — confirm with the user before pushing to `origin/master`, per this repository's own convention of committing directly on the default branch but treating pushes as a distinct, explicit step.
