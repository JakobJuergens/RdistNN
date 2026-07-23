# Meaningful Test-Coverage Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Repair installed-layout test gaps and add focused tests for public scientific contracts, user-facing result methods, and PLQR compiled/fallback equivalence.

**Architecture:** Keep all changes in `tests/testthat/` except for this planning record. Repair existing tests where the contract already has a natural home, add one focused NPR demo-data file, and extend the existing PLQR and GMM suites without changing package behavior. Use installed-package execution to distinguish public namespace behavior from `pkgload::load_all(export_all = TRUE)` behavior.

**Tech Stack:** R, `testthat` edition 3, `pkgload`, `rcmdcheck`, and `covr`.

## Global Constraints

- Coverage percentage is a diagnostic, not an acceptance target.
- Every added expectation must protect an active public, scientific, or engine-equivalence contract.
- Do not add tests for benchmark timing utilities, package-startup hooks, optional third-party benchmark integrations, or the unused `.validate_same_n()` helper solely to increase coverage.
- Tests must work both from the source tree and against an installed package.
- Restore `options(RdistNN.use_compiled)` after every engine-parity test.
- Make no production-code changes unless a new test exposes an actual defect.
- Execute this plan inline because the active workspace instructions do not authorize subagent delegation.

---

### Task 1: Repair installed-layout test execution

**Files:**
- Modify: `tests/testthat/test-user-cate-methods.R:1-8`
- Modify: `tests/testthat/test-user-cate-dnn-custom.R:11-48`
- Modify: `tests/testthat/test-user-plqr-dnn-custom.R:9-46`

**Interfaces:**
- Consumes: exported `cate_dnn()`, `plqr_dnn()`, `register_distance()`, `get_distance()`, `distance_eval_batch()`, `weights_dnn()`, `as_rank_weight()`, and `rank_weight_eval_batch()`.
- Produces: installed-layout tests that fail when public S3 dispatch or custom extension points fail, instead of skipping because an internal symbol is hidden.

- [ ] **Step 1: Record the current installed-layout failure evidence**

Run against an installed copy of the current package:

```bash
R_LIBS=/tmp/rdistnn-installed-audit.8ONKpd Rscript --vanilla -e \
  'library(testthat); test_dir("tests/testthat", filter = "user-cate-methods|user-cate-dnn-custom|user-plqr-dnn-custom", reporter = "summary", package = "RdistNN", load_package = "installed", stop_on_failure = FALSE)'
```

Expected current result: `test-user-cate-methods.R` reports one file-level skip and zero expectations; the CATE and PLQR custom tests report `required APIs not available`.

- [ ] **Step 2: Remove the invalid CATE S3 visibility guard**

Delete the file-level `skip_if_not(...)` block from `test-user-cate-methods.R`. Do not replace it with another skip. The existing calls through `summary()`, `confint()`, and `predict()` are the public dispatch contract.

- [ ] **Step 3: Remove the custom tests' internal-helper dependency**

In both custom test files, remove `ord_to_ranks` and the unused `weights_matrix` from the inner availability list. Replace:

```r
RKS <- apply(ORD, 2L, ord_to_ranks)
```

with:

```r
RKS <- apply(ORD, 2L, function(ord) {
  ranks <- integer(length(ord))
  ranks[ord] <- seq_along(ord)
  ranks
})
```

This conversion is test-local and reproduces the rank permutation without reaching into the package namespace.

- [ ] **Step 4: Run the repaired files against the installed package**

Run the Step 1 command again.

Expected: the three files execute expectations; neither namespace-related skip reason appears.

- [ ] **Step 5: Run the same files from the source tree**

```bash
Rscript --vanilla -e 'devtools::test(filter = "user-cate-methods|user-cate-dnn-custom|user-plqr-dnn-custom", reporter = "summary")'
```

Expected: zero failures and zero errors.

### Task 2: Add NPR demo-data and shared subset contracts

**Files:**
- Create: `tests/testthat/test-demo-data-npr.R`

**Interfaces:**
- Consumes: `np_demo_data(n, p, design, sigma, seed, x0, inject, idx_inject)`, `np_demo_subsets(n, K, s, seed)`, `cate_demo_subsets(n, K, s, seed)`, and `plqr_demo_subsets(n, K, s, seed)`.
- Produces: reproducibility, oracle, schema, injection, and subset-sampling regression tests.

- [ ] **Step 1: Add the schema and deterministic-oracle tests**

Create the file with focused tests of the documented return structure and every supported design:

```r
test_that("np_demo_data returns its documented schema reproducibly", {
  first <- np_demo_data(n = 24L, p = 3L, design = "sinusoidal", seed = 71L)
  second <- np_demo_data(n = 24L, p = 3L, design = "sinusoidal", seed = 71L)

  expect_named(first, c("data", "x0", "truth", "m_x0"))
  expect_named(first$data, c("X1", "X2", "X3", "Y", "m_oracle"))
  expect_equal(first$data, second$data)
  expect_equal(first$x0, rep(0, 3L))
  expect_equal(first$m_x0, second$m_x0)
})

test_that("np_demo_data zero-noise outcomes equal the oracle for every design", {
  x0 <- c(0.2, -0.1, 0.4)
  for (design in c("linear", "sinusoidal", "bump", "constant")) {
    demo <- np_demo_data(
      n = 30L, p = 3L, design = design, sigma = 0,
      seed = 72L, x0 = x0
    )
    X <- as.matrix(demo$data[paste0("X", 1:3)])
    expect_equal(demo$data$Y, demo$data$m_oracle)
    expect_equal(demo$truth$m_fn(X), demo$data$m_oracle)
    expect_equal(demo$m_x0, demo$truth$m_fn(matrix(x0, nrow = 1L)))
  }
})
```

- [ ] **Step 2: Add injection and validation tests**

Add assertions that `inject = "na"` and `inject = "inf"` alter exactly the requested `Y` position, and that invalid `sigma`, `x0`, `inject`, and `idx_inject` values return argument-specific errors.

```r
test_that("np_demo_data supports diagnostic injection and validates it", {
  na_demo <- np_demo_data(n = 12L, p = 2L, seed = 73L,
                          inject = "na", idx_inject = 4L)
  inf_demo <- np_demo_data(n = 12L, p = 2L, seed = 73L,
                           inject = "inf", idx_inject = 5L)
  expect_true(is.na(na_demo$data$Y[4L]))
  expect_true(is.infinite(inf_demo$data$Y[5L]))
  expect_error(np_demo_data(sigma = -1), "sigma")
  expect_error(np_demo_data(p = 2L, x0 = 0), "x0")
  expect_error(np_demo_data(inject = "bad"), "inject")
  expect_error(np_demo_data(n = 5L, inject = "na", idx_inject = 6L), "idx_inject")
})
```

- [ ] **Step 3: Add a table-driven subset-helper contract**

```r
test_that("public demo subset helpers return deterministic valid subsets", {
  helpers <- list(
    np = np_demo_subsets,
    cate = cate_demo_subsets,
    plqr = plqr_demo_subsets
  )
  for (helper in helpers) {
    first <- helper(n = 12L, K = 7L, s = 4L, seed = 74L)
    second <- helper(n = 12L, K = 7L, s = 4L, seed = 74L)
    expect_equal(dim(first), c(7L, 4L))
    expect_equal(first, second)
    expect_true(all(first >= 1L & first <= 12L))
    expect_true(all(apply(first, 1L, function(idx) length(unique(idx)) == 4L)))
    expect_true(all(apply(first, 1L, function(idx) identical(idx, sort(idx)))))
    expect_error(helper(n = 3L, K = 2L, s = 4L), "s")
  }
})
```

- [ ] **Step 4: Prove the oracle test is discriminating**

Temporarily change the local test expectation from `demo$data$m_oracle` to `demo$data$m_oracle + 1`, run the file, and confirm the zero-noise test fails with a value mismatch. Restore the correct expectation with `apply_patch` immediately.

- [ ] **Step 5: Run the completed file**

```bash
Rscript --vanilla -e 'devtools::test(filter = "demo-data-npr", reporter = "summary")'
```

Expected: zero failures and zero errors.

### Task 3: Cover the public PLQR score-matrix formula and validation

**Files:**
- Modify: `tests/testthat/test-user-plqr-dnn.R` near the existing `plqr_score_matrix` diagnostics.

**Interfaces:**
- Consumes: `plqr_score_matrix(data, y_col, D_cols, g_col, mD_cols, theta, tau)` and `.make_plqr_data()` from `helper-user-plqr-dnn.R`.
- Produces: direct formula, selector equivalence, shape, and input-contract tests.

- [ ] **Step 1: Add the direct formula and selector-equivalence test**

```r
test_that("plqr_score_matrix matches the orthogonal score formula", {
  d <- .make_plqr_data(n = 32L, p = 2L, p_d = 2L, m = 1L, tau = 0.35)
  score_names <- plqr_score_matrix(
    data = d$df, y_col = "Y", D_cols = d$D_cols,
    g_col = "g_oracle", mD_cols = d$mD_cols,
    theta = d$theta_eval, tau = d$tau
  )
  selected <- match(c("Y", d$D_cols, "g_oracle", d$mD_cols), names(d$df))
  score_indices <- plqr_score_matrix(
    data = d$df, y_col = selected[1L], D_cols = selected[2:3],
    g_col = selected[4L], mD_cols = selected[5:6],
    theta = d$theta_eval, tau = d$tau
  )
  threshold <- as.vector(d$D %*% d$theta_eval) + d$g
  expected <- (d$tau - as.numeric(d$Y <= threshold)) * (d$D - d$mD)

  expect_equal(dim(score_names), c(32L, 2L))
  expect_equal(unname(score_names), unname(expected), tolerance = 1e-12)
  expect_equal(score_indices, score_names)
})
```

- [ ] **Step 2: Add meaningful public-input failures**

Add separate expectations for non-data-frame input, invalid `tau`, mismatched `mD_cols`, wrong theta length, and a non-finite selected nuisance column. Match argument names rather than full message prose.

- [ ] **Step 3: Prove the formula test is discriminating**

Temporarily replace the expected indicator `d$Y <= threshold` with `d$Y < threshold - 1` and confirm the formula test fails. Restore the documented formula immediately.

- [ ] **Step 4: Run the focused PLQR wrapper tests**

```bash
Rscript --vanilla -e 'devtools::test(filter = "user-plqr-dnn", reporter = "summary")'
```

Expected: zero failures and zero errors.

### Task 4: Strengthen CATE and GMM result-method contracts

**Files:**
- Modify: `tests/testthat/test-user-cate-methods.R`
- Modify: `tests/testthat/test-user-dnn-gmm.R`

**Interfaces:**
- Consumes: `print.cate_result`, `summary.cate_result`, `print.summary.cate_result`, `print.gmm_result`, `summary.gmm_result`, and `print.summary.gmm_result` through public S3 dispatch.
- Produces: stable object-class, invisibility, coefficient-schema, SE-state, and overidentification-display tests.

- [ ] **Step 1: Add CATE print-state tests**

For fits constructed by `.make_cate_fit(se = TRUE/FALSE)`, capture `print(fit)` and assert `<cate_result>`, the correct SE-state phrase, and `expect_identical(invisible(print(fit)), fit)` using `capture.output(returned <- print(fit))`.

- [ ] **Step 2: Add GMM result and summary output tests**

Use the existing overidentified DR-IV fixture (`q = 2`) and a just-identified fixture (`q = 1`). Assert:

```r
summary_overid <- summary(fit_overid)
expect_s3_class(summary_overid, "summary.gmm_result")
expect_named(
  summary_overid$coefficients,
  c("target", "term", "estimate", "std.error", "statistic", "p.value", "conf.low", "conf.high")
)
expect_named(summary_overid$overid, c("target", "J", "df", "p.value", "neff"))
expect_true(any(is.finite(summary_overid$overid$J)))
```

Capture both print methods, check their identifying headers, verify the overidentified summary prints `Overidentification diagnostics:`, and verify the just-identified summary omits that heading because all `J` values are `NA`.

- [ ] **Step 3: Add the no-SE display branch**

Fit one small GMM object with `se = FALSE`, capture `print(fit_no_se)`, assert `SEs not requested`, and verify the coefficient table's inferential columns are `NA`.

- [ ] **Step 4: Run the focused method tests**

```bash
Rscript --vanilla -e 'devtools::test(filter = "user-cate-methods|user-dnn-gmm", reporter = "summary")'
```

Expected: zero failures and zero errors.

### Task 5: Add PLQR compiled/fallback parity

**Files:**
- Modify: `tests/testthat/test-compiled-fallback-parity.R`

**Interfaces:**
- Consumes: `score_plqr_component()`, `score_setup()`, `score_eval()`, `score_eval_batch()`, `.rdistnn_compiled_fun()`, and the `RdistNN.use_compiled` option.
- Produces: equality evidence for compiled and pure-R PLQR score construction and subset aggregation.

- [ ] **Step 1: Write the parity test**

```r
test_that("PLQR score contributions agree across compiled and R engines", {
  skip_if(
    is.null(.rdistnn_compiled_fun("rustat_plqr_ctx_create")) ||
      is.null(.rdistnn_compiled_fun("rustat_subset_row_sums_g")),
    "Compiled PLQR score symbols are unavailable."
  )
  old <- getOption("RdistNN.use_compiled")
  on.exit(options(RdistNN.use_compiled = old), add = TRUE)

  demo <- plqr_demo_data(n = 36L, p = 2L, p_d = 2L, tau = 0.4, seed = 75L)
  score <- score_plqr_component(
    x0 = demo$x0, x_cols = c("X1", "X2"),
    D_cols = c("D1", "D2"), g_col = "g_oracle",
    mD_cols = c("mD1", "mD2"), y_col = "Y",
    theta = demo$truth$theta_true, j = 2L, tau = 0.4, s = 3L
  )
  subsets <- rbind(c(1L, 2L, 3L), c(4L, 8L, 12L), c(20L, 21L, 36L))

  options(RdistNN.use_compiled = TRUE)
  compiled_ctx <- score_setup(score, demo$data)
  compiled <- list(
    contribution = compiled_ctx$score,
    single = score_eval(score, subsets[1L, ], demo$data, compiled_ctx),
    batch = score_eval_batch(score, subsets, demo$data, compiled_ctx)
  )

  options(RdistNN.use_compiled = FALSE)
  fallback_ctx <- score_setup(score, demo$data)
  fallback <- list(
    contribution = fallback_ctx$score,
    single = score_eval(score, subsets[1L, ], demo$data, fallback_ctx),
    batch = score_eval_batch(score, subsets, demo$data, fallback_ctx)
  )

  expect_equal(fallback, compiled, tolerance = 1e-12)
})
```

- [ ] **Step 2: Run the parity file**

```bash
Rscript --vanilla -e 'devtools::test(filter = "compiled-fallback-parity", reporter = "summary")'
```

Expected: zero failures; a skip is acceptable only if the compiled PLQR symbols are absent.

### Task 6: Verify the final snapshot and quantify the outcome

**Files:**
- Modify only if a verification failure identifies an in-scope defect.

**Interfaces:**
- Consumes: all changed test files and the package CI helpers.
- Produces: fresh local evidence for test execution, installed-layout behavior, package-check health, and post-change coverage.

- [ ] **Step 1: Run formatting and dependency guards**

```bash
git diff --check
Rscript tools/ci/lint.R TRUE
Rscript tools/ci/check_tests_declared_deps.R
```

Expected: all commands exit 0.

- [ ] **Step 2: Run the full source-tree suite once**

```bash
Rscript --vanilla -e 'devtools::test(reporter = "summary")'
```

Expected: zero failures and zero errors; optional-dependency skips must have explicit reasons.

- [ ] **Step 3: Run the suite against an installed package**

Install into a fresh temporary library, then run `testthat::test_dir(..., load_package = "installed")`. Confirm the repaired CATE method and custom-extension tests execute and that no namespace-layout skips remain.

- [ ] **Step 4: Run the local package check**

```bash
Rscript --vanilla -e 'res <- rcmdcheck::rcmdcheck(".", args = c("--no-manual", "--no-vignettes"), error_on = "never", check_dir = tempdir()); quit(status = as.integer(length(res$errors) > 0L || length(res$warnings) > 0L))'
```

Expected: zero errors and zero warnings. Record any pre-existing notes separately.

- [ ] **Step 5: Rerun coverage in an archived clean copy**

Create a temporary source copy from the final tracked tree plus the working test changes, run `covr::package_coverage()` there, and report overall and per-file changes. Do not run instrumented coverage with `clean = FALSE` in the live checkout because it leaves `.gcda`, `.gcno`, and `covr/` artifacts.

- [ ] **Step 6: Audit the remaining uncovered expressions**

Group remaining gaps into active public behavior, optional integrations, inactive/dead helpers, and low-value display/benchmark branches. Add another test only if this fresh audit reveals a material active contract omitted by Tasks 1-5.

- [ ] **Step 7: Inspect the final diff**

```bash
git status --short
git diff --stat
git diff -- tests/testthat markdown_notes/coverage
```

Expected: only the approved coverage note and focused test changes are present; no generated compiled artifacts are tracked or untracked.
