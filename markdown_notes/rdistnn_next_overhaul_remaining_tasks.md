# RdistNN Next Overhaul Note: Remaining Pre-Public Tasks Only

**Created:** 2026-07-09  
**Repository:** `JakobJuergens/RdistNN`  
**Scope:** Remaining work after the latest RdistNN cleanup commits  
**Audience:** author-facing and Codex-facing  
**Purpose:** This note lists only tasks still worth doing next. It intentionally omits items already completed by the recent cleanup commits.

---

## 0. Current assumptions

This note assumes the package has already adopted the paper-aligned PLQR public names (`g`, `g_col`, `theta_eval`, `theta_se`, `moment_se`, `theta_vcov`, `theta_omega`), has already removed the main old compatibility aliases, and has already pruned the most obvious raw distance/order/rank-table public exports.

The next pass should focus on remaining issues:

1. add a principled PLQR jackknife-meat sandwich route;
2. make PLQR supplied-`theta_eval` semantics impossible to misread;
3. fix PLQR plotting defaults and estimate-data extraction;
4. update the remaining stale documentation and design-note wording;
5. tighten direct-sandwich robustness and diagnostics;
6. decide what to do with remaining exported advanced APIs;
7. remove tracked generated CI logs/artifacts from the public branch.

---

## 1. Add a PLQR jackknife-meat sandwich route

### 1.1 Problem

The current PLQR implementation supports coefficient inference through a density-backed direct sandwich. It also computes `moment_se` as a jackknife diagnostic for localized PLQR moments when density is not supplied. What remains missing is a full jackknife covariance matrix for the localized vector moment that can be used as the **meat** of a coefficient sandwich.

Do not confuse these two objects:

- `moment_se`: componentwise jackknife standard errors for the moment; diagnostic only.
- `moment_vcov`: full covariance matrix of the localized vector moment; valid sandwich meat when paired with the density-backed bread.

### 1.2 Target formula

For target point `x`, supplied coefficient `theta_eval`, DNN weights `w_i(x)`, and PLQR score

\[
\psi_i(\theta)=\{\tau-1(Y_i\le D_i^\top\theta+g_i)\}\tilde D_i,
\qquad
\tilde D_i=D_i-m_D(U_i),
\]

the localized moment is

\[
\widehat m_n(\theta;x)=\sum_i w_i(x)\psi_i(\theta).
\]

The density-backed bread is

\[
\widehat J(x)=\sum_i w_i(x)\widehat f_i\tilde D_iD_i^\top.
\]

The current analytic meat route uses

\[
\widehat S_{analytic}(x)=\tau(1-\tau)\sum_i w_i(x)\tilde D_i\tilde D_i^\top.
\]

The jackknife meat route should compute

\[
\widehat S_{JK}(x)=\frac{n-1}{n}\sum_{\ell=1}^n
\{\widehat m_{n,-\ell}(\theta;x)-\bar m_{JK}(\theta;x)\}
\{\widehat m_{n,-\ell}(\theta;x)-\bar m_{JK}(\theta;x)\}^\top.
\]

Then coefficient covariance is

\[
\widehat V_\theta^{JK}(x)=\widehat J(x)^{-1}\widehat S_{JK}(x)\widehat J(x)^{-T}.
\]

### 1.3 API addition

Add an explicit coefficient-meat choice to `plqr_dnn()`:

```r
plqr_dnn(
  ...,
  density_col = "density_oracle",
  theta_eval = theta0,
  se = TRUE,
  theta_meat = c("analytic", "jackknife")
)
```

Default should remain:

```r
theta_meat = "analytic"
```

Use `theta_meat`, not a generic `meat`, because the option controls coefficient-level inference, not moment diagnostics.

### 1.4 Return-field design

For `theta_meat = "analytic"`:

```r
theta_se
theta_vcov
theta_omega
score_covariance      # analytic S
theta_meat_method = "analytic"
moment_se = NULL unless explicitly requested as a diagnostic
moment_vcov = NULL
```

For `theta_meat = "jackknife"`:

```r
moment_vcov           # list of p_d x p_d matrices, one per target
moment_se             # sqrt(diag(moment_vcov[[j]])), as p_d x m matrix
theta_se
theta_vcov
theta_omega
theta_meat_method = "jackknife"
moment_se_method = "jackknife"
```

For no density input:

```r
theta_se = NULL
theta_vcov = NULL
theta_omega = NULL
theta_meat_method = "none"
moment_se = optional diagnostic
moment_vcov = optional if implemented
```

### 1.5 Internal implementation sketch

Add an internal function that computes full vector covariance, not only componentwise SEs.

Possible names:

```r
.dnn_jackknife_vcov_loo_matrix()
.plqr_dnn_jackknife_vcov_loo()
```

Expected output:

```r
list(
  vcov = list(target_1 = matrix(p_d, p_d), ...),
  se = matrix(p_d, m)
)
```

The existing fast/loo jackknife paths should be reused where possible. The first implementation should prioritize transparent correctness and tests. Optimization can follow after the route is mathematically locked down.

### 1.6 Required tests

Use `p_d = 2` or larger so off-diagonal covariance is visible.

Minimum tests:

1. `moment_vcov` is a named list of length `m`, each entry `p_d x p_d`.
2. `moment_se[, j]^2` equals `diag(moment_vcov[[j]])`, up to the package’s jackknife scaling convention.
3. A small-n fixture manually computes leave-one-out moment vectors and verifies `moment_vcov`.
4. `theta_vcov[[j]]` equals `solve(J) %*% moment_vcov[[j]] %*% t(solve(J))` for nonsingular `J`.
5. A correlated-score fixture produces nonzero off-diagonal covariance.
6. `meta$theta_meat_method` or the chosen metadata field records the meat route.
7. The analytic and jackknife routes can both be requested without changing the localized moment itself.

---

## 2. Revise PLQR variance-design note

### 2.1 Problem

The current design note still states that coefficient-level inference is gated on the direct sandwich and that no theorem supports translating PLQR moment jackknife SEs into coefficient intervals. That statement is too broad unless it is restricted to **componentwise moment SEs**. A full jackknife covariance of the localized vector moment is a valid candidate meat, once sandwiched with the density-backed bread.

### 2.2 Required replacement language

Replace the overly broad wording with:

```markdown
A componentwise jackknife SE for the localized PLQR moment is not by itself a coefficient-level SE. Coefficient inference requires a bread/meat/bread sandwich. The bread remains density-backed in all coefficient-inference routes. The meat can be estimated either by the analytic PLQR formula S = tau(1-tau) sum_i w_i tilde_D_i tilde_D_i' or by a full vector jackknife covariance of the localized PLQR moment. The package must distinguish moment diagnostics from jackknife-meat coefficient inference.
```

### 2.3 Route table to add

Add a table to design note 15:

| Route | Bread | Meat | Coefficient inference? | Package fields |
|---|---|---|---|---|
| Analytic direct sandwich | density-backed `jacobian` | analytic `score_covariance` | yes | `theta_vcov`, `theta_se` |
| Jackknife-meat sandwich | density-backed `jacobian` | full `moment_vcov` | yes | `moment_vcov`, `theta_vcov`, `theta_se` |
| Moment diagnostic only | none | componentwise `moment_se` | no | `moment_se` |

### 2.4 Scope note update

The scope note should not say future contributors should avoid jackknife-based coefficient CI in all forms. It should say future contributors should avoid **componentwise moment-SE shortcuts**, but may implement a full jackknife-meat sandwich route with appropriate tests and theory-facing documentation.

---

## 3. Fix PLQR supplied-`theta_eval` semantics in S3 methods

### 3.1 Problem

`plqr_dnn()` currently records `meta$theta_source = "supplied"`, but `summary()`, `confint()`, `coef()`, and `predict()` still resemble ordinary coefficient-estimation methods. This can make `theta_eval` look like an estimated local root.

### 3.2 `summary.plqr_result()`

Add `theta_source` to the summary object and printed summary.

Suggested printed header:

```text
theta_source: supplied; no local PLQR root solving was performed
```

Consider suppressing `statistic` and `p.value` whenever `theta_source = "supplied"`.

Preferred behavior:

```r
if (identical(object$meta$theta_source, "supplied")) {
  statistic <- NA_real_
  p.value <- NA_real_
}
```

Rationale: `estimate / std.error` for a supplied evaluation value is easy to misread as a test of an estimated coefficient. Confidence intervals centered at supplied `theta_eval` are already delicate; p-values are more misleading.

If p-values are kept, add an unavoidable note in `summary()` output:

```text
statistics are centered at supplied theta_eval; this object does not solve a local PLQR root
```

### 3.3 `coef.plqr_result()`

Document that `coef()` returns the coefficient value represented by the object. For current `plqr_dnn()` fits, this is the supplied `theta_eval` repeated across targets.

Optionally attach:

```r
attr(out, "theta_source") <- object$meta$theta_source
```

### 3.4 `confint.plqr_result()`

Add a `theta_source` column or an attribute. The method should clearly document that intervals are centered at the supplied `theta_eval` unless a future root-solving mode produces `theta_source = "local_root"`.

### 3.5 `predict.plqr_result()`

Update docs to say explicitly:

```text
The returned fit matrix repeats the supplied theta_eval. Target-specific variation enters only through DNN weights and resulting direct-sandwich standard errors. No target-specific PLQR root is solved.
```

Add tests checking this language or checking an output attribute such as `theta_source = "supplied"`.

---

## 4. Fix PLQR plotting semantics

### 4.1 Problem

`plot_data_estimates.plqr_result()` currently prioritizes `object$theta` whenever present. Since current PLQR fits always store `theta`, estimate plots show supplied `theta_eval`, not localized moments. That conflicts with zero-reference line logic that is meaningful for PLQR moments, not for arbitrary supplied coefficients.

### 4.2 Required design change

Add a quantity/type argument to PLQR estimate extraction and plotting.

Suggested API:

```r
plot_data_estimates.plqr_result <- function(
  object,
  level = 0.95,
  quantity = c("moment", "theta"),
  ...
)
```

and similarly in plotting wrappers if needed:

```r
plot_estimates(fit, quantity = "moment")
plot_estimates(fit, quantity = "theta")
```

### 4.3 Default

Until root solving exists, the default for `plqr_result` should be:

```r
quantity = "moment"
```

### 4.4 Zero-reference rule

Use zero reference lines for:

- CATE effects;
- PLQR moments.

Do **not** automatically use zero reference lines for supplied `theta_eval` plots unless explicitly requested.

### 4.5 Candidate-grid behavior

Candidate-grid PLQR plots should default to moment coordinates because grid mode stores moment estimates. Do not force candidate-grid PLQR plots into supplied-theta space.

### 4.6 Tests

Add tests that:

1. PLQR `plot_data_estimates()` defaults to `moment`.
2. `quantity = "theta"` returns repeated `theta_eval`.
3. zero-reference lines appear for moment plots but not for theta plots.
4. candidate-grid PLQR plot data uses moment estimates by default.

---

## 5. Harden `.plqr_direct_sandwich()` internally

### 5.1 Problem

The public wrapper rejects direct coefficient inference with `normalize = FALSE`, but `.plqr_direct_sandwich()` itself still accepts arbitrary `W`. Since `predict.plqr_result()` and future call sites invoke this function directly, the internal function should assert normalized weights.

### 5.2 Required guard

Inside `.plqr_direct_sandwich()`, add:

```r
col_sums <- colSums(W)
if (any(abs(col_sums - 1) > 1e-8)) {
  .abort_bad_state(
    "Direct PLQR sandwich requires normalized weight columns."
  )
}
```

or define:

```r
.plqr_validate_normalized_weights(W, tol = 1e-8)
```

### 5.3 Tests

Add tests that rescale `W` or corrupt a stored model object and verify coefficient SE computation fails clearly rather than silently producing scale-dependent values.

---

## 6. Polish PLQR diagnostics

### 6.1 Rename density diagnostics

Current diagnostic names `density_min_weighted` and `density_max_weighted` are misleading if they are support extrema rather than weighted extrema. Rename:

| Current | Proposed |
|---|---|
| `density_min_weighted` | `density_min_support` |
| `density_max_weighted` | `density_max_support` |
| `density_median_weighted` | keep |

### 6.2 Add weighted density mean

Add:

```r
density_mean_weighted = sum(w * density)
```

per target.

### 6.3 Singular Jacobian tests

Add tests for:

1. singular `J` produces `theta_se = NA` and `theta_vcov_available = FALSE`;
2. diagnostics report `jacobian_available`, rank, and condition number correctly;
3. a nearly singular `J` produces a large/infinite condition number.

### 6.4 Vignette update

Update PLQR diagnostics examples to include the Jacobian diagnostics:

```r
d_plqr[, c(
  "target_id",
  "gram_min_eig",
  "gram_condition_number",
  "jacobian_min_singular_value",
  "jacobian_condition_number",
  "theta_vcov_available"
)]
```

---

## 7. Add true PLQR oracle density to demo data

### 7.1 Problem

The README currently sets a toy constant density column by hand. This is acceptable for a smoke example but is weak for public documentation.

### 7.2 Implementation

`plqr_demo_data()` uses asymmetric Laplace errors with known `tau` and `sigma`. Add the density at the zero quantile to the demo output.

Possible fields:

```r
data$density_oracle
truth$density_fn
```

If the density is constant under the current demo, document that.

### 7.3 Tests

Add tests that:

1. `density_oracle` is positive and finite;
2. `plqr_dnn(..., density_col = "density_oracle")` produces finite `theta_se` in a standard demo;
3. README/vignette examples no longer manually create `density_oracle`.

---

## 8. Clean remaining documentation wording

### 8.1 `DESCRIPTION`

The description still says variance estimation is via delete-one jackknife. Replace with estimator-specific wording:

```text
Uses estimator-specific inference routes, including delete-one jackknife diagnostics and direct sandwich covariance estimates.
```

### 8.2 `NEWS.md`

The older PLQR section still mentions `se`, `vcov`, and `omega` as coefficient fields, and the plotting section still mentions a `beta/tau` subtitle. Update or remove stale historical wording before public release.

Because this is pre-public NEWS, edit history for clarity rather than preserving outdated names.

### 8.3 Conditional Z vignette

Update text that still says coefficient inference is only the direct analytic route. Once jackknife meat is implemented, describe both coefficient-meat routes. If jackknife meat is not yet implemented, clarify that a full jackknife-meat route is a planned coefficient route distinct from componentwise `moment_se`.

### 8.4 README quick-start examples

Use `tuning_dnn_paper()` consistently in NP, CATE, and PLQR quick-start examples, or label explicit small `s` as toy. At present, PLQR uses paper tuning while NP/CATE quick-start examples still use toy fixed `s` values.

---

## 9. Decide remaining public API tiers

### 9.1 Remaining exported advanced APIs

The most obvious raw distance/order/rank-table exports have been pruned. The namespace still exports several advanced or experimental APIs:

- score constructors and score operations: `new_score()`, `register_score()`, `score_np_y()`, `score_cate_rho()`, `score_plqr_component()`, `score_diff_path()`, `score_eval()`, `score_eval_batch()`, `score_order()`, `score_setup()`;
- jackknife primitives: `jackknife_fast_rank()`, `jackknife_slow_rank()`, `jackknife_rank_full()`;
- TDNN user-facing APIs: `np_tdnn()`, `weights_tdnn()`, `weights_tdnn_exact()`, `weights_tdnn_asymp()`, `rank_weight_tdnn_exact()`, `rank_weight_tdnn_asymp()`, `plot_tdnn_scales()`;
- experimental IV/GMM APIs: `gmm_dnn()`, `gmm_model_dr_iv()`, `sim_iv_*`, `summarise_iv_results()`, `iv_demo_data()`.

### 9.2 Required decision

Create an explicit public API tier table:

| Tier | Meaning | Examples | Public docs location |
|---|---|---|---|
| Stable user API | expected to survive public release | `weights_dnn()`, `cate_dnn()`, `plqr_dnn()` | README + main reference |
| Developer API | intended for custom extensions | distance/rank-weight/score registries | advanced article |
| Paper replication API | tied to simulations | `sim_cate_*` | replication article |
| Experimental companion API | may change | GMM/IV, TDNN if not central | experimental article |
| Internal | not exported | low-level helpers | no public docs |

### 9.3 Package documentation change

Do not leave all exported functions in a flat reference index. Move advanced score/jackknife/TDNN/GMM functions into clearly labelled pkgdown groups.

---

## 10. Reposition experimental IV-GMM in README

### 10.1 Problem

The README labels IV-GMM as experimental, but still gives it a full quick-start block. That can make it look comparable in maturity to CATE/PLQR.

### 10.2 Recommended change

Move the full IV-GMM example out of the README quick start into a dedicated experimental vignette, for example:

```r
vignette("experimental-gmm", package = "RdistNN")
```

In README, keep only a short note:

```text
Experimental local IV-GMM wrappers are included for companion-paper development. Their API may change before the GMM paper is released. See the experimental GMM vignette for details.
```

---

## 11. Remove tracked generated CI logs/artifacts from the public branch

### 11.1 Problem

Current HEAD is a `GH_logs` update commit marked `[skip ci]`. Generated CI logs should not be part of the public source history going forward. They create noisy diffs and obscure substantive changes.

### 11.2 Required cleanup

1. Stop tracking `GH_logs/` in git.
2. Add to `.gitignore`:

```gitignore
GH_logs/
archives/*.zip
archives/*.tar.gz
```

3. Keep CI logs in GitHub Actions artifacts, not the source tree.
4. If logs are needed for debugging, store them outside the public branch or in release artifacts.

### 11.3 Release rule

No public-release candidate should have a latest commit that only updates generated logs.

---

## 12. CI and release hygiene

### 12.1 Required checks after next pass

Run and record:

```r
devtools::document()
devtools::test()
devtools::check()
```

Then trigger GitHub Actions manually on the final code commit, not on a log-update commit.

### 12.2 Acceptance criteria

- `devtools::check()` has 0 errors and 0 warnings.
- GitHub Actions are attached to the public-release candidate commit.
- Generated docs are in sync.
- No tracked generated CI logs change in the release commit.
- README examples use current API names and do not create misleading inference examples.

---

## 13. Cross-repository compatibility sweep

Before finalizing the package API, search the paper repositories for old names and changed columns.

Search terms:

```text
q_oracle
q_col
q =
beta_true
beta =
plqr_moment
return = "g"
fit$se
fit$vcov
fit$omega
kappa
kappa_grid
tau_x0
truth$tau_fn
tau_oracle
register_rank
get_rank
rank_names
```

Update package-consuming code in:

- `Cond_Z_Estimation`;
- `Cond_GMM_with_DNN`;
- any simulation scripts or analysis notebooks;
- any website/software examples.

For historical simulation outputs with old column names, prefer explicit one-time migration or reader-side translation over reintroducing package-level aliases.

---

## 14. Suggested Codex prompt for the next overhaul

```text
Update RdistNN’s next pre-public overhaul note and implement the remaining tasks only.

1. Add PLQR jackknife-meat coefficient inference:
   - Add theta_meat = c("analytic", "jackknife") to plqr_dnn(), default analytic.
   - Implement full vector jackknife covariance for the localized PLQR moment, not only componentwise moment SEs.
   - Store moment_vcov and use theta_vcov = J^{-1} moment_vcov J^{-T} for theta_meat = "jackknife".
   - Preserve the existing analytic meat route.
   - Add small-n manual covariance tests, off-diagonal covariance tests, and sandwich propagation tests.

2. Revise design note 15:
   - Replace the too-broad “jackknife not valid for PLQR coefficient inference” framing.
   - Distinguish componentwise moment_se from full moment_vcov used as sandwich meat.
   - Add a route table: analytic meat, jackknife meat, moment diagnostic only.

3. Fix PLQR supplied-theta semantics:
   - summary.plqr_result() must show theta_source.
   - Suppress statistic/p.value for theta_source = "supplied", or add an unavoidable warning note.
   - coef(), confint(), and predict() documentation must say current PLQR objects repeat supplied theta_eval and do not solve local roots.

4. Fix PLQR plotting:
   - Add quantity/type = c("moment", "theta") to PLQR plot data or plotting methods.
   - Default PLQR plots to moment while no root solver exists.
   - Use zero reference lines for PLQR moments, not supplied theta plots.

5. Harden direct sandwich:
   - Add an internal normalized-weight check inside .plqr_direct_sandwich().
   - Add tests for rescaled/corrupted weights.

6. Polish diagnostics:
   - Rename density_min_weighted/density_max_weighted to support-extrema names.
   - Add density_mean_weighted.
   - Add singular/near-singular Jacobian tests.

7. Add real PLQR oracle density:
   - Add density_oracle to plqr_demo_data() using the ALD DGP.
   - Update README/vignettes to use it rather than manually setting density_oracle <- 2.

8. Clean docs and release hygiene:
   - Update DESCRIPTION’s variance sentence.
   - Update stale NEWS references to se/vcov/omega and beta/tau subtitle.
   - Move full IV-GMM example out of README quick start or shorten it to an experimental note.
   - Stop tracking GH_logs/ and generated archives; update .gitignore.
   - Run devtools::document(), devtools::test(), devtools::check(), then trigger GitHub Actions on the final code commit.
```

---

## 15. Final next-step summary

The next overhaul should not redo completed renaming and namespace-pruning work. The remaining high-value changes are:

1. implement full jackknife-meat PLQR sandwich inference;
2. revise design note 15 accordingly;
3. make supplied-`theta_eval` semantics explicit in all PLQR S3 methods;
4. make PLQR plots default to moments rather than supplied coefficients;
5. harden direct sandwich internals;
6. polish density/Jacobian diagnostics;
7. add real oracle density to PLQR demo data;
8. finish documentation cleanup;
9. tier the remaining advanced public APIs;
10. remove tracked generated CI logs/artifacts before public release.
