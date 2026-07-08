# RdistNN Public-Release Cleanup: Legacy Compatibility, Naming, and Paper-Notation Alignment

**Created:** 2026-07-08  
**Repository:** `JakobJuergens/RdistNN`  
**Scope:** Pre-public API cleanup before broader release / website exposure  
**Primary goal:** Remove avoidable legacy compatibility, align package names and symbols with the conditional Z-estimation paper, and reduce confusion before users begin relying on public APIs.

---

## 0. Executive summary

The package is close to being public-useful, but it is not yet public-clean. The current codebase reflects several development histories layered on top of each other:

1. the generalized U-statistics / jackknife project;
2. the early DNN wrapper implementation;
3. the current conditional Z-estimation paper;
4. the emerging conditional GMM and PLQR companion functionality.

This is understandable internally. It is risky publicly. Once the package is linked from a website, arXiv paper, or documentation site, exported names and examples become sticky. Before that happens, the package should undergo a deliberate breaking cleanup pass.

The main issue is not computational correctness. The recent PLQR work has materially improved validation, documentation, tests, and S3 method parity. The issue is public semantics: several exported names, result fields, and examples still encode legacy conventions that now conflict with the paper’s notation. The largest offenders are the PLQR API: `beta` versus `theta`, `q` / `q_col` / `q_oracle` versus the paper’s nonparametric shape nuisance \(g_0(U)\), overloaded `se`, public helper functions tied to demo-only column names, and a PLQR public surface that reads closer to coefficient estimation than the current supplied-coefficient implementation actually supports.

The package has essentially no external compatibility burden yet. That is an advantage. Use it. Do not ship public version 0.1.0 with deprecated aliases and “legacy compatibility” wrappers unless there is a strong reason. The first public version should reflect the paper’s notation, not the package’s development history.

---

## 1. Guiding principles for the cleanup

### 1.1 Prefer breaking cleanup now over deprecation later

Because the package is not yet public, compatibility with prior private usage should not dominate API design. If a name is wrong, rename it. If an exported function is not meant for users, unexport it. If a field changes meaning across cases, split it.

A compatibility layer makes sense only when:

- an external paper, vignette, or public code already relies on the old name;
- the old name is genuinely clearer for users;
- or removing it would create excessive internal churn relative to the benefit.

Most current legacy issues do not satisfy these criteria.

### 1.2 Public names should match the paper’s notation

The package is intended to support the conditional Z-estimation paper. Public names should therefore mirror the paper as much as possible:

| Concept | Paper notation | Preferred package name |
|---|---:|---|
| localization covariates | \(X_i\) | `X`, `x_cols` |
| target point(s) | \(x\), \(x_0\) | `x0`, `X0` |
| DNN order | \(s\) | `s` |
| fixed nearest-neighbor count | \(k\) | `k` |
| PLQR coefficient | \(\theta\) | `theta`, `theta_eval`, `theta_start` |
| PLQR shape nuisance | \(g_0(U)\) | `g`, `g_col`, `g_oracle` |
| residualizing nuisance | \(m_0(U)=E[D\mid U]\) | `mD`, `mD_cols` |
| PLQR density at threshold | \(q_{0,x}(D,Z)\) or density object | `density`, preferably `density_at_threshold` |
| observation-level score | \(\psi(W_i;\theta,\eta)\) | `score`, `psi` |
| localized moment | \(\sum_i w_i\psi_i\) | `moment` |
| coefficient covariance | \(\widehat V_\theta\) | `theta_vcov` |

Where the package deviates, it should be because the package name is unambiguously better, not because of historical artifact.

### 1.3 Distinguish estimator objects from diagnostic objects

The current CATE wrapper estimates a local scalar target from a pseudo-outcome. The current PLQR wrapper does **not** yet solve the local PLQR Z-estimation problem. It evaluates a nonsmooth PLQR moment at a supplied coefficient value and, when a density is supplied, attaches the direct coefficient sandwich around that value.

This is a valid object, but it is not the same maturity level as CATE in the econometric sense. Public names should not imply otherwise.

### 1.4 Keep the stable public API small

Before public release, every exported function should be deliberately classified:

- stable public API;
- advanced developer API;
- paper-replication API;
- experimental API;
- internal, not exported.

The current namespace is much too broad. Public release should not expose every internal layer.

---

## 2. Package identity and metadata

### 2.1 Current issue

The package metadata still describes the package as:

> “Generalized U-Statistics and Conditional Z-Estimation with Pluggable Scores”

and the description foregrounds generalized U-statistics, pluggable scores, local nonparametric regression, CATE, PLQR, and delete-one jackknife variance estimation.

This is too broad for the package’s intended public role. It makes the package look like a general research workbench rather than a focused companion to the conditional Z-estimation paper.

### 2.2 Why this matters

A user arriving from the conditional Z-estimation paper should immediately understand the package’s role:

- implements DNN localization;
- supports local conditional moment / Z-estimation workflows;
- provides CATE and PLQR wrappers;
- includes diagnostics and replication helpers.

If the package leads with generalized U-statistics, the public story becomes diffuse. It also risks mixing the job-market paper with the jackknife paper in a way that weakens the website/software narrative.

### 2.3 Recommended change

Keep the package name `RdistNN` if changing it would be too disruptive. But change the package title and description.

Suggested title:

```text
Distributional Nearest-Neighbor Inference for Local Conditional Moment Models
```

Alternative:

```text
Distributional Nearest-Neighbor Methods for Local Z-Estimation
```

Suggested description:

```text
Implements distributional nearest-neighbor weights and diagnostics for local
conditional moment and Z-estimation problems. Includes user-facing wrappers for
localized conditional means, conditional average treatment effects, and partially
linear quantile-regression moment evaluation, with tuning-grid diagnostics and
paper-replication tools. Lower-level generalized U-statistic and jackknife
utilities are included as supporting infrastructure.
```

The generalized U-statistic layer should be mentioned, but not as the main identity.

---

## 3. README and public documentation

### 3.1 Current issue: README PLQR example is stale

The README’s PLQR quick-start example fits `plqr_dnn()` without supplying `density` / `density_col`, then prints `fit_plqr$moment`, `fit_plqr$se`, and diagnostics. In the current implementation, this means `fit_plqr$se` is a moment-level jackknife diagnostic, not coefficient-level inference.

The newer vignette correctly explains that coefficient-level PLQR SEs and confidence intervals require the direct sandwich route with a supplied density. The README should match that standard.

### 3.2 Why this matters

The README is the first thing many users will read. If it shows `fit_plqr$se` without density, users may infer that PLQR coefficient inference works like CATE jackknife inference. That is exactly the misunderstanding the new design note is trying to prevent.

### 3.3 Recommended change

Rewrite the PLQR README example to include density and coefficient-level methods:

```r
demo_plqr <- plqr_demo_data(n = 250, p = 2, p_d = 2, tau = 0.5, seed = 1)
demo_plqr$data$density_oracle <- 2  # toy/oracle density for example only

fit_plqr <- plqr_dnn(
  data = demo_plqr$data,
  x_cols = c("X1", "X2"),
  D_cols = c("D1", "D2"),
  g_col = "g_oracle",
  mD_cols = c("mD1", "mD2"),
  y_col = "Y",
  density_col = "density_oracle",
  theta_eval = demo_plqr$truth$theta_true,
  tau = 0.5,
  x0 = c(0, 0),
  method = "exact",
  k = 1L,
  tuning = tuning_dnn_paper(),
  se = TRUE,
  return = c("W", "score", "ess")
)

summary(fit_plqr)
confint(fit_plqr)
predict(fit_plqr, x0 = c(0.25, -0.25), interval = "confidence")
diagnostics(fit_plqr)
```

If density is omitted in any example, label the example explicitly as a **moment diagnostic** example.

### 3.4 Current issue: examples use toy `s = 3L`

The README examples often use small explicit `s` values, such as `s = 3L`. That is acceptable for toy examples, but the package now has a paper-tuning convention via `tuning_dnn_paper()`.

### 3.5 Recommended change

Use paper-tuning examples in the README:

```r
tuning = tuning_dnn_paper()
```

If small fixed `s` values are retained for speed, label them:

```r
# Toy choice for a tiny example; use tuning_dnn_paper() in paper-style runs.
s = 3L
```

---

## 4. PLQR API: the largest cleanup target

### 4.1 Current issue: `beta` input versus `theta` output

`plqr_dnn()` currently accepts a `beta` argument and returns `theta`. Internally, the output `theta` is just the supplied `beta` repeated across targets.

This creates three problems:

1. The paper’s coefficient notation is \(\theta\), not \(\beta\).
2. Users may infer that `theta` was estimated.
3. Prediction appears to produce a target-specific coefficient surface, when it currently repeats the supplied coefficient vector.

### 4.2 Recommended change

Rename the input before public release.

Preferred option while no root solver exists:

```r
theta_eval
```

Use:

```r
plqr_dnn(..., theta_eval = ...)
```

Return:

```r
fit$theta
fit$theta_source = "supplied"
```

or:

```r
fit$theta_eval
```

If you expect root solving soon, design now for both modes:

```r
plqr_dnn(
  ...,
  theta = NULL,
  theta_start = NULL,
  solve = FALSE
)
```

But if no solver exists yet, the cleanest immediate public name is `theta_eval`.

### 4.3 Implementation consequences

Rename throughout:

| Current | Proposed |
|---|---|
| argument `beta` | `theta_eval` |
| internal `beta` | `theta_eval` |
| `fit$beta` | remove or rename `fit$theta_eval` |
| `fit$theta` | keep, but document as supplied/evaluated |
| demo `beta_true` | `theta_true` |
| `plqr_moment(..., beta=)` | `plqr_score_matrix(..., theta=)` |

### 4.4 Do not keep `beta` as a public alias

Since the package is not public yet, do not keep `beta` as a deprecated alias. It will only preserve confusion.

---

## 5. PLQR nuisance naming: remove `q`

### 5.1 Current issue

The package still uses `q`, `q_col`, `q_oracle`, `q_fn`, and `q_x0`. The new design note correctly says these are legacy labels and that the object is really the nonparametric shape nuisance \(g_0(U)\), not a conditional quantile function of \(X\) alone.

This is the most important notation cleanup.

### 5.2 Why this matters

The paper’s PLQR setup has a partially linear threshold:

\[
D_i^\top \theta_0 + g_0(U_i).
\]

The package’s `q` object is this \(g_0(U_i)\). Calling it `q_tau(X)` is misleading because:

- it suggests a quantile function rather than a shape nuisance;
- it suggests dependence only on \(X\), while the paper permits \(U=(X,Z)\);
- it conflicts with the package’s GMM use of `q` for the moment/instrument dimension.

### 5.3 Recommended rename

| Current | Proposed |
|---|---|
| `q_col` | `g_col` |
| `q` | `g` |
| `q_oracle` | `g_oracle` |
| `q_fn` | `g_fn` |
| `q_x0` | `g_x0` |
| `q_density` | `density_oracle`, `f_oracle`, or `density_at_threshold_oracle` |

Update:

- `plqr_dnn()`;
- `plqr_demo_data()`;
- `plqr_moment()` or replacement helper;
- `score_plqr_component()`;
- tests;
- README;
- vignettes;
- roxygen documentation;
- pkgdown examples.

### 5.4 Do not keep `q_col` as alias

This is a pre-public breaking change. Keep the public API clean.

---

## 6. PLQR demo data

### 6.1 Current issue

`plqr_demo_data()` still describes the model as

\[
Q_\tau(Y\mid X,D)=D^\top\beta_0(\tau)+q_\tau(X),
\]

and returns `q_oracle`, `q_fn`, and `q_x0`.

This conflicts with the corrected theory mapping.

### 6.2 Recommended change

Rewrite the demo as:

\[
Q_\tau(Y\mid D,U)=D^\top\theta_0+g_0(U),
\]

with the demo specializing to \(U=X\).

Suggested output fields:

```r
list(
  data = df,               # includes X, D, Y, g_oracle, mD*, rD*
  x0 = x0,
  truth = list(
    theta_true = theta_true,
    g_fn = g_fn,
    mD_fn = mD_fn
  ),
  g_x0 = g_fn(matrix(x0, nrow = 1))
)
```

### 6.3 Column names

Use:

- `g_oracle`;
- `density_oracle` or `f_oracle` if density is added;
- keep `mD1`, `mD2`, etc. unless you want to move to `m_D1` style, which is less R-friendly.

---

## 7. PLQR moment helper

### 7.1 Current issue

`plqr_moment()` is exported and hardcodes the demo-data naming scheme:

- `D1`, `D2`, ...
- `mD1`, `mD2`, ...
- `q_oracle`
- `Y`

It is therefore not a clean user-facing helper. It is a demo helper.

### 7.2 Recommended change

Either remove it from the public namespace or replace it with a general helper.

Recommended replacement:

```r
plqr_score_matrix <- function(
  data,
  y_col,
  D_cols,
  g_col,
  mD_cols,
  theta,
  tau
)
```

Return an \(n\times p_D\) matrix of observation-level scores.

If you keep a convenience helper for demo data, make it internal or name it explicitly:

```r
.plqr_demo_score_matrix()
```

---

## 8. `score_plqr_component()`

### 8.1 Current issue: hardcoded `Y`

`score_plqr_component()` accepts `q_col`, `D_cols`, `mD_cols`, etc., but hardcodes the outcome column as `"Y"`.

This is inconsistent with `plqr_dnn()`, which now correctly requires `y_col` in the data-frame interface.

### 8.2 Recommended change

If the function remains public, add:

```r
y_col
```

and use it inside setup.

### 8.3 Current issue: ignored `k`

The function documents `k` as “kept for API compatibility with other scores; ignored here.”

That is a pure compatibility artifact. It should not be public.

### 8.4 Recommended change

Remove `k` from `score_plqr_component()` if the function remains public.

If some internal score API requires `k`, make a private internal wrapper and keep the public constructor clean.

### 8.5 Current issue: old notation

The score constructor still uses `q_col`, \(q_\tau(X)\), `m_D(X)`, and error messages involving `q_tau(X)`.

Rename consistently:

| Current | Proposed |
|---|---|
| `q_col` | `g_col` |
| `q_tau(X)` in errors | `g` or `g_col` |
| `beta` | `theta` |
| `g` as computed score | `score` or `psi` |

### 8.6 Current issue: malformed roxygen

There is a malformed roxygen tag:

```r
#' @examplesIf TRUE
```

Fix this before public documentation builds.

### 8.7 Export status

If the score registry is not central to the public package, unexport `score_plqr_component()` entirely. If it remains exported, make it clean, general, and notation-aligned.

---

## 9. PLQR result fields and standard errors

### 9.1 Current issue: `se` is overloaded

For PLQR, `se` can mean two different things:

1. coefficient-level direct-sandwich SE when density is supplied;
2. moment-level jackknife SE when density is not supplied.

This is dangerous. A user may inspect `fit$se` without checking `meta$se_route` and misinterpret it.

### 9.2 Recommended split

Use explicit fields:

| Field | Meaning |
|---|---|
| `theta_se` | coefficient-level direct-sandwich SE |
| `theta_vcov` | coefficient-level covariance matrix/list |
| `theta_omega` | asymptotic covariance \(\Omega_\theta\) before dividing by \(n_{\mathrm{eff}}\) |
| `moment` | localized score/moment |
| `moment_se` | moment-level jackknife diagnostic |
| `score` or `psi` | observation-level score matrix if requested |
| `score_covariance` | \(S\) object in the sandwich |

Avoid storing PLQR `se` unless it has a single fixed meaning.

### 9.3 Method behavior

Methods should use `theta_se` and `theta_vcov` only:

- `summary.plqr_result()` reports coefficient table only when `theta_se` exists.
- `confint.plqr_result()` requires `theta_se` or `theta_vcov`.
- `predict.plqr_result(se.fit=TRUE)` requires direct sandwich.

Moment SEs should be shown in diagnostics, not in coefficient summaries.

---

## 10. `se_route` and inference metadata

### 10.1 Current issue

`meta$se_route` is a single string with values like:

- `"direct_sandwich"`;
- `"moment_jackknife"`;
- `"none"`.

This is ambiguous because PLQR has two distinct objects:

- coefficient inference;
- moment diagnostics.

### 10.2 Recommended replacement

Use structured metadata:

```r
meta$theta_inference <- list(
  method = "direct_sandwich",
  available = TRUE,
  requires_density = TRUE
)

meta$moment_inference <- list(
  method = "jackknife",
  available = TRUE
)
```

or simpler:

```r
meta$theta_se_method = "direct_sandwich"
meta$moment_se_method = "jackknife"
```

Do not use a single `se_route` to describe both.

---

## 11. Prediction semantics

### 11.1 Current issue

`predict.plqr_result()` currently validates stored model state, computes weights at new target points, and returns the supplied coefficient vector repeated across targets with optional target-specific direct-sandwich SEs.

This is technically consistent with the current implementation. But public users will interpret `predict()` as producing estimated local coefficients at new targets.

### 11.2 Recommended change before root solving

Either:

1. keep `predict.plqr_result()` but document very sharply that it evaluates the stored supplied coefficient at new targets; or
2. do not advertise `predict.plqr_result()` in the README front page until root solving exists.

A better method name internally would be:

```r
evaluate.plqr_result()
```

but S3 `predict()` is useful once root solving exists.

### 11.3 Add source flag

Add:

```r
meta$theta_source = "supplied"
```

For future root-solving fits:

```r
meta$theta_source = "local_root"
```

Then `print.plqr_result()` can say:

```text
theta : supplied evaluation value
```

rather than implying estimation.

---

## 12. Direct sandwich normalization

### 12.1 Current issue

The direct sandwich uses the weight matrix as supplied and computes:

```r
neff <- 1 / colSums(W^2)
```

If `normalize = FALSE`, coefficient SEs can depend on arbitrary weight scaling.

### 12.2 Recommended change

For direct PLQR coefficient inference, require normalized weights:

```r
if (!normalize && se && !is.null(density_vec)) {
  abort("Direct PLQR coefficient inference requires normalized weights.")
}
```

Alternatively, internally normalize each weight column inside `.plqr_direct_sandwich()` before computing \(S\), \(J\), and \(n_{\mathrm{eff}}\).

I prefer the guard because it makes the theorem convention explicit.

---

## 13. Direct sandwich Jacobian diagnostics

### 13.1 Current issue

PLQR diagnostics currently include residualized Gram diagnostics based on \(D-m_D\). This is useful, but coefficient inference depends on the density-weighted Jacobian:

\[
J(x)=\sum_i w_i f_i \tilde D_iD_i'.
\]

A target may have acceptable residualized Gram behavior but an ill-conditioned \(J(x)\).

### 13.2 Recommended additions

When direct sandwich inference is available, diagnostics should include:

- `jacobian_available`;
- `jacobian_min_singular_value`;
- `jacobian_condition_number`;
- `jacobian_rank`;
- `density_min_weighted`;
- `density_median_weighted`;
- `density_max_weighted`;
- `theta_vcov_available`.

This should be target-level.

### 13.3 Naming

If the stored `jacobian` uses the positive-bread convention rather than the literal derivative sign, document this. The sign does not affect \(J^{-1}SJ^{-T}\), but it matters for future root solving.

---

## 14. PLQR grid mode

### 14.1 Current issue

Direct PLQR sandwich inference is unavailable in `tuning_mode = "grid"`. Density with grid mode currently errors. This is acceptable for the current stage, but it limits PLQR/CATE parity.

### 14.2 Recommendation

Do not implement grid-mode coefficient inference before root solving unless needed. But document the limitation as intentional:

- grid mode currently evaluates moment sensitivity;
- coefficient-level direct sandwich is default-mode only;
- future work: candidate-level direct sandwich after local root solving.

### 14.3 Once root solving exists

Extend grid payload with candidate-level:

- `theta`;
- `theta_se`;
- `theta_vcov` summary;
- `jacobian_condition_number`;
- `moment_norm`;
- convergence status.

---

## 15. CATE versus PLQR maturity

### 15.1 Current state

CATE estimates a target directly from a pseudo-outcome and uses jackknife SEs. PLQR evaluates a nonsmooth moment at a supplied coefficient and uses direct sandwich coefficient inference only when density is supplied.

### 15.2 Public communication

Do not say PLQR is at CATE maturity in the estimator sense. It is closer in:

- API surface;
- tests;
- diagnostics;
- documentation;
- method consistency.

It is not yet equivalent in:

- target solving;
- simulation support;
- root-finding diagnostics;
- candidate-level coefficient inference.

### 15.3 Suggested public phrase

```text
PLQR support currently evaluates localized orthogonal PLQR moments and
direct-sandwich coefficient uncertainty at supplied coefficient values.
Local root-solving for target-specific PLQR coefficients is under development.
```

Use this until root solving is implemented.

---

## 16. Rank-weight API naming

### 16.1 Current issue

The object is called a `rank_weight`, but registry functions are:

- `register_rank()`;
- `get_rank()`;
- `rank_names()`.

These names suggest ranking algorithms, not rank-weight functions.

### 16.2 Recommended rename

Before public release:

| Current | Proposed |
|---|---|
| `register_rank()` | `register_rank_weight()` |
| `get_rank()` | `get_rank_weight()` |
| `rank_names()` | `rank_weight_names()` |

Do not keep aliases unless required internally. Update examples, tests, and documentation.

### 16.3 Related functions

The evaluator names are already good:

- `rank_weight_eval()`;
- `rank_weight_eval_batch()`;
- `rank_weight_dnn_exact()`;
- `rank_weight_dnn_asymp()`.

The registry should match this family.

---

## 17. Distance API naming

### 17.1 Current issue

The package exports both:

- raw calculators: `dist_l2()`, `dist_l1()`, `dist_l2_scaled()`;
- spec constructors: `distance_l2()`, `distance_l1()`, etc.

This creates two public naming families for distances.

### 17.2 Recommendation

Make `distance_*()` the main public API. Consider unexporting raw `dist_*()` functions, or label them clearly as low-level helpers.

If raw calculators stay public, documentation should say:

- use `distance_l2()` in estimator functions;
- use `dist_l2()` only for direct numeric distance calculation.

### 17.3 Scaled distances

The registry names include entries like `l1_scaled_sd`, while spec names use strings like `l2_scaled[sd]`. This is not fatal, but it is slightly messy. Public examples should use one convention consistently.

---

## 18. Low-level exports

### 18.1 Current issue

The namespace exports many low-level functions:

- low-level DNN weight constructors from ranks;
- TDNN rank tables;
- jackknife internals;
- score constructors;
- simulation functions;
- IV/GMM helpers.

This is too much for a first public release.

### 18.2 Classification proposal

#### Stable public API

- `weights_dnn()`
- `np_dnn()`
- `cate_dnn()`
- `plqr_dnn()` after cleanup
- `diagnostics()`
- `summary()`, `coef()`, `confint()`, `predict()` methods
- `tuning_dnn_paper()`
- `distance_l2()`, `distance_l1()`, perhaps `distance_mahalanobis()`

#### Developer API

- `new_distance()`
- `as_distance()`
- `register_distance()`
- `new_rank_weight()`
- `as_rank_weight()`
- `register_rank_weight()`
- possibly `new_score()` if custom scores are part of the public story

#### Paper replication API

- `sim_cate_*`
- `summarise_sim_results()`

#### Experimental API

- `gmm_dnn()`
- `gmm_model_dr_iv()`
- `sim_iv_*`
- `summarise_iv_results()`

#### Internal / not exported

- `dnn_weights_*_from_ranks`
- `tdnn_rank_table_*`
- `jackknife_*`
- `score_*` constructors unless explicitly part of developer API
- raw compiled-path wrappers
- demo-only helpers

### 18.3 Public documentation structure

Pkgdown reference should reflect this classification. Do not present everything as equally public.

---

## 19. GMM / IV surface

### 19.1 Current issue

The README advertises local IV-GMM estimation and simulation. The GMM companion paper is not yet public-ready. If the package goes public first, this creates a mismatch: users may interpret the GMM layer as theorem-complete.

### 19.2 Recommended change

For the first public release, either:

1. remove GMM/IV from the README front page; or
2. label it explicitly as experimental companion functionality.

Example wording:

```text
Experimental companion functionality: local IV-GMM wrappers and simulation
helpers are included for ongoing companion-paper development. Their API and
theorem alignment may change before the GMM paper is released.
```

### 19.3 Do not delete code unnecessarily

You can keep GMM code in the repository. The issue is public presentation and export stability, not necessarily code presence.

---

## 20. TDNN / jackknife support

### 20.1 Current issue

The package exports TDNN and jackknife-related functions that are more directly tied to the jackknife paper than the conditional Z paper.

### 20.2 Recommendation

If `RdistNN` is meant to be a joint companion to both papers, keep them but document the split. If it is meant mainly for conditional Z, hide TDNN/jackknife details from the front page and main reference groups.

Possible reference structure:

- “DNN localization”
- “Conditional Z wrappers”
- “Diagnostics”
- “Paper replication”
- “Advanced: rank-weight and score APIs”
- “Advanced: TDNN / jackknife infrastructure”

---

## 21. Observation-level score naming

### 21.1 Current issue

`plqr_dnn()` uses return extra `"g"` and stores `out$g`. This conflicts with the proposed use of `g` for the PLQR shape nuisance \(g_0(U)\).

### 21.2 Recommended change

Rename:

| Current | Proposed |
|---|---|
| `return = "g"` | `return = "score"` or `return = "psi"` |
| `out$g` | `out$score` or `out$psi` |
| `g_mat` internally | `score_mat` or `psi_mat` |

I prefer `score` for user-facing code and `psi` internally if the paper uses \(\psi\).

---

## 22. Moment versus score convention

### 22.1 Recommended convention

Use the following consistently:

| Object | Meaning | Name |
|---|---|---|
| observation-level contribution | \(\psi_i(\theta,\eta)\) | `score` or `psi` |
| weighted local average | \(\sum_i w_i(x)\psi_i\) | `moment` |
| derivative / bread | \(J(x)\) | `jacobian` or `bread` |
| score covariance | \(S(x)\) | `score_covariance` |
| coefficient covariance | \(J^{-1}SJ^{-T}\) | `theta_omega` |
| finite-sample covariance | \(\theta_\omega/n_{\mathrm{eff}}\) | `theta_vcov` |

This will make the package easier to explain in the paper and documentation.

---

## 23. Compiled helper names

### 23.1 Current issue

Compiled helpers use names like:

- `rustat_plqr_ctx_create`;
- `rustat_plqr_g_all`;
- `rustat_subset_row_sums_g`.

These are internal, but the prefix `rustat` reflects older branding.

### 23.2 Recommendation

Not urgent unless exposed. Before public release, confirm:

- no user-facing docs mention `rustat_*`;
- no error messages expose `rustat_*`;
- generated docs do not include these helpers.

If you later clean compiled code, rename the internal C++ symbols too, but this is lower priority than public R API cleanup.

---

## 24. Generated artifacts and repository hygiene

### 24.1 Current issue

The repository has historically included or compared against large generated / debug artifacts such as `GH_logs` and `archives/source.zip`.

The current `.gitignore` ignores only `GH_logs/run_*.log`, not all `GH_logs` directories. It also does not ignore `archives/source.zip`.

### 24.2 Why this matters

Public repositories should not contain CI debug snapshots or generated source archives unless they are part of an explicit release process. They make diffs noisy and reviews harder.

### 24.3 Recommendation

Before public release:

- remove tracked debug logs from the public branch;
- ignore `GH_logs/`;
- ignore generated archives unless intentionally distributed;
- keep release artifacts in GitHub Releases, not in the main repo tree.

Suggested `.gitignore` additions:

```gitignore
GH_logs/
archives/*.zip
archives/*.tar.gz
```

Be careful if some `GH_logs` files are deliberately used for CI tests. If not, remove.

---

## 25. Public notation page

### 25.1 Current issue

The PLQR design note contains a useful notation mapping, but it is buried in markdown notes. A public user needs a concise notation page.

### 25.2 Recommendation

Add a vignette or pkgdown article:

```text
Notation and package conventions
```

Include:

| Package | Paper | Meaning |
|---|---|---|
| `X` | \(X_i\) | localization covariates |
| `x0`, `X0` | \(x\), target points | local target(s) |
| `s` | \(s_n\) | DNN order |
| `k` | \(k\) | fixed nearest-neighbor count |
| `theta` | \(\theta\) | PLQR/GMM target coefficient |
| `g` | \(g_0(U)\) | PLQR nonparametric shape |
| `mD` | \(m_0(U)\) | residualizing nuisance |
| `density` | local threshold density | PLQR Jacobian nuisance |
| `score` | \(\psi_i\) | observation contribution |
| `moment` | \(\sum_i w_i\psi_i\) | weighted local moment |

This will reduce the burden on every individual help page.

---

## 26. Immediate implementation plan

### Phase 1: PLQR breaking rename

1. Rename `beta` to `theta_eval`.
2. Rename `q` / `q_col` to `g` / `g_col`.
3. Rename demo columns `q_oracle` to `g_oracle`.
4. Rename `q_fn` to `g_fn`, `q_x0` to `g_x0`.
5. Rename `return = "g"` to `return = "score"` or `"psi"`.
6. Update all tests, docs, vignettes, README.

Expected result: PLQR package notation matches the paper.

### Phase 2: PLQR result-field cleanup

1. Split `se` into `theta_se` and `moment_se`.
2. Split `vcov` into `theta_vcov`.
3. Add `theta_omega` if useful.
4. Add `meta$theta_source = "supplied"`.
5. Replace `meta$se_route` with explicit coefficient/moment inference metadata.
6. Update `print`, `summary`, `confint`, `predict`, `plot_data_estimates`.

Expected result: no ambiguous PLQR standard-error fields.

### Phase 3: Public API pruning

1. Decide stable/developer/experimental/internal function groups.
2. Unexport functions not meant for public use.
3. Rename rank registry functions to `*_rank_weight`.
4. Decide whether raw `dist_*` functions stay exported.
5. Move GMM/IV to experimental reference group or hide from README.
6. Hide TDNN/jackknife infrastructure from front-page docs unless needed.

Expected result: small, intentional public namespace.

### Phase 4: Documentation rewrite

1. Rewrite README with paper notation.
2. Update conditional Z vignette to use renamed PLQR arguments.
3. Add notation/conventions vignette.
4. Add public status note: experimental package, breaking changes until first public release.
5. Update NEWS to describe the breaking cleanup.

Expected result: public docs tell one coherent story.

### Phase 5: Diagnostics and guardrails

1. Add direct-sandwich `normalize = FALSE` guard.
2. Add Jacobian conditioning diagnostics.
3. Add tests for singular/ill-conditioned Jacobian.
4. Add tests for renamed API and removed aliases.
5. Regenerate docs and run full tests/checks.

Expected result: safer PLQR inference and better public reliability.

---

## 27. Priority table

| Priority | Issue | Action |
|---:|---|---|
| 1 | `q` / `q_col` / `q_oracle` legacy PLQR naming | Rename to `g` / `g_col` / `g_oracle` |
| 1 | `beta` input versus `theta` output | Rename to `theta_eval` or `theta` |
| 1 | PLQR `se` overloaded | Split into `theta_se` and `moment_se` |
| 1 | README PLQR no-density example | Replace with density-gated coefficient inference example |
| 1 | `score_plqr_component()` hardcodes `Y` | Add `y_col` or unexport |
| 1 | ignored `k` in `score_plqr_component()` | Remove if public |
| 2 | rank registry names | Rename to `register_rank_weight()`, etc. |
| 2 | `return = "g"` conflict | Rename to `return = "score"` or `"psi"` |
| 2 | GMM/IV advertised too strongly | Mark experimental or remove from README |
| 2 | raw distance and low-level exports | Prune namespace |
| 2 | direct sandwich normalization | Forbid or normalize internally |
| 2 | missing Jacobian diagnostics | Add target-level diagnostics |
| 3 | malformed roxygen | Fix `@examplesIf TRUE` |
| 3 | generated artifacts/debug logs | Remove/ignore before public release |
| 3 | central notation page | Add vignette/article |

---

## 28. Suggested codex prompt for the first cleanup pass

```text
Perform a breaking pre-public PLQR naming cleanup in RdistNN. Rename the PLQR nonparametric shape nuisance from q/q_col/q_oracle/q_fn/q_x0 to g/g_col/g_oracle/g_fn/g_x0 throughout the user-facing API, demo data, score helpers, tests, README, vignettes, roxygen docs, and generated Rd files. Rename the PLQR coefficient input beta to theta_eval throughout plqr_dnn(), score helpers, tests, and documentation, while keeping result theta as the matrix of evaluated coefficients and adding meta$theta_source = "supplied". Do not preserve q_col/q/beta aliases; this is a pre-public breaking change. Regenerate documentation and update examples so they use the new names.

While making this change, also rename PLQR return extra "g" to "score" or "psi" consistently, so that g is reserved for the shape nuisance g_0(U). Update tests to check that old argument names fail clearly rather than silently aliasing. Keep the numerical behavior unchanged: the localized PLQR moment, direct sandwich J/S formulas, and prediction behavior should match the pre-cleanup outputs up to renamed fields. Run the full test suite and update NEWS with a short breaking-change note.
```

---

## 29. Final assessment

The package is in a good position for a cleanup pass precisely because the recent PLQR maturity work improved the tests and documentation. That means a breaking rename can be done safely: the tests should catch unintended numerical changes.

The main risk is not code failure. The main risk is public semantic drift. If released as-is, users will see:

- `beta` as input but `theta` as output;
- `q_tau(X)` used for a nuisance that is really \(g_0(U)\);
- `se` meaning either coefficient SE or moment SE depending on route;
- `score_plqr_component()` hardcoding `Y`;
- rank-weight registries called `rank`;
- GMM/IV modules presented alongside theorem-complete CATE functionality;
- README examples that do not reflect density-gated PLQR inference.

All of these are fixable. The best time to fix them is before public release. Once the package is public, every one of these becomes a deprecation problem rather than a simple cleanup.
