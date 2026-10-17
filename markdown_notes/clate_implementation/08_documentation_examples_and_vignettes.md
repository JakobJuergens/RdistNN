# 08 — Help, executable examples, vignettes and pkgdown

[Index](README.md) · [Previous](07_unit_integration_and_regression_tests.md) ·
[Next: implementation sequence](09_implementation_sequence_and_release.md)

**Status:** documentation plan for a future feature. No Rmd article, help
topic or exported example is added in the present planning pass. All code
below uses proposed APIs and becomes executable only with their
implementation.

## 1. Reader tasks and help-topic ownership

The documentation should let a user answer four questions: which effect is
being reported, what each prediction array means, how uncertainty is
computed, and how to inspect an unstable first stage without changing the
estimator. Begin with the economic roles and supplied-prediction workflow;
put proof bookkeeping in the conceptual article.

| Future topic | Roxygen source | Required content |
| --- | --- | --- |
| `iv_augmented_signals` | `R/score-clate.R` | Both signal formulas, arm indexing, clipping, one-sided restriction, IDs and returned arrays |
| `affine_ratio_dnn` | `R/user-dnn-affine-ratio.R` | Dual interfaces, common weights, ratio of averages, covariance fields/scales and exact-zero behavior |
| `clate_dnn` | `R/user-dnn-clate.R` | Identification assumptions, reporting versus adjustment, five predictions, training declaration and experimental lifecycle |
| `affine_ratio_result` methods | `R/user-dnn-affine-ratio-methods.R` | coef/vcov/confint/summary/predict, target indices, missing covariance and per-target statuses |
| CLATE diagnostics | `R/user-dnn-clate-diagnostics.R` | Component SEs, N, denominator diagnostics, clipping and concentration |
| Later `clate_demo_data` | `R/demo-data-clate.R` | Known synthetic target, explicit oracle columns, latent audit payload and seed |
| Optional `ratio_confidence_set` | `R/clate-confidence-sets.R` | Set-valued output, rays/empty/all-real geometries and pointwise coverage qualifications |

Document low-level numeric inputs separately from user-facing data-frame
extraction. In each high-level help page show both matrix/vector and
data-frame usage and explain that tibbles use the same tabular path.
Do not imply a formula/NSE interface that is not implemented.

## 2. Essential terminology and variance table

Use A for instrument and D for treatment throughout this collection.
In prose use instrument assignment and participation/treatment explicitly.
The adjustment set V may be larger than reporting X, and the output
averages over the remaining adjustment distribution at each reporting
profile. The five functions are instrument propensity, two outcome-arm means
and two treatment-arm means.

Include this scale table in help and the inference article:

| Object | Interpretation |
| --- | --- |
| `components[, "numerator"]` | Reduced-form estimate for CLATE |
| `components[, "denominator"]` | First-stage estimate for CLATE |
| `component_vcov[[j]]` | Estimator-scale paired component covariance |
| `theta_vcov[[j]]` | Estimator-scale ratio variance |
| `theta_omega[[j]]` | N_j times ratio variance, local-limit covariance estimate |
| `se[j]` | sqrt(theta_vcov[[j]][1,1]) |
| `neff[j]` | Rank-weight concentration, not the number of compliers |

Specify that frozen-score covariance holds nuisances, fitted target,
geometry and s fixed. It recomputes retained weights. The theoretical
deleted nuisance fits do not appear in routine standard-error evaluation.
No second ESS division is performed after theta_vcov is returned.

## 3. A small standalone example

Use an explicit perfect-compliance model for a fast first example. It
demonstrates the signatures and reduction to CATE without a suggested
learner dependency:

~~~r
set.seed(20261017)
n <- 32L
X <- matrix(seq(-1, 1, length.out = n), ncol = 1L,
            dimnames = list(NULL, "reporting_x"))
A <- stats::rbinom(n, 1L, .5)
D <- A
baseline <- .2 * X[, 1L]
Y <- baseline + .5 * D + stats::rnorm(n, sd = .1)
pred <- list(e = rep(.5, n), m_y0 = baseline,
             m_y1 = baseline + .5,
             m_d0 = rep(0, n), m_d1 = rep(1, n))
sig <- iv_augmented_signals(
  Y, D, A, pred$e, pred$m_y0, pred$m_y1,
  pred$m_d0, pred$m_d1, training = "oracle"
)
fit <- affine_ratio_dnn(
  X = X, x0 = 0, numerator = sig$phi_y,
  denominator = sig$phi_d, s = 6L
)
fit$components
coef(fit)
vcov(fit)
confint(fit)
~~~

The regressions in this synthetic example are known by construction.
Replace them with aligned supplied predictions in an application;
changing the training label does not change their mathematical status.
The target is .5, but a single simulated interval is not evidence of
coverage. Use tolerance assertions in tests rather than hard-coded
realized estimates in narrative text.

Follow with a data-frame `clate_dnn()` example using explicit columns:

~~~r
d <- data.frame(reporting_x = X[, 1L], Y, D, A,
                e = pred$e, m_y0 = pred$m_y0, m_y1 = pred$m_y1,
                m_d0 = pred$m_d0, m_d1 = pred$m_d1)
fit_iv <- clate_dnn(
  data = d, x_cols = "reporting_x",
  y_col = "Y", d_col = "D", a_col = "A", e_col = "e",
  m_y0_col = "m_y0", m_y1_col = "m_y1",
  m_d0_col = "m_d0", m_d1_col = "m_d1",
  x0 = 0, s = 6L, training = "oracle"
)
summary(fit_iv)
diagnostics(fit_iv)
~~~

The example's object names are local and self-contained. Roxygen may use
an even shorter version; a topic must not depend on an object created by
another topic's example.

## 4. Vignette one: conditional-late.Rmd

Proposed title: **Conditional LATE from supplied nuisance predictions**.
Audience: an applied researcher who already has predictions and wants a
correct local estimate with transparent components.

The article should proceed through:

1. A binary encouragement/participation question, the complier estimand and
   the difference between adjustment V and reporting X.
2. A small fully specified synthetic instrument/treatment example. An
   inline base-R generator can support the initial feature; use the later
   exported demo only once it exists and has tests.
3. Construction of the five oracle/supplied prediction arrays, with their
   roles displayed in a table. Explain which are estimated in real data.
4. Matrix and data-frame calls, including one tibble equivalence example
   guarded by requireNamespace("tibble",quietly=TRUE).
5. A two- or three-target component table with ratio, SE, interval, N and
   status. Do not print only the final ratio.
6. A selected one-sided restriction with explicit observed coding checks
   and explanation that the structural assumption comes from the design.
7. New reporting targets using predict(), with no nuisance refitting.
8. A deliberately zero-denominator fixed-array example showing retained
   components and unavailable ratio/Wald interval.

Keep evaluated default chunks to small base-R calculations. Optional
estimated-learner chunks must be clearly marked, use Suggested packages
only when available, and avoid costly forest fitting during package
checks. A vignette should remain useful without optional learners.
Do not access live URLs or the sibling manuscript at build time.

## 5. Vignette two: conditional-late-inference.Rmd

Proposed title: **Paired frozen-score inference for conditional ratios**.
Audience: users/developers who need to understand covariance and diagnostics.

Use the four-row fixture from Note 04 to calculate weights, retained
averages, pseudo-values, paired covariance and ratio variance line by line.
Show that a deleted denominator can be zero while the frozen-score
calculation remains finite. Include a counterexample where ignoring the
off-diagonal covariance changes the standard error.

Then show perfect compliance's singular component covariance, signed
augmented first-stage signals, numerical versus statistical denominator
diagnostics, and the two-stratum estimand comparison (3/2,1,18/13).
Explain the existing local IV-GMM article as a separate score/workflow.
Cross-link rather than replace it.

Discuss orthogonal mean bias and deletion stability as separate
requirements. A package result is a defined computation; the full-sample
inference theorem still has learner/design assumptions. Use the
bounded-series example as one future sufficient construction, while raw
financial outcomes and practical forests have their own scope.

An optional final section can demonstrate ratio_confidence_set() after
that stage lands. Preserve unbounded/ray output in tables. Explain the
fixed-candidate pointwise justification and the separate uniform weak-IV
question. Before that API exists, describe the extension in prose without
unevaluated code presented as an installed feature.

## 6. Public integration and lifecycle

Edit `README.Rmd`, not generated README.md, to add a short experimental
CLATE pointer after the established CATE/regression introduction. Keep
those initial-release workflows prominent. Update
`vignettes/getting-started.Rmd` with a concise route to the CLATE article,
and `vignettes/notation-and-conventions.Rmd` with A,D,phi_Y,phi_D,
theta_vcov/theta_omega and the status contract.

Update `R/RdistNN-package.R` and NEWS.md when user-visible functions land.
Add a reference group **Experimental — conditional LATE** and the two
articles to `_pkgdown.yml`. Each exported function/method must have
generated help and appear in reference coverage. New constructor/engine
exports should be presented by role rather than buried in an unrelated
GMM group. Do not emit runtime warnings solely for experimental lifecycle.

Proposed method topics include print/summary/coef/vcov/confint/predict and
diagnostics. Group inherited base methods clearly; do not register
nonexistent CLATE methods to obtain a long reference listing.
Do not change DESCRIPTION's supported-workflow claims to imply theorem
coverage for every practical learner.

## 7. Generation, checks and accessibility

Use the declared roxygen2 version (8.0.0 at the inspected baseline).
Regenerate NAMESPACE and man/*.Rd with devtools::document(); regenerate
README.md with devtools::build_readme(). Do not hand-edit generated output.
Update the declared generator version only through a separate deliberate
toolchain change if it becomes necessary.

Execute examples and build both articles with the package installed.
Inspect single-target and multi-target printed tables, missing intervals,
long column names and units. Figures need useful captions and alt text,
and an adjacent numeric table for statuses and components. Use color plus
line/shape differences; do not encode weak/undefined output in color alone.

Run existing reference-coverage and generated-drift checks. Validate local
article links in rendered HTML as well as Rmd source. Check optional
dependency absence and a clean build without sibling repositories or
network access. Re-run examples after any signature or field change.

## 8. Documentation acceptance evidence

Completion requires real generated topics, executable matrix/data-frame
examples, built articles, correct reference navigation and a reviewed
rendered result. A grep for the word CLATE does not establish that users
can execute the workflow. Conversely, prose-only repository planning notes
do not require a package/site rebuild today: their links, source claims,
numeric examples and syntax are the changed surface.

Publishing package sources is distinct from deploying pkgdown. The
inspected pkgdown workflow runs on a published release or manual dispatch.
Only request that deployment when the feature is implemented and the
author asks for the public documentation refresh; this planning push does
not require or authorize a release tag.
