# Conditional LATE Implementation Plan

> **For agentic workers:** use superpowers:executing-plans for future native
> task-by-task implementation. Delegation is not required by this plan.
> Steps use unchecked boxes because no implementation has started.

**Goal:** deliver a documented, tested conditional-LATE workflow using paired
exact-DNN component averages and frozen-score covariance.

**Architecture:** a numeric ratio engine owns locality and paired deletions;
an IV signal constructor owns augmentation and an outer wrapper owns input
resolution. Classed results share reporting methods, while learner adapters
and optional confidence sets extend the same stable component contract.

**Tech Stack:** existing R >=4.1.0, base/stats, Rcpp/C++17 infrastructure,
testthat edition 3, knitr/rmarkdown and the declared roxygen2 toolchain.
No new runtime dependency is needed for the initial feature.

**Spec:** [01 scientific scope](01_estimator_scope_and_source_map.md),
[02 public contracts](02_public_interfaces_and_results.md),
[03 signals](03_augmented_signals_and_validation.md),
[04 computation](04_localization_and_deletion_engine.md),
[05 inference](05_inference_confidence_sets_and_diagnostics.md),
[06 studies](06_learners_bounded_model_and_studies.md),
[07 tests](07_unit_integration_and_regression_tests.md) and
[08 documentation](08_documentation_examples_and_vignettes.md).

## Global constraints

- Scope is a scalar conditional target, binary instrument/treatment, exact
  without-replacement DNN, method="exact",k=1L and common full/deleted s.
- User-facing estimators support matrices/vectors and data frames/tibbles;
  internals use canonical numeric arrays.
- Nuisance training stays external to the initial clate_dnn() wrapper.
- Freeze nuisance predictions, fitted target and geometry for covariance.
- Preserve signed first-stage signals and finite negative denominators.
- Keep paired estimator covariance separate from theta_omega=N*theta_vcov.
- Do not hand-edit NAMESPACE, man/*.Rd or RcppExports.*.
- Preserve existing CATE/NP, PLQR/IV-GMM and frozen replication behavior.
- Optional dependency paths are guarded; no network or sibling checkout
  is required by tests/examples.
- These are planning documents, not authorization to implement the feature
  during the present notes-only task.

## Review focus

The highest-risk inputs have explicit owners in Note 07: independently
misaligned prediction rows (signal contract and adapters); zero/negative
denominators (core/status tests); singular paired covariance with positive
ratio variance (perfect compliance); tied or data-dependent geometry
(engine parity); and m=1 versus m>1 covariance/target shapes (interfaces
and methods). Each must have a public or independent-reference test before
the corresponding task is accepted.

## 1. File ownership

| Proposed file | Single responsibility |
| --- | --- |
| `R/score-clate.R` | IV signal constructor and augmentation-specific validation |
| `R/user-dnn-affine-ratio.R` | Generic fixed-array wrapper, paired computation and statuses |
| `R/user-dnn-clate.R` | Dual-interface IV input resolution and delegation |
| `R/user-dnn-affine-ratio-methods.R` | Stable coefficient/covariance/interval/summary/predict methods |
| `R/user-dnn-clate-diagnostics.R` | Generic-ratio and IV-specific diagnostic adapters |
| `R/demo-data-clate.R` | Later bounded synthetic benchmark |
| `R/sim-clate-learners.R` | Later internal ridge prediction adapter |
| `R/clate-confidence-sets.R` | Optional quadratic set inversion |
| `tests/testthat/helper-clate-reference.R` and suites in Note 07 | Independent numerical and behavior acceptance |
| Two CLATE vignettes, roxygen and public doc owners in Note 08 | User workflow and scientific interpretation |

Existing shared code changes should be minimal. Reuse
user-dnn-core.R, utils-validate.R, exact weights and full-rank deletion
dispatch. If a shared signature must change, preserve old call semantics
and expand regression scope to affected estimators. Do not make a broad
core refactor a prerequisite for the new wrapper.

## 2. Task 1 — Signal constructor and independent algebra

**Files:** create R/score-clate.R,
tests/testthat/helper-clate-reference.R, test-clate-signals.R and the
signal-specific validation cases. Generate the constructor help/namespace
from roxygen.

**Consumes:** supplied numeric prediction/observation vectors.
**Produces:** iv_augmented_signals() and the exact iv_signals payload in
Notes 02–03.

- [ ] Write "augmented responses match explicit formulas" with explicit
  arrays, plus the one-sided negative-signal and perfect-compliance cases.
- [ ] Run the focused new test; confirm failure because the proposed
  constructor is absent, rather than from a broken fixture.
- [ ] Implement the exact Note 02 signature and Note 03 validation sequence,
  shared clipped propensity, row IDs and descriptive training metadata.
- [ ] Add the finite-distribution product-remainder test and all
  constructor misuse cases from Note 07.
- [ ] Regenerate help; run the focused signals/validation suite and commit
  the complete constructor/tests/docs slice.

**Accept when:** both signals match an independent formula, one-sided
equivalence holds and invalid arrays fail before localization.

## 3. Task 2 — Generic ratio point estimation and shapes

**Files:** create R/user-dnn-affine-ratio.R and core/validation suites.
Add focused shared-helper regression cases only where reuse needs them.

**Consumes:** canonical two-column signals and existing weight/tuning APIs.
**Produces:** affine_ratio_dnn() point estimates, components, neff,
target matrix, metadata, optional returns and per-target statuses.

- [ ] Write the four-row fixture's components/theta/shape assertions and
  batch zero/negative-denominator assertions with covariance disabled.
- [ ] Run the focused tests and confirm absence/contract failures.
- [ ] Implement both input routes, exact/k=1 restriction, default tuning
  resolution, target/ID validation and ratio status rules.
- [ ] Test n=1,s=1 point estimation, unsupported controls, optional-return
  warnings and no dropped batch targets.
- [ ] Generate help, verify focused tests and commit this point-estimation
  slice without advertising covariance or S3 methods not yet present.

**Accept when:** matrix/data-frame paths produce identical core fields and
all valid-target status cases are explicit.

## 4. Task 3 — Paired frozen covariance and engine parity

**Files:** extend the ratio engine and core suite; add
test-clate-engine-parity.R. Existing full-rank engine files change only if
a demonstrated shared bug or missing contract requires it.

**Consumes:** full components, common ranks, two frozen signals and resolved s.
**Produces:** component_vcov, theta_vcov, theta_omega, se, optional deletion/
pseudo-value arrays and covariance statuses.

- [ ] Extend the four-row test to require V_theta=34,
  N=18/7 and Omega=612/7, including its zero deleted denominator.
- [ ] Confirm the covariance assertions fail before implementation.
- [ ] Route fast computation through existing matrix-signal deletions;
  implement "loo" as a frozen-rank reference, with no geometry refit.
- [ ] Form paired covariance and the stable centered-score ratio variance;
  retain exact-zero and numerical-failure statuses.
- [ ] Run exhaustive small-subset checks, tied-rank cases, R/available-C++
  parity and existing prefix/scalar-reference suites.
- [ ] Generate help and commit the paired inference slice.

**Accept when:** all arrays and scales agree with independent references,
including m=1 and m>1. No deleted nuisance models or ratios are used.

## 5. Task 4 — Public CLATE wrapper and dual interfaces

**Files:** create R/user-dnn-clate.R,
test-user-clate-dnn.R, test-user-clate-dnn-data-interface.R and
test-user-clate-dnn-validation.R.

**Consumes:** signal constructor and generic ratio engine.
**Produces:** clate_dnn(), class c("clate_result","affine_ratio_result"),
component-role and observation/prediction metadata.

- [ ] Write wrapper-versus-explicit-composition tests, matrix/data-frame/
  optional-tibble equivalence and unrestricted/one-sided cases.
- [ ] Confirm failure at the missing wrapper and then implement the exact
  Note 02 signature without duplicating formulas or covariance.
- [ ] Test perfect-compliance equality to CATE, full instrument relabeling,
  wrong-arm nuisance mistakes and column/row-ID errors.
- [ ] Validate training declarations and optional nuisance returns; retain
  A and used/raw propensity needed by IV diagnostics.
- [ ] Generate help, run all CLATE/core suites and commit.

**Accept when:** every route computes the same supplied-prediction
functional and both sample interfaces are documented and tested.

## 6. Task 5 — Reporting, prediction and diagnostics

**Files:** create the methods/diagnostics files and their test suites;
register only implemented S3 methods through roxygen.

**Consumes:** the stable Note 02 result and retained frozen model state.
**Produces:** coef/vcov/confint/summary/print/predict plus diagnostics.

- [ ] Write method shape and missing-covariance tests, including a batch
  with undefined ratios and zero score variance.
- [ ] Implement marginal covariance access and pointwise Wald reporting;
  confidence bounds remain NA when wald_available is false.
- [ ] Implement prediction by relocalizing the same paired signals and s.
  Compare original/new targets with explicit engine calls and serialized
  state; confirm no nuisance fitting occurs.
- [ ] Implement diagnostic reconstruction of paired variation, N,
  concentration and available instrument information.
- [ ] Run the relevant methods/target-matrix/diagnostic regressions,
  generate help and commit the reporting slice.

**Accept when:** users can inspect components and every target's status
through ordinary methods, with correctly scaled uncertainty.

## 7. Task 6 — User documentation and initial feature acceptance

**Files:** create the two Rmd articles and update README.Rmd,
vignettes/getting-started.Rmd, notation-and-conventions.Rmd,
R/RdistNN-package.R, _pkgdown.yml, NEWS.md and tests/README.md.

**Consumes:** Tasks 1–5 and the teaching contract in Note 08.
**Produces:** a complete supplied-prediction user feature with experimental
lifecycle, runnable examples and discoverable help.

- [ ] Execute the standalone matrix/data-frame examples and add their
  equivalence assertions to the behavior tests where needed.
- [ ] Write and build both articles, keeping optional learner dependencies
  guarded and network/sibling access unnecessary.
- [ ] Generate help and README from their sources with the declared tools;
  inspect rendered tables, equations, status output and article links.
- [ ] Run focused tests, the full package suite, alignment checks and an
  installable-source check; build pkgdown locally.
- [ ] Review all warnings/notes and supported-platform CI evidence before
  declaring the implemented feature release-ready.
- [ ] Commit the documentation/release-note slice. Tagging, announcing or
  deploying is a separate author-requested action.

**Accept when:** Tasks 1–6 deliver all initial interfaces, validation,
numeric contracts, tests and user documentation. Do not require a new SIPP
analysis, confidence-set API or practical-forest theorem to accept this
well-defined first feature.

## 8. Task 7 — Bounded generator and designated learner

This is a later study/example stage, not a hidden missing piece of
supplied-prediction estimation.

**Files:** proposed demo/learner files and their dedicated tests.
**Consumes:** the structural model, exact objective and schedule in Note 06.
**Produces:** clate_demo_data() and internal
.clate_ridge_predictions(data,M,lambda,normalizer_n,delete).

- [ ] Write exact structural-moment and oracle-constant expectations.
- [ ] Implement generator/return schema and verify both core and
  heterogeneous target definitions.
- [ ] Write explicit small retained-refit tests before adding the ridge
  adapter; verify the same ambient normalizer in every fit.
- [ ] Implement the adapter, projection records and seed/row metadata.
- [ ] Run deterministic checks and small execution studies, document scope
  and commit. Production simulation remains a separate study decision.

**Accept when:** generated truths and fitted algorithms match the stated
construction; no claim of practical forest validity or coverage follows
solely from these tests.

## 9. Task 8 — Optional ratio confidence sets

**Files:** R/clate-confidence-sets.R, its suite, help and inference-article
extension.
**Consumes:** existing components and component_vcov.
**Produces:** ratio_confidence_set() and private .ratio_quadratic_set().

- [ ] Write every algebraic geometry and coherent PSD-input case from
  Note 07, including exact-zero-denominator degeneracies.
- [ ] Implement the Note 05 numeric contract, stable roots and explicit
  near-boundary/unresolved statuses.
- [ ] Test membership by substitution, rescaling, serialization and
  infinite endpoints.
- [ ] Generate help and article examples; state the fixed-candidate
  justification and uniform weak-IV limitation precisely.
- [ ] Run focused and affected full checks and commit this separate feature.

**Accept when:** set geometry and numeric statuses are correct without
turning all outputs into finite Wald endpoints.

## 10. Task 9 — Optional tuning, plots and learner adapters

Only begin after a concrete user need and a contract update:

- Reuse the integrated candidate engine and candidate_fit() conventions for
  grid sensitivity. Preserve one common selected order across deletions.
- Reuse plot-data extractors; add component/ratio/status views and table
  equivalents. Add numerical plot-data tests before optional visual snapshots.
- Add practical prediction adapters with explicit all-tree/OOB/cross-fit
  semantics. Guard Suggested dependencies and test ID alignment.
- Profile before target chunking or C++ changes. Measure allocation and
  runtime outside brittle unit-test thresholds.
- Coordinate Stata parity through shared small fixtures if companion work
  is requested; do not mutate that repository as a side effect here.

These tasks cannot be treated as already implemented by a collection of
names or a successful array smoke test. Scope, files and acceptance tests
must be added to their own execution slice before coding.

## 11. Verification commands for the implemented feature

From the package root, select checks by changed surface and use installed
development dependencies:

~~~bash
Rscript -e 'devtools::test(filter = "^(clate-|user-clate|user-affine-ratio)")'
Rscript -e 'devtools::test()'
Rscript tools/ci/check_tests_declared_deps.R
Rscript tools/ci/check_pkgdown_reference_coverage.R
Rscript tools/ci/check_generated_file_drift.R
Rscript tools/ci/check_repository_hygiene.R
~~~

For a feature release, follow [CONTRIBUTING](../../CONTRIBUTING.md) and
[RELEASING](../../RELEASING.md), refreshing historical version-specific
instructions rather than copying the old 0.1.0 tag:

~~~r
devtools::document()
devtools::build_readme()
dir.create("archives", showWarnings = FALSE)
src <- devtools::build(path = "archives", vignettes = TRUE, manual = FALSE)
rcmdcheck::rcmdcheck(src, args = "--no-manual", error_on = "warning",
                   check_dir = "archives/check")
pkgdown::build_site()
~~~

Require no check errors/warnings; inspect any notes against the baseline.
Confirm compiled/fallback parity on Linux x86-64, Windows x86-64 and both
macOS architectures, with R release/devel/oldrel as configured. Never
substitute a local check for terminal hosted-platform evidence when
claiming that matrix passed.

For **this planning collection**, verify all local links and source-owner
claims, parse proposed R examples, check the independent exact fixture,
review the nine-note requirement and run whitespace/hygiene checks.
Package tests or builds cannot validate unimplemented APIs and are not
required merely to commit repository-only planning prose.

## 12. Handoff and completion record

Before future implementation, recheck live package/paper snapshots and
resolve any changed shared signatures. Update a task's status with the
actual verified commit and commands when it lands; do not check a box
because code was merely written. Keep a requirements-to-evidence table
covering formulas, interfaces, failure semantics, methods, examples,
vignettes and generated files.

Use small coherent commits on the repository's default branch unless the
author requests another workflow. Preserve user-owned changes and stage
only the relevant files. Reuse successful verification while content is
unchanged; fetch/rebase only when upstream movement is plausible or a push
rejects. A successful push confirms publication of that commit. Publishing
these notes does not publish the proposed estimator.
