# Test-suite map

Tiny-Beta and close-order fixtures cover bounded cases. They do not certify
all subnormal tail support, order-statistic mixtures or signed TDNN blending;
the unresolved limitation is recorded in the
[numerical reference](../markdown_notes/design/19_weighted_accumulation.md).
Fixed-table engine agreement does not establish correct table construction.

Run tests from the repository root with `devtools::test()`. The R package
suite lives directly in `tests/testthat/`; `tests/ci/` contains separate
Python checks for repository and workflow policy.

Replication subprocesses and multisession workers need a matching installed
RdistNN on their library path. For a full candidate run, install the checked
archive into a private library and expose it through `R_LIBS` before running
the source suite; source loading alone is confined to the current R process.
The SIPP smoke integration runs during checks of the installed archive. Its
test skips under pkgload source loading, even when another installation is
visible on the library path.

## Finding a test

Filenames describe the behavior under test. Use the matching `R/` module as
the first place to look, then the relevant integration or numerical suite.

| Area | Files in `tests/testthat/` |
| --- | --- |
| Structured errors, bounded warning copy, base warning controls and selective recovery | `test-error-reporting.R`, `test-reporting-recovery.R` |
| Prepared group facts, deletion-refit cache retention, named distances and serialized reuse | `test-distance-tie-localization.R` |
| Saved geometry and shared inputs | `test-fitted-geometry-contract.R`, `test-input-route-contract.R`, `test-prediction-selectors.R` |
| Weight control and scalar/grid consumers | `test-weight-control-contract.R`, `test-tuning-consumer-contract.R` |
| Shared grid callbacks, failure rows, duplicate keys and cached target/inference reuse | `test-tuning-grid-engine.R` |
| Local origins, offsets, translated SD and covariance | `test-centering-translation.R` |
| Neighbor extras and CATE nuisance contracts | `test-result-neighbors.R`, `test-cate-nuisance-contract.R`, `test-score-cate-accumulation.R` |
| Marginal prediction allocation limits, explicit-deletion accuracy and filtered-profiler positive/negative controls | `test-np-prediction-memory.R` |
| CATE signal provenance, original creation version and current serialized reuse | `test-cate-provenance.R` |
| Demo subset matrix dimensions, batch scores and preserved RNG draws | `test-demo-subsets.R` |
| Test fixture state, public seed scopes and bounded oracle simulation backend/dependency contract | `test-fixture-rng-state.R`, `test-seeded-rng-scope.R`, `test-simulation-rng-backends.R`, `test-demo-target-rng.R` |
| Exact-table represented parity and near-full support | `test-exact-table-portability.R` |
| Weighted CDF boundary indices, physical radius and outcome SD ranges | `test-weighted-distribution-range.R`, `fixtures/weighted-cdf.csv` |
| Component contracts | `test-api-*`, `test-distance-*`, `test-order-*`, `test-score-*`, `test-utils-validate.R` |
| Captured factory configuration, group callbacks and package-owned S3 identities | `test-distance-factory-state.R`, `test-score-factory-state.R`, `test-gmm-factory-state.R`, `test-class-names.R` |
| Experimental GMM callback recovery and current result arrays | `test-gmm-recovery-state.R` |
| Plot warning controls and reuse of retained caption facts | `test-plot-recovery-contract.R` |
| Singleton rank-weight batches and scalar score callback validation/data binding | `test-api-ranks-singleton-batch.R`, `test-api-scores-scalar-batch-validation.R` |
| Scalar-only distance batches, singleton weight/estimator consumers and shared setup/context | `test-distance-scalar-batch-consumers.R` |
| Compensated weighted sums, native/R parity, retained deletion corrections and buffer boundaries | `test-compensated-accumulation.R` |
| Native rank/score index bounds, missing indices and serial/parallel error reporting | `test-native-index-boundaries.R` |
| Numerical references and optimized paths | `test-rank-weight-*`, `test-jackknife-rank-*`, `test-dnn-scalar-reference.R`, `test-dnn-jackknife-prefix.R`, `test-tdnn-*`, `test-compiled-fallback-parity.R`, `test-errors-fastpaths.R` |
| Exact comparison semantics, bounded failure messages and closure-safe diagnostics | `helper-comparisons.R`, `test-comparison-assertions.R` |
| Native/R distance arithmetic, large-offset and accumulation-sensitive ranks, scalar/batch agreement and FMA-sensitive ties | `test-distance-scaled-arithmetic.R` |
| Scaled L1/L2 products below the coordinate range, final-rounded subnormal norms, shared range boundaries, ordinary arithmetic and resulting distance-block weights | `test-distance-scaled-underflow.R` |
| Huge/tiny finite L2 norms, shared exceptional boundaries, exact ranks and serialized prediction | `test-distance-range-stability.R` |
| Close-order TDNN coefficients and tiny Beta interval masses against independent references | `test-weight-numerical-stability.R` |
| Support-aware jackknife centering across estimators, compiled/R and fast/naive paths | `test-jackknife-support-centering.R` |
| Exact-distance block averaging, survivor-weight redistribution, actual R/C++ deletions, paired covariance, permutation invariance, reuse and target-specific diagnostics/plots | `test-duplicate-block-averaging.R`, `helper-duplicate-blocks.R`, `test-candidate-weight-targets.R` |
| Target-specific exact-distance groups, symmetric and support-edge ties, signed mass, retained policy/groups, actual deletions and legacy-fit rejection | `test-distance-tie-localization.R`, `test-canonical-covariate-ties.R` |
| Original-signal block prefix/suffix states, finite survivors after cancellation, signed zero-mass covariance, genuine subnormals and interleaved rank masses | `test-block-numerical-regression.R` |
| Weighted origins with subnormal coefficients or nonzero mass, raw normalization and ordinary tiny-tail stability | `test-weighted-underflow-origin.R` |
| Shared prediction targets, stored selectors and candidate views | `test-prediction-targets-common.R`, `test-prediction-selectors.R`, `test-candidate-fit.R`, `test-candidate-retained-state-validation.R` |
| Complete current reporting/model state, legitimate NULL fields and numerical-range failure at fitting/prediction | `test-retained-report-state.R`, `test-retained-model-fields.R`, `test-fitted-estimate-range.R` |
| CATE public, nuisance, provenance and supported numerical routes | `test-user-cate-*`, `test-cate-*`, `test-score-cate-*` |
| NP/TDNN identities, data interfaces, public validation, methods and diagnostics | `test-user-np-*`, `test-user-tdnn-np*`, `test-demo-data-npr.R` |
| Conditional LATE and paired ratios | `test-clate-*`, `test-user-clate-*`, `test-user-affine-ratio-*`, `test-weights-single-row.R` |
| Tuning and stored candidate views | `test-api-tuning.R`, `test-user-dnn-tuning.R`, `test-tuning-*`, `test-candidate-fit.R`, `test-user-*-tuning-grid.R`, `test-candidate-retained-state-validation.R` |
| Shared diagnostics validation and printing | `test-user-dnn-diagnostics.R` |
| Shared pointwise confidence levels across result, prediction, candidate and plot consumers, plus simulation interval metrics | `test-user-confidence-levels.R` |
| Estimator diagnostics: NP outcomes, CATE overlap, PLQR Gram/Jacobian | `test-user-np-diagnostics.R`, `test-user-cate-diagnostics.R`, `test-user-plqr-diagnostics.R` |
| Candidate diagnostics, metadata and plot integration | `test-tuning-grid-diagnostics.R`, `test-tuning-grid-diagnostics-plots.R`, `test-plot-*` |
| CATE simulations | `test-sim-cate-data.R`, `test-sim-cate-one-rep.R`, `test-sim-cate-run.R`, `test-sim-cate-summary.R`, `test-sim-cate-tuning*.R` |
| Experimental PLQR and GMM/IV | `test-user-plqr-*`, `test-user-dnn-gmm*`, `test-user-gmm-methods.R`, `test-gmm-*.R`, `test-score-plqr.R`, `test-demo-data-plqr.R`, `test-demo-data-iv.R`, `test-sim-iv-*` |
| Integration, replication, and documentation | `test-integration-*`, `test-paper-*`, `test-result-*`, `test-*-documentation-*`, `test-notation-*`, `test-readme-*`, `test-citation-*`, `test-cate-paper-diagnostics.R`, `test-s-default-paper-rule.R` |

The larger wrapper families separate estimator identities, data interfaces,
validation, and result methods where those contracts warrant separate files.
Generic full rank jackknife contracts belong to `test-jackknife-rank-full.R`;
TDNN weight contracts stay in `test-tdnn-weights.R`. PLQR direct-sandwich
arithmetic, result methods, prediction and diagnostics have separate
`test-user-plqr-{inference,methods,prediction,diagnostics}.R` owners. Density
validation and matrix/data-frame equivalence live in the PLQR validation and
data-interface suites. The generic score derivative API is in
`test-api-scores-diff.R`; the public PLQR score matrix is in
`test-user-plqr-score-matrix.R`. The integration pipeline suite retains one
complete geometry/weight/score scenario and checks independent extension
contracts in separate test blocks.

## Focused runs

`filter` is a regular expression matched against filenames after removing
`test-` and `.R`. Anchor both ends to select one file:

```r
devtools::test(filter = "^sim-cate-summary$")
devtools::test(filter = "^user-dnn-gmm-validation$")
devtools::test(filter = "^jackknife-support-centering$")
devtools::test(filter = "^distance-(scaled-arithmetic|scaled-underflow|range-stability)$")
devtools::test(filter = "^(retained-(report-state|model-fields)|fitted-estimate-range)$")
```

Select the relevant family when a change crosses contracts. NP, CATE and
PLQR families include the shared demo-subset and target-validation RNG suites;
simulation families also include the shared backend tests:

```r
devtools::test(filter = "^(user-cate|cate-|score-cate|prediction-(targets-common|selectors)$|candidate-(fit|retained-state-validation|weight-targets)$|demo-(subsets|target-rng)$|seeded-rng-scope$|fixture-rng-state$|user-confidence-levels$)")
devtools::test(filter = "^(user-(np-|tdnn-np)|demo-data-npr$|prediction-(targets-common|selectors)$|candidate-(fit|retained-state-validation|weight-targets)$|demo-(subsets|target-rng)$|seeded-rng-scope$|fixture-rng-state$|user-confidence-levels$)")
devtools::test(filter = "^(prediction-(targets-common|selectors)|candidate-(fit|retained-state-validation|weight-targets)|user-confidence-levels)$")
devtools::test(filter = "^(api-tuning|tuning-|candidate-|user-dnn-tuning|user-.*-tuning-grid|sim-cate-tuning|demo-(subsets|target-rng)$|seeded-rng-scope$|simulation-rng-backends$|fixture-rng-state$|user-confidence-levels$)")
devtools::test(filter = "^(plot-|user-(dnn|np|cate|plqr)-diagnostics$|tuning-grid-diagnostics|user-confidence-levels$)")
devtools::test(filter = "^(sim-cate|demo-(subsets|target-rng)$|seeded-rng-scope$|simulation-rng-backends$|fixture-rng-state$|user-confidence-levels$)")
devtools::test(filter = "^(clate-|user-clate|user-affine-ratio|weights-single-row|seeded-rng-scope$|fixture-rng-state$|user-confidence-levels$)")
devtools::test(filter = "^(user-plqr|user-dnn-gmm|user-gmm|gmm-|score-plqr|demo-data-(plqr|iv)|sim-iv|demo-(subsets|target-rng)$|seeded-rng-scope$|simulation-rng-backends$|fixture-rng-state$|user-confidence-levels$)")
```

For changes to fitting, retained state, prediction, reporting or signal ranges,
compose the family run with the focused retained-state/range group above:
`retained-(report-state|model-fields)|fitted-estimate-range`. These shared
contracts are not selected by the family filters alone.

Use the full suite for changes shared across families. The release check
also runs these tests against the installed package; see
[RELEASING.md](../RELEASING.md).

## References, engines and optional checks

Domain-specific `helper-*.R` files define reusable fixtures without fitting
models or seeding eagerly. `helper-test-rng.R` owns the independent test-only
RNG scope; fixture seeds are explicit and restore the caller state. The
simulation backend suite locally extends that scope to options and future
plans. Moved random test blocks scope their own draws so family selection
does not depend on earlier files. Single-suite fixtures stay beside their tests.
`helper-diagnostics-reference.R` supplies scalar empirical diagnostics;
`helper-dnn-reference.R` supplies explicit scalar deletions,
`helper-sim-cate.R` independent AIPW calculations, `helper-user-dnn-gmm.R`
custom models, and `helper-clate-reference.R` complete-subset enumeration.
CLATE tests retain complete-subset, explicit-deletion covariance and
confidence-set geometry references. Small simulations check execution, not
coverage or manuscript claims.

Force and establish each engine; test automatic fallback separately. Use
`options(RdistNN.use_compiled = FALSE)` for the R path and require actual native
symbols for native acceptance. Exact ranks, support, statuses and missingness
accompany quantity-specific numeric tolerances. The
[testing policy](../markdown_notes/design/09_testing_and_regression_policy.md)
and [compiled-code policy](../markdown_notes/design/10_rcpp_and_compiled_code_policy.md)
own these acceptance rules. Shared primitive changes also require their
estimator, candidate, covariance and saved-prediction consumers.

Use `.expect_identical_fast()` for whole fitted-state snapshots, byte vectors,
simulation tables and extreme numeric identity fixtures where recursive
failure diffs are costly. It preserves base `identical()` checks and gives a
bounded message without traversing retained callback environments. Detailed
testthat differences remain useful for small ordinary values. Reporting
immutability compares object/data state directly; serialization round trips
and frozen-file byte identity keep their dedicated coverage.

Tests must run independently through public interfaces. Exercise S3 methods
through their generics on fitted objects; missing unexported aliases must not
skip installed-package coverage. Helper loading follows testthat's
[special-files guide](https://testthat.r-lib.org/articles/special-files.html).
Optional packages and compiled capabilities have explicit guards; native
requirements in native checks are failures rather than silent skips.
The Linux R-release job runs installed tests in an isolated library containing
Rcpp, testthat and their strong dependencies. Reproduce it with
`Rscript tools/ci/check_minimal_dependencies.R path/to/RdistNN_VERSION.tar.gz`.

Visual snapshots require a deliberate opt-in:

```sh
RDISTNN_VDIFFR=true Rscript -e 'devtools::test(filter = "^plot-vdiffr$")'
```

Inspect changed SVGs in `tests/testthat/_snaps/plot-vdiffr/`. Ordinary runs
retain numerical plot-data checks without platform-dependent visual snapshots.

## Coverage

Use `Rscript tools/ci/run_coverage.R` for full installed-package coverage. It
copies current package inputs into a private workspace before covr instruments
them and requires both R and handwritten C++ execution. Reports and provenance
remain under ignored `archives/coverage/`; the existing manual coverage workflow
uses the same runner.

For a focused installed suite, pass a testthat filename filter, for example
`--filter '^(native-threads|compiled-fallback-parity)$' --require-native`.
Pure R selections can omit `--require-native`. An empty selection, failed test
or missing required native hit fails the run. Do not substitute `test_local()`
or `devtools::test()` inside covr: those helpers can reload uninstrumented code.
See [developer tools](../CONTRIBUTING.md#developer-analysis-tools) for dependencies,
report destinations and optional HTML output. Coverage identifies gaps;
independent references and actual engine assertions establish correctness.

## Package-maintenance checks

[CONTRIBUTING.md](../CONTRIBUTING.md#editing-and-checking) lists the commands.
`tests/ci/` uses isolated positive/negative fixtures for repository/workflow
policy, generated-file drift, candidate hygiene and function ownership.
Release checks run the installed-package suite; see
[RELEASING.md](../RELEASING.md).
