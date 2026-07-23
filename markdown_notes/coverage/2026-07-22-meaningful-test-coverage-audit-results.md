# Meaningful Test-Coverage Audit Results

Date: 2026-07-22

## Summary

The audit raised overall measured coverage from 86.38% to 88.50%, but the
percentage change is not the main result. The added tests close specific gaps
in installed-package execution, scientific score definitions, public demo-data
contracts, user-facing result methods, and compiled/fallback equivalence.

Coverage was measured with `covr::package_coverage()` in archived copies of
the relevant Git snapshots. Instrumented coverage builds were not run in the
live checkout.

## Material improvements

| Contract | Before | After | What is now protected |
|---|---:|---:|---|
| `np_demo_data()` | 50.00% | 100.00% | schema, reproducibility, all designs, oracle equality, injection, validation |
| `np_demo_subsets()` | 0.00% | 100.00% | dimensions, reproducibility, valid distinct indices, invalid order |
| `cate_demo_subsets()` | 0.00% | 100.00% | shared subset-sampling contract |
| `plqr_demo_subsets()` | 0.00% | 100.00% | shared subset-sampling contract |
| `plqr_score_matrix()` | 14.29% | 95.24% | defining formula, selector equivalence, schema and value validation |
| `summary.cate_result()` | 0.00% | 100.00% | installed-package S3 dispatch and summary structure |
| `confint.cate_result()` | 0.00% | 100.00% | installed-package S3 dispatch and selectors |
| `predict.cate_result()` | 44.00% | 96.00% | installed-package prediction paths and validation |
| `print.cate_result()` | 0.00% | 100.00% | visible identity, SE state, invisible return |
| `print.gmm_result()` | 0.00% | 100.00% | visible identity, SE state, invisible return |
| `summary.gmm_result()` | 33.33% | 100.00% | coefficient and overidentification schemas |
| `print.summary.gmm_result()` | 0.00% | 100.00% | conditional overidentification display |
| `score_plqr_component()` | 66.20% | 71.83% | compiled and pure-R construction/aggregation equivalence |

The CATE method gains came primarily from repairing a test-suite defect. The
methods were correctly registered, but the file-level guard searched for
unexported method symbols and therefore skipped the entire file against an
installed package. The repaired suite uses the public generics and fails if S3
registration is missing.

The CATE and PLQR custom-distance tests had a similar installed-layout issue:
their public contract was unnecessarily coupled to the internal
`ord_to_ranks()` helper. A test-local order-to-rank conversion now lets both
extension-point tests run against the installed namespace.

## Why some file percentages remain low

`R/demo-data-npr.R` is reported at 34.91% even though its two exported
functions, `np_demo_data()` and `np_demo_subsets()`, are now fully covered. The
same file also contains the larger internal `np_y_demo_sets()` fixture factory,
which is not part of the public API and is not used by installed-package tests.
The file-level number therefore understates coverage of its active public
surface.

This is one reason package and file percentages should be treated as routing
signals rather than objectives. Function-level ownership and behavioral risk
are more informative.

## Remaining gaps and disposition

### Intentionally not targeted

- `R/sim-cate-benchmarks.R` is dominated by optional comparison-method
  adapters. Exercising those paths requires optional third-party integrations
  and belongs in explicit benchmark/integration tests, not ordinary coverage
  padding.
- `R/utils-bench.R` contains benchmark harness utilities rather than estimator
  behavior. Timing or helper smoke tests would add hits without materially
  improving inferential confidence.
- `R/zzz.R::.onLoad()` only registers built-ins at package load. Public
  registry behavior is already exercised through the API tests.

### Inactive or internal-only

- `.validate_same_n()` has no call sites in the package or tests. A direct unit
  test would protect unused code rather than active behavior. A separate
  cleanup could remove it after confirming no downstream use is intended.
- `np_y_demo_sets()` is an internal fixture factory. Its lack of installed
  coverage does not leave `np_demo_data()` or `score_np_y()` untested.

### Residual defensive branches

The remaining uncovered expressions in distance, rank-weight, score, plotting,
and grid-validation helpers are mostly alternate error branches or short
dispatch fallbacks. They should receive tests when a bug, interface change, or
new supported input makes the branch behavior consequential. Adding one test
per uncovered line now would primarily optimize the percentage.

Small unexecuted S3 fallback wrappers such as `plot.cate_result()`,
`plot.gmm_result()`, and default candidate-fit error methods are legitimate
future contract candidates, but their underlying plot builders and primary
methods already have broad coverage. They do not outweigh the completed gaps
in this audit.

## Verification evidence

- Controlled mutation of the zero-noise NPR oracle expectation produced four
  failures, one for each supported design, before restoration.
- Controlled mutation of the PLQR threshold indicator made the direct formula
  test fail before restoration.
- Focused tests passed after every change.
- The full source-tree suite passed through `devtools::test()` with the one
  opt-in visual-baseline skip and two existing near-uniform-weight warnings.
- The full suite also passed through `testthat::test_dir(...,
  load_package = "installed")`; the former namespace-related skips disappeared.
- The lint and declared-test-dependency guards passed.
- Local `rcmdcheck` finished with 0 errors, 0 warnings, and 1 existing installed
  package-size note.
- The final archived coverage run reported 88.50%.

## Conclusion

There is no reason to set 100% as the next target. The suite now provides
substantially better protection of the contracts most likely to affect users
or scientific correctness. Further test work should be triggered by active
behavioral risk, observed regressions, or newly supported integrations rather
than by the remaining percentage alone.
