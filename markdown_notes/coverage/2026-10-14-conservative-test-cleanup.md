# Conservative Unit-Test Cleanup Plan and Record

> For agentic workers: implement this plan inline with the executing-plans workflow. Read-only family audits inform decisions; the parent owns all edits and acceptance.

**Goal:** Remove demonstrated redundancy and obsolete development scaffolding while preserving distinct regression, scientific, validation, and fallback contracts.

**Architecture:** Change tests and this record only. Preserve the existing release-polish work. Delete unused test helpers, reuse identical fixture calls and assertions, and retain every independent numerical or boundary scenario.

**Tech stack:** R, testthat, pkgload, R CMD build/check, Python repository-policy checks.

**Spec:** The user's conservative unit-test cleanup objective, `AGENTS.md`, and design policy 09. The July 2026 coverage audit is historical evidence, not an instruction to expand coverage.

## Constraints and review focus

- Do not change package computation, exports, documentation, defaults, or dependencies.
- A fixed bug is not a reason to drop its regression. Removal needs a surviving equivalent contract or proof that a helper/file is unused.
- Keep independent oracles, exact finite-population support, ties, signed weights, moment/coefficient distinctions, malformed-input branches, and compiled/fallback parity.
- Keep genuine optional-dependency, compiled-capability, and opt-in visual skips. Required package APIs should fail tests if missing.
- Preserve exact fixture defaults and inputs when moving assertions or sharing helpers.
- Verify the source and installed test layouts. Compare package/doc hashes against the pre-cleanup snapshot.

## Baseline

The test tree contains 22,450 R lines. A fresh source run recorded 835 cases,
4,276 passed assertions, no failures/errors, two uniform-weight fixture
warnings, and one opt-in visual skip. Results are saved in
`/tmp/rdistnn-test-baseline-results.rds`; tests and non-test hashes are backed
up under `/tmp/rdistnn-conservative-test-cleanup-baseline/`.

## Planned changes and surviving protection

- [x] Delete the zero-byte `test-api-jackknife.R` and unused
  `helper-fastpaths.R`, `helper-registry.R`, and `helper-scores.R`. Their
  definitions have no consumers; active engine checks use the package dispatch
  helpers directly.
- [x] Remove unused `.make_X_helper` and `.make_x0_jk`; remove the duplicate
  `.make_X0` in `test-tdnn-weights.R`, retaining the identical shared definition
  in `helper-jackknife.R`.
- [x] Remove the generated distance battery's first-column-only equality
  test: the next battery test uses identical inputs/context and compares every
  batch column with scalar evaluation. Keep the other numerical distance tests.
- [x] Consolidate NP and TDNN interface metadata checks into their existing
  exact/no-SE interface tests using the same seeded fixture and calls. Collapse
  separate NA/Inf index blocks into the stronger distinct-error case, retaining
  explicit error expectations, message matches, and distinctness.
- [x] Capture each malformed NP outcome error once instead of calling the same
  failing estimator twice. Preserve both regex checks and message distinction.
- [x] Reuse one identical geom extractor across five plotting suites; remove
  one repeated build of the same plot and one immediately overwritten seed.
  Route the PLQR candidate moment regression through `plot_estimates()` and
  retain its numeric assertions against the public plot's data.
- [x] Keep one shared `tuning_mode` signature check instead of three identical
  layer-file checks. Remove the completed-rollout source-filename/grep check:
  actual grid/scalar equivalence and no-public-refit plotting tests protect the
  behavioral contracts in both source and installed layouts.
- [x] Remove the two TDNN repeated-table equality blocks already covered by
  the table-contract battery. They did not establish cache reuse. Keep their
  separate evaluator parity fixtures, all formulas, and independent oracles.
- [x] Remove obsolete required-API availability guards after checking their
  dependencies against the active namespace/source. Keep optional capabilities.
  Correct the flag-only rank-weight fixture so it does not set a class first.
- [x] Capture the two intentional uniform-weight CATE warnings explicitly;
  keep the seven-observation inputs and their uniform-weight mean/quantile
  identities unchanged.
- [x] Run focused tests, full source and installed suites, lint/dependency
  checks, a package check, and an independent review. Record deletions, retained
  coverage, metrics, and non-test preservation evidence.

## Retained after audit

The plotting phase/contract files, layer-0 oracle, scalar ties/custom-distance
suites, default-order guards, scientific documentation checks, simulation and
replication contracts, generalized/TDNN jackknife engines, and compiled/fallback
suites remain. Similar tests with different meaningful inputs or branches are
not treated as duplicates. Python workflow-policy checks remain unchanged.

The audit deliberately declines broad file deletion, consolidation of unlike
plot fixtures, removal of random numerical parity fixtures merely because the
same path appears elsewhere, and a package helper/framework refactor.

## Verification and results

Completed on 2026-10-14. Tests/helper R files went from 92 to 89;
lines went from 22,450 to 21,682, a reduction of 768
lines (3.42%). This is not a target or a claim
about coverage percentages: removals are justified by the surviving contracts
listed above.

| Run | Cases | Passed assertions | Failures/errors | Unhandled warnings | Skips |
|---|---:|---:|---:|---:|---:|
| Before, source | 835 | 4,276 | 0 | 2 | 1 |
| After, source | 818 | 4,255 | 0 | 0 | 1 |
| After, installed package | reported by R CMD check | 4,220 | 0 | 0 | 3 |

The source skip is the existing opt-in visual-baseline suite. Installed skips
are that visual suite plus two source-only documentation checks; the removed
rollout source-grep check no longer contributes an installed-layout skip.
The previously unhandled warnings are now explicit expectations on unchanged
uniform-weight fixtures. Compiled/fallback and optional-package tests retain
their genuine capability guards.

- Focused affected suites passed, followed by `devtools::test(reporter =
  "summary")` across the complete source suite.
- A source build with all 11 vignettes and `rcmdcheck::rcmdcheck(src, args =
  c("--no-manual", "--no-vignettes"), error_on = "warning")` passed with
  0 errors, 0 warnings, and the existing 6.7 MB installed-size note (4.4 MB
  compiled libraries). Vignette inputs were unchanged, so their already
  established runtime checks were not repeated in the package check.
- Two isolated mutation controls failed once each, without harness errors:
  an incorrect flag on an unclassed rank-weight fixture and an incorrect
  expected NP metadata count. This confirms that the corrected flag branch
  and moved metadata assertions discriminate errors. No live package/test
  sources were mutated for these controls.
- An independent review found no substantive coverage loss and independently
  confirmed all 238 baseline package/documentation hashes are unchanged.
  Its minor obsolete optional-integration wording was corrected, as were
  two custom-rank-weight titles that named an API those tests do not call.
- After those label/comment changes, parsing all 89 test/helper files proved
  the executable test bodies identical to the installed-check snapshot.
  Fresh focused runs of the three relabeled suites passed 102 assertions.
  Final whitespace cleanup also preserves those parsed bodies.
- The refreshed installable archive contains byte-identical copies of all
  current 89 R test/helper files, excludes deleted helpers and private/CI
  files, and retains all 11 built HTML vignettes. The last rebuild reused
  those unchanged built vignettes rather than executing them again.
  `archives/SHA256SUMS` records `4876b24df17ec3d0711894a21cdb1d48816cbf7705d2c29f6d3ecf01bd60d12c`.
- Dependency declarations, repository hygiene, reference coverage, all five
  Python policy checks, lint, and `git diff --check` pass. The Python suite
  and every package/documentation file in the baseline remain unchanged.

Verification receipts: `/tmp/rdistnn-test-cleanup-baseline.log`,
`/tmp/rdistnn-test-cleanup-focused.log`, `/tmp/rdistnn-test-cleanup-final.log`,
`/tmp/rdistnn-test-cleanup-mutation.log`, and
`archives/test-cleanup-check/RdistNN.Rcheck/00check.log`. The last package
check and full test runs remain valid for the final bodies; subsequent
changes only refreshed test labels/comments and whitespace, with parsing
and focused execution evidence recorded above. No estimator code, public
API, documentation, dependency, or Python policy test was changed.
