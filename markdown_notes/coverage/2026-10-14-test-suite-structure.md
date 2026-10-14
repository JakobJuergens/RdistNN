# Test-suite structure before the public release

Status: implemented and verified locally; 2026-10-14. Not committed or pushed.

## Scope

Organize the existing tests by behavior while preserving their assertions,
fixtures, numerical tolerances, and optional-dependency guards. This follows
the conservative cleanup recorded in
[2026-10-14-conservative-test-cleanup.md](2026-10-14-conservative-test-cleanup.md).
Earlier release edits remain outside this pass.

## Implementation

- Replace rollout-stage filenames and private fixture prefixes with names for
  scalar references, prefix jackknife, tuning grids, candidate views, and
  diagnostic plots. Keep test files directly in `tests/testthat/`.
- Split the mixed score API, CATE simulation, PLQR wrapper, and GMM wrapper
  files along their existing data, estimation, validation, and method
  contracts. Move fixtures used by the resulting files to domain helpers.
- Add `tests/README.md` with the domain map, focused commands, optional-test
  controls, and fixture guidance; route contributors there from
  `CONTRIBUTING.md`.

## Acceptance

- Compare parsed expressions against the post-cleanup baseline, allowing
  only private identifier renames and relocation.
- Check the split files separately to expose dependencies on test ordering.
- Run the full source suite and an installed-package check of a refreshed
  source archive. Preserve the eleven already-built vignettes because no
  vignette inputs change in this pass.
- Verify unrelated package and documentation content against the baseline.

The baseline contains 89 R test/helper files and records 818 cases, 4,255
passing assertions, no failures or warnings, and one opt-in visual skip in
the complete source run. Baseline files and parsed expressions are stored
under `/tmp/rdistnn-test-structure-baseline/` for this session.

## Resulting layout

Four mixed suites now have the following boundaries:

| Former file | Current files in `tests/testthat/` |
| --- | --- |
| `test-api-scores.R` | `test-api-scores.R`, `test-api-scores-diff.R` |
| `test-sim-cate.R` | `test-sim-cate-data.R`, `test-sim-cate-one-rep.R`, `test-sim-cate-run.R`, `test-sim-cate-summary.R`; shared references in `helper-sim-cate.R` |
| `test-user-dnn-gmm.R` | `test-user-dnn-gmm.R`, `test-user-dnn-gmm-data-interface.R`, `test-user-dnn-gmm-validation.R`, `test-user-gmm-methods.R`; model fixtures in `helper-user-dnn-gmm.R` |
| `test-user-plqr-dnn.R` | `test-user-plqr-dnn.R`, `test-user-plqr-dnn-data-interface.R`, `test-user-plqr-dnn-validation.R`, `test-user-plqr-score-matrix.R`, `test-demo-data-plqr.R` |

Fourteen stage-based files were renamed. The scalar oracle helper is now
`helper-dnn-reference.R`; DNN scalar and prefix tests are
`test-dnn-scalar-reference.R` and `test-dnn-jackknife-prefix.R`. Candidate
schema, grid engine, candidate views, and grid diagnostic plots now use
`test-tuning-candidate-schema.R`, `test-tuning-grid-engine.R`,
`test-candidate-fit.R`, and `test-tuning-grid-diagnostics-plots.R`.
Estimator grid tests live beside their corresponding user-wrapper families;
the simulation grid test is `test-sim-cate-tuning-grid.R`. Plotting phase
files are now `test-plot-diagnostics.R`,
`test-plot-diagnostics-contracts.R`, and `test-plot-specialized.R`.
Private fixture prefixes also describe these contracts. Historical audit
records retain their original filenames as evidence of those snapshots.

The final tree has 102 R test/helper files and 21,683 lines, compared with
89 files and 21,682 lines before this pass. The extra files separate
contracts; this pass adds no assertions or estimator behavior.
[tests/README.md](../../tests/README.md) provides the current maintenance map.

## Verification receipts

- Parsed conservation: all 946 top-level expressions, including 788 literal
  test blocks and all reference helpers, match the baseline after applying
  the documented private-prefix renames. No assertion, tolerance, fixture
  parameter, seed, or optional guard changed.
- Full source run: 818 cases, 4,255 passing assertions, zero failures,
  errors, or warnings; one existing opt-in visual skip.
- Independent runs: all 28 renamed or split test files pass in fresh R
  processes, totaling 218 cases and 1,433 passing assertions, with no
  failures, warnings, or skips.
- Installed-package check: zero errors, zero warnings, one existing size
  note (6.7 MB installed, including 4.4 MB of libraries). Installed tests
  report 4,220 passing assertions and three expected skips: visual
  snapshots and two documentation tests requiring a source checkout.
- The test dependency declaration and repository hygiene checks pass.
  New local documentation links and whitespace checks pass. Package code,
  generated help, vignettes, and earlier release content remain unchanged;
  the only existing non-test files edited here are `CONTRIBUTING.md` and
  `.Rbuildignore` (which excludes the repository-only test map).
- Independent read-only review found no substantive issue. Its suggestion
  to include simulation tuning tests in the documented tuning filter was
  applied. That excluded Markdown file does not change the checked package.
- The source archive retains eleven built HTML vignettes and exactly the
  current 102 test/support files, checked byte for byte. The new archive
  supersedes the earlier cleanup archive, whose receipt remains historical.
  SHA256: `82190dd2cb099366dbd77501c2a6b582e6bec2826430c552099cbd8ac237d1b2`.

Session receipts are `/tmp/rdistnn-test-structure-full.log`,
`/tmp/rdistnn-test-structure-isolated.log`,
`/tmp/rdistnn-test-structure-archive.log`, and
`/tmp/rdistnn-test-structure-check.log`. The persistent check log is
`archives/test-structure-check/RdistNN.Rcheck/00check.log` (ignored).
The expression comparator and move manifest are
`/tmp/rdistnn-verify-test-structure.R` and
`/tmp/rdistnn-test-structure-mapping.json`.

No website rebuild is needed for this pass: the changed map and contributor
guidance are repository-only documents, and no package/site inputs changed
apart from tests and the source-archive exclusion for that map.
