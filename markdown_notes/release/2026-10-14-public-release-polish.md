# Initial Public Release Polish Implementation Plan

> For agentic workers: use the executing-plans workflow for this documentation and packaging pass.

**Goal:** Prepare the existing 0.1.0 package for a public introduction centered on CATE and nonparametric regression, with PLQR explicitly experimental.

**Architecture:** Preserve every estimator, export, default, and computational path. Edit documentation sources, regenerate help and README output, and validate an installable source archive and the pkgdown site.

**Tech stack:** R, roxygen2 8.0.0, knitr/rmarkdown, pkgdown, R CMD build/check.

**Spec:** The user's active release-polishing objective and `AGENTS.md`.

## Constraints and review focus

- Work on `master`; preserve any unrelated work. Do not add functionality.
- CATE inference depends on the documented nuisance-rate, stability, overlap, and localization conditions; demonstration nuisances are oracle inputs.
- NP objects expose `mu`/`se`, diagnostics, and plotting; do not promise CATE-specific S3 inference methods for NP.
- Support matrix and data-frame onboarding, with executable deterministic examples.
- Mark PLQR estimators, scores, demo helpers, dedicated methods and plots, and shared PLQR dispatch as experimental.
- Keep generated files generated and local build outputs out of version control and source distributions.
- Prepare website copy and publication instructions; the task is preparation, not a request to publish a release or edit a live website.

## Tasks

- [x] Refocus `README.Rmd`, `vignettes/getting-started.Rmd`, package help, `DESCRIPTION`, and `_pkgdown.yml` on CATE/regression; preserve advanced reference coverage.
- [x] Add experimental lifecycle text to PLQR roxygen sources and PLQR article sections; update the lifecycle policy and live roadmap to match this release scope.
- [x] Add `CONTRIBUTING.md`, `RELEASING.md`, and `WEBSITE.md`, record the scope in `NEWS.md`, and exclude repository-only publication material from the R source package.
- [x] Regenerate help and README, verify generated-file drift and reference coverage, build/check the package with vignettes, build pkgdown, inspect archive contents and rendered onboarding, and verify code expressions match the starting snapshot.

## Verification record

Completed on 2026-10-14 against the local polished source. Publication is a separate step described in `RELEASING.md`.

- `devtools::document()` and `devtools::build_readme()` regenerated help and README using roxygen2 8.0.0.
- `devtools::test(filter = "notation-vignette-contract|readme-scientific-contract|dnn-layer4-grid-engine")`: 121 assertions passed, no failures/warnings/skips.
- `devtools::test(reporter = "summary")`: no failures; one opt-in visual-baseline skip. Two near-uniform-weight diagnostics arose in the existing seven-observation CATE propensity fixtures; they are expected from those inputs, not a changed estimator path.
- `devtools::build(path = "archives", vignettes = TRUE, manual = FALSE)` and `rcmdcheck::rcmdcheck(src, args = "--no-manual", error_on = "warning", check_dir = "archives/check-final-clean")`: 0 errors, 0 warnings, 1 installed-size note (6.7 MB, including 4.4 MB compiled libraries). Installed-package tests passed 4,238 assertions; source-only checks were additionally run locally as above.
- Full pkgdown site built successfully; changed pages were rebuilt after review fixes using `lazy = TRUE`. NEWS now renders as a versioned changelog. The final generated-output cleanup also adds navigation from PLQR methods back to their estimator and removes README output whitespace. A documented pkgdown redirect from `README.html` to `index.html` fixes generated links in the agent entry pages; `pkgdown::build_redirects()` regenerated that final addition.
- HTML inspection: all relative links/assets resolve across 196 pages; the home page, executable onboarding, experimental reference group/method labels, and versioned NEWS have the required content.
- Source archive audit: 385 entries and all 11 built HTML vignettes; excludes private notes, agent settings, CI-only tests, repository-only publication docs, logs, prior archives, and compiled binaries. Archived README matches the generated current README. `archives/SHA256SUMS` records the final archive digest (`f6d06e251c4b8e99608fa4e8b8e0d21c60d7c32b98a1c54b0a2b1619fbdd845b`).
- All 124 public help topics are indexed. Repository hygiene, test dependency declarations, generated-file drift, five Python repository-policy checks, and lint passed. The final plot-tuning and PLQR-method documentation also passes focused lint; `git diff --check` is clean.
- Parsed expressions in all 74 R source files match starting commit `40060f1`; NAMESPACE, Rcpp exports, and C++ sources are unchanged. All 14 exported PLQR S3 registrations map to experimental-labeled help, and all dedicated PLQR topics carry lifecycle sections.
- An independent read-only review found the missing PLQR tuning-method label and inaccessible companion-repository link. Both were resolved; the latter now links to the public research entry on jakobjuergens.com.

Check logs are in `/tmp/rdistnn-release-*.log`; the final R check receipt is `archives/check-final-clean/RdistNN.Rcheck/00check.log`. Local artifacts and this internal record are excluded from the installable package. The full hosted platform matrix must pass on the publication commit before tagging, as required by the release policy.
