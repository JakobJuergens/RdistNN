# Release and CI Policy

This note standardizes merge, release, and CI expectations for `RdistNN`.

## 1. General rule

CI is part of the package contract, not an afterthought.

A code change is not complete just because it works locally. It is complete when
it fits the package’s documented and tested cross-platform behavior.

## 2. Merge standard

Before merging a nontrivial package change, the following should hold:

- package checks pass on the supported CI matrix;
- tests for new behavior are present;
- docs are updated if user-visible behavior changed;
- generated files are up to date;
- no obvious interface/design-note conflicts remain.

## 3. Release standard

Before tagging or calling a state release-ready, ensure:

- R-CMD-check is green on the current matrix;
- pkgdown builds successfully;
- coverage does not regress unexpectedly without reason;
- vignettes and examples reflect the preferred current API;
- no known generated-file drift remains.

## 4. CI as regression protection

CI should guard the package against:
- cross-platform breakage;
- stale generated files;
- fallback/compiled drift when possible;
- accidental documentation rot.

## 5. Generated-file drift checks

Where practical, CI should detect when generated artifacts are stale.

Examples:
- `NAMESPACE` out of sync with roxygen comments;
- `RcppExports.*` out of sync with compiled-interface source.

Generated-file drift should fail fast.

## 6. Coverage expectations

Coverage is a tool, not a fetish.

High-priority areas for coverage:
- exported functions;
- validation branches;
- object methods;
- compiled/fallback equivalence points;
- simulation defaults and summarizers.

Low-value coverage should not crowd out high-value regression tests.

## 7. CI and optional dependencies

Suggested-package paths should be tested where valuable, but the package should
remain robust when optional dependencies are absent.

CI configuration should reflect the intended dependency contract, not a hidden one.

## 8. Documentation CI

Pkgdown and vignette build health should be treated as part of package quality.

If docs no longer build, the package has regressed in a user-facing way.

## 9. Branch review checklist

Before approving a package PR, ask:

- Are the design notes still being followed?
- Are tests present for new validation or interface behavior?
- Are generated files current?
- Does CI cover the risk introduced by this change?
- Would a new user understand the updated behavior from the docs?

## 10. Standard

A stable package workflow is one where:
- code,
- docs,
- generated files,
- tests,
- and CI expectations

move together, not separately.
