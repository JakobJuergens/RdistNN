# Release and CI policy

[CONTRIBUTING.md](../../CONTRIBUTING.md) owns local development/check commands;
[RELEASING.md](../../RELEASING.md) owns candidate identity and publication gates.
This note records the cross-platform acceptance policy.

Follow [version progression](../../RELEASING.md#version-progression) for each
development commit and production merge; read the current version from
`DESCRIPTION` and the current `NEWS.md` heading and keep them synchronized.
The first development commit after production `0.1.x` establishes `0.1.x.9000`;
subsequent commits increment to `.9001`, `.9002`, and so on.
Drafting an uncommitted version does not advance that sequence.

## Merge and release gates

Nontrivial changes need relevant regression coverage, current documentation and
generated files, consistency with design contracts, and required package checks.
Match local checks to the changed surface; use broad suites for shared/high-risk
changes and final release validation. Reuse successful results for unchanged
content rather than repeating checks for staging, commits or pushes.

Before tagging or declaring release readiness, require green R CMD checks on
the supported matrix, a successful pkgdown build, current examples/vignettes,
no generated-file drift and an explanation for unexpected coverage regressions.
The compiled matrix covers Linux x86-64, Windows x86-64, macOS arm64 and macOS
x86-64; Linux exercises R devel, release and oldrel. Jobs verify actual runner
architecture against their labels.

Package-build steps allow two concurrent compiler jobs through `MAKEFLAGS`.
Keep this bound separate from runtime OpenMP/BLAS or future-worker settings.
Compiler optimization and floating-point flags, platform coverage and numerical
acceptance requirements stay the same. Do not cache package test results or
reuse native objects across incompatible compiler/R configurations.
Apply the same compiler-job bound to pkgdown's local installation and coverage's
instrumented installation; these are package builds as well.

## What CI protects

CI checks portability, compiled/fallback parity, generated-file consistency and
executable documentation. Test optional-package paths where valuable and retain
useful behavior when those packages are absent. Prioritize public formulas,
validation, methods, engine parity, simulation defaults and summarizers over a
coverage percentage or defensive assertions without a supported failure mode.

Regenerate namespace, Rd and Rcpp interfaces with their declared tooling;
stale/missing/orphan outputs must fail the drift gate. Pkgdown and vignette
failures are user-facing regressions.

## Source and output hygiene

Local agent settings, paper PDFs, check directories, tarballs, rendered README
HTML, harvested logs and generated archives are not repository source. Keep
them ignored and untracked. Workflows use finite artifact retention and must
not commit logs/archives back to the default branch. The alignment workflow's
repository-hygiene gate enforces candidate-file and ownership rules.

## Developer analysis and CI ownership

Keep fast script correctness checks in the existing alignment workflow. Pin
analyzer versions and select meaningful correctness rules; formatting remains
an explicit file selection, not an automatic repository-wide rewrite. Apply
exclusions to explicit paths as well as recursive discovery, including frozen
scientific source, generated files and private output. See
[Ruff configuration](https://docs.astral.sh/ruff/configuration/).

Use the same coverage runner locally and in the existing manual coverage
workflow. Native Clang, compilation-database, tidy and sanitizer profiles are
explicit developer commands. They complement required portability/release
gates without adding a duplicate platform matrix. Optional analysis tools stay
outside runtime package dependencies. Fail requested checks when tools or
eligible inputs are missing, keep diagnostic exits visible, and retain bounded
failure artifacts. Source archives exclude analyzer configuration and caches.
