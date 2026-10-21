# Release and CI policy

[CONTRIBUTING.md](../../CONTRIBUTING.md) owns local development/check commands;
[RELEASING.md](../../RELEASING.md) owns candidate identity and publication gates.
This note records the cross-platform acceptance policy.

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
