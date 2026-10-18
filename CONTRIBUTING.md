# Contributing to RdistNN

The initial public release focuses on CATE, nonparametric regression, and their
supporting DNN weights, tuning, diagnostics, and plots. All PLQR functionality
and local IV-GMM are experimental. Please open an issue before proposing a new
estimator or a change to an exported interface.

## Local setup

Use R 4.1.0 or newer and a C++17 toolchain. OpenMP is optional. Install the
development tools and package dependencies from the repository root:

```r
install.packages(c("devtools", "pkgdown", "rcmdcheck", "yaml"))
devtools::install_deps(dependencies = TRUE)
```

Documentation is generated with **roxygen2 8.0.0**, as declared in
`DESCRIPTION`; use that version to avoid namespace drift. If needed:

```r
install.packages("remotes")  # if not already installed
remotes::install_version("roxygen2", version = "8.0.0", upgrade = "never")
```

## Editing and checking

- Document exported functions in their R source. Regenerate help with
  `devtools::document()`; do not hand-edit `NAMESPACE`, `man/*.Rd`, or
  `RcppExports.*`.
- Edit `README.Rmd`, then run `devtools::build_readme()`.
- Write small, deterministic examples. Mark oracle quantities explicitly and
  keep simulation illustrations distinct from manuscript evidence.
- Preserve matrix and data-frame interfaces, returned fields, and defaults
  during documentation or packaging changes.
- Add focused testthat coverage for changes to behavior or validation. Read
  `AGENTS.md` and the linked design policies for the detailed contracts.

For package behavior changes, run `devtools::test()` and relevant focused
checks. [The test-suite map](https://github.com/JakobJuergens/RdistNN/blob/master/tests/README.md) lists the behavior-based files,
shared helpers, focused filters, and optional visual tests. Before a release,
run the full package and site builds described in
[RELEASING.md](RELEASING.md). Documentation changes should at least execute
their changed examples and verify their links and generated output.

The existing alignment checks can be run from the repository root:

```bash
Rscript tools/ci/check_tests_declared_deps.R
Rscript tools/ci/check_repository_hygiene.R
Rscript tools/ci/check_pkgdown_reference_coverage.R
Rscript tools/ci/check_generated_file_drift.R
```

## Maintainer reference

The [design index](markdown_notes/design/README.md) routes interface,
validation, numerical and engineering contracts. The
[result architecture](markdown_notes/design/13_estimator_and_result_architecture.md)
explains shared generic dispatch and method-specific adapters; the
[paired-ratio contract](markdown_notes/design/17_clate_ratio_inference.md)
records CLATE signal, covariance and numerical-status requirements.
[Package scope](markdown_notes/rdistnn_roadmap.md) distinguishes implemented
capabilities from extensions requiring further justification.

Public maintainer documentation should explain current decisions, scientific
limits and reproducible development procedures. Keep personal source-library
configuration, session instructions, completed task lists and local execution
receipts in private working records. Local `archives/` is ignored and excluded
from R source packages; check that private records are untracked before committing.

## Reporting problems

Open an [issue](https://github.com/JakobJuergens/RdistNN/issues) with a minimal
reproducible example, the expected and actual results, and `sessionInfo()`.
For numerical discrepancies, also report whether they persist with
`options(RdistNN.use_compiled = FALSE)`. That selects supported R algorithms;
it does not remove the compiler requirement for source installation.

For a pull request, describe the concrete change, its scope, and the checks
run. Update `NEWS.md` for user-visible changes. Exclude generated check
directories, source archives, rendered site files, and harvested CI logs from
commits.
