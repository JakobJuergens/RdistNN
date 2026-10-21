# Contributing to RdistNN

Read [AGENTS.md](https://JakobJuergens.github.io/RdistNN/AGENTS.md)
before editing. The [design
index](https://github.com/JakobJuergens/RdistNN/blob/master/markdown_notes/design/README.md)
routes interface, validation and numerical contracts; the [package
scope](https://github.com/JakobJuergens/RdistNN/blob/master/markdown_notes/rdistnn_roadmap.md)
distinguishes supported and experimental workflows. Open an issue before
proposing a new estimator or exported-interface change.

## Development and production workflow

Work on local `development` and keep it unpublished until publication is
authorized. Branches in a public GitHub repository are public when
pushed. For an authorized production update, publish a review branch,
open a PR to `master`, complete review and required checks, then
squash-merge. Direct pushes to `master` require an explicit override.
Keep development current with the reviewed production state.

## Local setup

Use R 4.1.0 or newer and a C++17 toolchain; OpenMP is optional. From the
root:

``` r

install.packages(c("devtools", "pkgdown", "rcmdcheck", "yaml"))
devtools::install_deps(dependencies = TRUE)
```

Use **roxygen2 8.0.0**, as declared in `DESCRIPTION`, to avoid
generated-file drift. Install the declared version if needed:

``` r

install.packages("remotes")
remotes::install_version("roxygen2", version = "8.0.0", upgrade = "never")
```

## Editing and checking

| Changed surface | Edit and verify |
|----|----|
| Function help or exports | Edit roxygen in `R/`; run `devtools::document()`. Never hand-edit `NAMESPACE`, `man/*.Rd` or `RcppExports.*`. |
| README | Edit `README.Rmd`; run `devtools::build_readme()`. |
| Examples or articles | Use small seeded examples, guard optional packages, execute changed code and check links/generated output. Refer to the paper for theory. |
| Behavior or validation | Add focused testthat regressions. Select suites using the [test map](https://github.com/JakobJuergens/RdistNN/blob/master/tests/README.md); use the full suite for shared or broad changes. |
| Release candidate | Run the full package/site and cross-platform gates in [RELEASING.md](https://JakobJuergens.github.io/RdistNN/RELEASING.md). |

Preserve public inputs, defaults and result schemas during editorial
changes. Identify oracle quantities and distinguish demonstrations from
manuscript evidence. Update `NEWS.md` for user-visible changes.

Run the alignment checks relevant to the changed surface:

``` sh
Rscript tools/ci/check_tests_declared_deps.R
Rscript tools/ci/check_repository_hygiene.R
Rscript tools/ci/check_pkgdown_reference_coverage.R
Rscript tools/ci/check_generated_file_drift.R
python3 -m unittest discover -s tests/ci -v
```

The generated-file gate regenerates namespace, Rd topics and Rcpp glue
in an isolated tree using the declared roxygen2 version. Hygiene checks
staged/new candidate files and duplicate top-level R function ownership;
local closures may reuse names. Keep check output, archives and private
session records under ignored `archives/`, excluded from source
packages. Reuse successful checks while the verified content is
unchanged. Reserve `build/` for R-generated package metadata, including
the vignette index.

## R and Stata parity

R is the numerical reference for shared methods. Compare identical
inputs, row order, targets, fitted geometry, signals, weight controls
and covariance definitions. Require exact discrete decisions and
justified numeric tolerances; force and establish each engine, then test
automatic fallback separately. Keep independent analytic/deletion
references. Do not regenerate golden values merely to make a failure
pass. The [compiled-code
policy](https://github.com/JakobJuergens/RdistNN/blob/master/markdown_notes/design/10_rcpp_and_compiled_code_policy.md)
owns the arithmetic contract.

The Stata checkout’s `tests/R_TEST_CROSSWALK.md` maps shared and R-only
contracts. Its `tools/generate_r_*.R` scripts own cross-language
fixtures, including reporting and saved-fit references. Run matching
Stata regressions after an intentional shared-contract change. PLQR
parity concerns supplied-coefficient moments; Stata has no
density-backed coefficient inference.

With both checkouts present, verify mirrored native headers from the R
root:

``` sh
python3 ../RdistNN_stata/tools/ci/package_hygiene.py --r-root .
```

## Reporting problems

Open a [GitHub issue](https://github.com/JakobJuergens/RdistNN/issues)
with a minimal example, expected/actual results and
[`sessionInfo()`](https://rdrr.io/r/utils/sessionInfo.html). For
numerical discrepancies, report whether they persist with
`options(RdistNN.use_compiled = FALSE)`. This selects R algorithms;
source installation still requires compilation.

PRs should describe the final behavior, scope and verification.
Maintainer documentation records current contracts and procedures; keep
personal paths, session instructions, completed task lists, generated
logs and execution receipts in private working records.
