# Contributing to RdistNN

Read [AGENTS.md](AGENTS.md) before editing. The
[design index](https://github.com/JakobJuergens/RdistNN/blob/master/markdown_notes/design/README.md)
routes interface, validation and numerical contracts; the
[package scope](https://github.com/JakobJuergens/RdistNN/blob/master/markdown_notes/rdistnn_roadmap.md)
distinguishes supported and experimental workflows. Open an issue before
proposing a new estimator or exported-interface change.

## Development and production workflow

Work on local `development` and keep it unpublished until publication is
authorized. Branches in a public GitHub repository are public when pushed.
For an authorized production update, publish a review branch, open a PR to
`master`, complete review and required checks, then squash-merge. Direct pushes
to `master` require an explicit override. Keep development current with the
reviewed production state.

The first development commit after production `0.1.x` establishes
`0.1.x.9000`. Each subsequent development commit increments the fourth component,
including small changes: `0.1.x.9001`, `0.1.x.9002`, and so on. Uncommitted drafts
do not advance this sequence. Read the current version from `DESCRIPTION` and
the current `NEWS.md` heading; update both together. By default, each production
merge becomes `0.1.(x+1)`; the maintainer may choose another version. After that
release, append `.9000` to the selected production version to start development.
For example, the cycle after production `0.1.2` starts at `0.1.2.9000`, then
`0.1.2.9001`; its next default production version is `0.1.3`.

## Current development compatibility

Use the current documented interfaces and complete current fit schemas.
During this early development stage, remove compatibility layers for older
package versions: retired aliases and wrappers, saved-state upgraders and
migration shims. Keep one current implementation unless the maintainer
explicitly requests an exception. Older, incomplete or inconsistent fits
require an actionable refit error; do not silently reconstruct or repair them.
Current fits must still survive serialization and support documented reuse.

Document intentional interface changes and update their regression tests.
Preserve historical changelogs, released tags/assets and frozen scientific
evidence. Compatibility and deprecation commitments can be revisited when the
packages stabilize; prior production releases do not change the current policy.

See the [lifecycle policy](https://github.com/JakobJuergens/RdistNN/blob/master/markdown_notes/design/07_lifecycle_and_deprecation_policy.md) for the package lifecycle and release context.

## Local setup

Use R 4.1.0 or newer and a C++17 toolchain; OpenMP is optional. From the root:

```r
install.packages(c("devtools", "pkgdown", "rcmdcheck", "yaml"))
devtools::install_deps(dependencies = TRUE)
```

Use **roxygen2 8.0.0**, as declared in `DESCRIPTION`, to avoid generated-file
drift. Install the declared version if needed:

```r
install.packages("remotes")
remotes::install_version("roxygen2", version = "8.0.0", upgrade = "never")
```

## Editing and checking

| Changed surface | Edit and verify |
| --- | --- |
| Function help or exports | Edit roxygen in `R/`; run `devtools::document()`. Never hand-edit `NAMESPACE`, `man/*.Rd` or `RcppExports.*`. |
| README | Edit `README.Rmd`; run `devtools::build_readme()`. |
| Examples or articles | Use small seeded examples, guard optional packages, execute changed code and check links/generated output. Refer to the paper for theory. |
| Behavior or validation | Add focused testthat regressions. Select suites using the [test map](https://github.com/JakobJuergens/RdistNN/blob/master/tests/README.md); use the full suite for shared or broad changes. |
| Release candidate | Run the full package/site and cross-platform gates in [RELEASING.md](RELEASING.md). |

Preserve public inputs, defaults and result schemas during editorial changes.
Identify oracle quantities and distinguish demonstrations from manuscript
evidence. Update `NEWS.md` for user-visible changes.

Run the alignment checks relevant to the changed surface:

```sh
Rscript tools/ci/check_tests_declared_deps.R
Rscript tools/ci/check_repository_hygiene.R
Rscript tools/ci/check_pkgdown_reference_coverage.R
Rscript tools/ci/check_generated_file_drift.R
python3 -m unittest discover -s tests/ci -v
```

The generated-file gate regenerates namespace, Rd topics and Rcpp glue in an
isolated tree using the declared roxygen2 version. Hygiene checks staged/new
candidate files and duplicate top-level R function ownership; local closures
may reuse names. Keep check output, archives and private session records under
ignored `archives/`, excluded from source packages. Reuse successful checks
while the verified content is unchanged.
Reserve `build/` for R-generated package metadata, including the vignette index.

## Developer analysis tools

These optional checkout tools do not change package installation requirements.
On Ubuntu/WSL, the native profiles use Clang 18 (`clang++-18`,
`clang-tidy-18`, `llvm-symbolizer-18`), LLVM OpenMP, Bear, ccache and Make.
They also need Rcpp and testthat in the selected R installation. Run from this
repository root:

```sh
python3 tools/dev/cpp_tooling.py doctor
python3 tools/dev/cpp_tooling.py build
python3 tools/dev/cpp_tooling.py compdb
python3 tools/dev/cpp_tooling.py tidy
python3 tools/dev/cpp_tooling.py sanitize --sanitizer-smoke
```

`build` installs a fresh private source copy and runs focused installed native
and portable-parity tests. `compdb` repeats that build with caching disabled
under Bear, then retains the source and `compile_commands.json` for `tidy`.
Run `compdb` again after native inputs, tooling, compiler or R changes; an
ordinary `build` invalidates the previous database. Analysis covers handwritten
translation units and headers; generated Rcpp glue still enters the build and
database. Use `--jobs 1` to reduce compiler concurrency; the maximum is two.

`tidy` also verifies the retained native/configuration bytes and database against
their successful capture. A changed retained copy, edited database or incomplete
receipt requires a fresh `compdb`; analysis never silently repairs that proof.

The default `clang` profile uses LLVM OpenMP; `--profile clang-serial` explicitly
disables it and `--profile gcc` selects the normal GCC comparison. Sanitizers
require a Clang profile. These Linux profiles complement the release matrix.

Logs, private installations, source hashes and command provenance stay under
`archives/cpp-tooling/PROFILE/`. The project cache lives under
`archives/cpp-tooling/cache/`, is bounded to 500 MB per compiler/R/profile key,
and caches compilation only. Tests always rerun. No global Makevars or compiler
settings are edited. AddressSanitizer and UndefinedBehaviorSanitizer instrument
the package; the test process preloads the matching Clang runtime. The host R
binary is not instrumented, and leak detection is disabled for that host.
`--sanitizer-smoke` proves detection with an intentionally invalid private
fixture; its failure must be recognized before package checks can succeed.

Coverage uses a separate source copy because covr cleans and instruments its
input. Install `covr` and `xml2`, then use the shared local/manual-CI runner:

```sh
Rscript tools/ci/run_coverage.R
Rscript tools/ci/run_coverage.R --filter '^(native-threads|compiled-fallback-parity)$' --require-native
```

The default runs the full installed suite and requires actual R and owned C++
coverage. A filter runs installed testthat files matching that expression;
add `--require-native` when the selection promises native execution. Pure R
selections may omit it. Reports, source and instrumented installation remain
under a fresh `archives/coverage/` directory. RDS, CSV, JSON and Cobertura XML
include execution/source provenance. `--html` additionally requires DT and
htmltools. `--output-dir`, `--workspace`, `--install-path` and `--cobertura`
accept fresh destinations; existing outputs are rejected rather than replaced.
Do not run a pkgload/devtools test helper inside coverage: it can replace the
instrumented installed package. See the [coverage test guidance](tests/README.md#coverage).

Python analysis requires Ruff 0.17.2; the configuration selects correctness
rules rather than enforcing a repository-wide formatting conversion:

```sh
python3 tools/ci/check_scripts.py
python3 tools/ci/check_scripts.py tools/dev/cpp_tooling.py
```

Install the pin with `python3 -m pip install ruff==0.17.2` in the Python
environment you use for maintenance. Formatting is opt-in: `--format FILE.py`
changes only explicitly selected files. Both formatting modes reject directories.
Lint defaults to maintained `tools/` and `tests/` Python files. Generated,
vendor, private output and frozen replication sources stay excluded, including
explicit selections. The runner never installs tools, auto-fixes lint or uses
cached test results. Existing static CI runs this lint; native profiles remain
explicit commands and coverage retains its manual workflow.

## R and Stata parity

R is the numerical reference for shared methods. Compare identical inputs,
row order, targets, fitted geometry, signals, weight controls and covariance
definitions. Require exact discrete decisions and justified numeric tolerances;
force and establish each engine, then test automatic fallback separately.
High-level distance ties require exact represented selected-metric groups
per target, rank-weight averaging over observations and actual-survivor deletion
references. Bare supplied orders remain strict unless distances accompany them.
Keep independent analytic/deletion references. Do not regenerate golden values
merely to make a failure pass. The
[compiled-code policy](https://github.com/JakobJuergens/RdistNN/blob/master/markdown_notes/design/10_rcpp_and_compiled_code_policy.md)
owns the arithmetic contract.

The Stata checkout's `tests/R_TEST_CROSSWALK.md` maps shared and R-only contracts.
Its `tools/generate_r_*.R` scripts own cross-language fixtures, including
reporting and saved-fit references. Run matching Stata regressions after an
intentional shared-contract change. PLQR parity concerns supplied-coefficient
moments; Stata has no density-backed coefficient inference.

With both checkouts present, verify mirrored native headers from the R root:

```sh
python3 ../RdistNN_stata/tools/ci/package_hygiene.py --r-root .
```

## Reporting problems

Open a [GitHub issue](https://github.com/JakobJuergens/RdistNN/issues) with a
minimal example, expected/actual results and `sessionInfo()`. For numerical
discrepancies, report whether they persist with
`options(RdistNN.use_compiled = FALSE)`. This selects R algorithms; source
installation still requires compilation.

PRs should describe the final behavior, scope and verification. Maintainer
documentation records current contracts and procedures; keep personal paths,
session instructions, completed task lists, generated logs and execution
receipts in private working records.

## Checkout benchmarks

Install the development-only `bench` package and the current RdistNN, then
run from the repository root:

```sh
Rscript tools/bench/run_bench.R l1 n=5000L reps=3L
Rscript tools/bench/run_bench.R l2 n=5000L reps=3L
Rscript tools/bench/run_bench.R cate n=1000L s=8L n_subsets=128L reps=3L
Rscript tools/bench/run_bench.R jackknife n=200L s=8L reps=3L
```

The same runner accepts the experimental `plqr` score benchmark. Shared
scaffolding belongs under `tools/bench/`, outside installed code. Score runs
compare current batch/scalar APIs and native reductions with explicit subset
sums. Distance runs compare current public providers, R and native distances.
Reporting benchmarks use `tools/bench/bench_reporting.R`; counters include a
positive control before measurement. Timings are review evidence.
