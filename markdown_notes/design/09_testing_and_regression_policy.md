# Testing and regression policy

Nontrivial public behavior changes need deterministic regression tests. Cover
valid, edge and invalid inputs and the intended class/shape. Exercise new
validation branches directly, including missing columns, unknown extras,
infeasible deletion orders, malformed metric/weight outputs and invalid ranks.
A bug regression must fail for the original defect at the most stable public
boundary available; use internals only when necessary.

Set seeds explicitly and locally. Reuse domain helpers rather than duplicate
random setup. Test shape, invariants, interface equivalence and known identities,
not incidental execution order. Match consequential errors/warnings on meaningful
substrings unless exact wording is part of the contract. For structured package
conditions, assert identity and useful fields as described below.

Keep the flat `tests/testthat/` layout and assign test blocks to their observable
contract: primitive arithmetic, estimator behavior, input validation, result
methods, retained-state reuse or execution infrastructure. Numerical consumer
regressions may intentionally span these layers. Shared fixtures and independent
references belong in `helper-*.R`; single-suite fixtures stay beside their tests.
Loading helpers must not consume or reseed the caller's RNG stream. Fixture
generators use explicit scoped seeds, and test-only state guards stay independent
of the production restoration code they help test. Keep future-plan/options
cleanup with the backend tests that own those resources.

This applies to suite-local generators and file loading as well as shared
helpers. A file-level seed does not make individual tests independently
reproducible. Use the existing independent test RNG scope for fixture seeds,
including restoration after failure and restoration of an absent seed. See
[R Packages, self-contained tests](https://r-pkgs.org/testing-design.html)
and its [state-restoration guidance](https://r-pkgs.org/code.html).
Ordinary test blocks that seed or consume the process RNG also restore the
incoming state. Testthat's test environment does not do this automatically.
Use standard [withr local seed preservation](https://withr.r-lib.org/reference/with_seed.html)
at the owning block, preserving intentional draw sequences within it. Keep the
independent test-only scope for direct restoration-contract checks.

When splitting suites, preserve reference calculations, assertion strength,
tolerances and native-engine requirements. Update the suite map and complete-family
filters together, then run each moved suite independently. Cross-estimator suites
such as subset shapes and simulation backends remain cohesive; include them in
the relevant family selections even when their filenames have a shared prefix.

## Error and warning contracts

Follow [the reporting design](05_error_messages_and_warnings.md) when
adding or consolidating conditions. Shared constructors need a small reviewed
set of complete-message examples, including expected/actual facts and an
actionable remedy. Consumer tests should assert the condition class/code,
relevant argument/location fields and meaningful remedy text. Do not repeat
complete English snapshots across every estimator or make punctuation an
accidental recovery protocol.

Test common input routes and shared choice validators through public calls,
including the actual selecting argument and bounded public `operation` context.
Retained notice flags/counts must survive serialization and candidate selection,
remain accessible without optional weights, and reject missing or inconsistent
current state. Success-path instrumentation must reject eager message formatting.
Keep opaque callback/native/resource failures distinct from invalid callback
returns; broad error wrapping is not evidence of condition coverage.

Test valid input, invalid input and relevant boundaries. Include missing
optional extras separately from inconsistent retained state. Exercise malformed
extension output and explicit invalid plot options; a convenience fallback
must not hide them. Prove that a recognized optional-panel/candidate failure
recovers, while an unexpected error, interrupt or allocation failure propagates.
The default candidate must not quietly become a missing result.

Cover `options(warn = -1/0/1/2)`, calling handlers and warning muffling, with
options restored after success and error. A warning promoted by `warn = 2`
must not become a failed-candidate row or disappear in a dashboard. Shared test
warning filters should identify only the intended package condition family
once structured conditions exist. Inject an unrelated warning to prove it is
still observable; blanket suppression is not a substitute for this coverage.

Verify one warning per issue family per public call, including nested panels
and grids, and a fresh warning in a separate independent call. Test early
ignored-input warnings, bounded lists with omitted counts, long/unusual names
and dynamic literal text. Mixed-target and all-unavailable ratio results must
preserve every row/status; deliberately disabled uncertainty is not a warning.
For shared Stata methods, cover ordinary, quiet, capture and capture-noisily
execution, exact standard return codes, literal rendering and promised caller
state preservation.

## Group reuse and reporting cost

Prove actual grouping reuse through scalar/batch fits, target reorderings,
grids, predictions, serialized current fits and multi-panel consumers. After
group preparation, reporting must cause no distance evaluation, sorting or
group reconstruction. Instrument those operations rather than merely asserting
that a cache field exists. Counts refer to actual unchanged geometry, not an
incidental helper-call sequence. Include corruption rejected before reuse,
independent calls, and a deletion-refitted geometry whose tie pattern differs
from the full sample.

Check that success paths do not construct messages or allocate per-observation
conditions, and that summaries retain at most bounded example lists. Keep
structural/performance tests with the behavior they protect; share fixtures
and counters when genuinely reused. Benchmarks measure before/after scalar,
batch, grid and reuse cost and allocations. They are review evidence, not
flaky wall-clock CI thresholds or a reason to skip scientific safety checks.

## Numerical assertions

Use explicit absolute/relative tolerances appropriate to quantity and scale.
Ranks, stable tie order, support, target identities, statuses and missingness
require exact equality after documented representation mappings. Where operation
order/rounding is controlled, require exact numeric equality, as in near-tie
distance tests. A small distance error cannot excuse a changed neighbor or
status. Investigate discrepancies before changing tolerances or fixtures.

For retained-state snapshots, long byte/numeric vectors and bit-sensitive
extreme fixtures, use `.expect_identical_fast()` from
`tests/testthat/helper-comparisons.R`. It uses base `identical()` and reports a
bounded failure without recursively diffing closure environments or formatting
every element. Exact values, types, dimensions and attributes remain checked;
this changes diagnostics rather than accepting a numeric tolerance. Ordinary
small assertions can retain testthat's detailed differences.

Test helpers must not catch an evaluator error and substitute the expected
answer, or skip a required compiled path after an unexpected failure. Such
adapters can turn a broken implementation into an apparent reference match.
Skip only a genuinely optional capability; test capability absence separately
from an available engine failing. Use direct supported APIs and independent
references for successful score and distance evaluations.

Required native tests exercise the installed DLL even without local build
tools. Do not skip them because a source file has a newer timestamp; timestamps
do not establish executable identity. Force the portable comparator and count
actual native execution so a native-labelled comparison cannot run the same
engine twice. Preserve native-boundary assertions rather than adding a
stale-install forgiveness adapter.

Namespace call counters use the shared instrumentation helper, which captures
a tracer closure and cleans up on success/error without fixed global bindings
or replacement constructor environments. Its current consumers are non-nested;
see [base R tracing](https://search.r-project.org/R/refmans/base/html/trace.html).
Declare the test dependency floor needed by APIs in use: the existing
`test_file(desc=)` lifecycle checks require testthat 3.2.0, as recorded in its
[release notes](https://testthat.r-lib.org/news/index.html#testthat-320).

Allocation-size regressions may filter `Rprofmem()` below the failure limit
to avoid writing a stack trace for every small allocation. Keep the original
input dimensions and rejection limit, require a nonempty trace, and retain
positive/negative allocation controls showing that both permitted work buffers
and forbidden quadratic buffers remain observable. Do not substitute timing
thresholds for the memory contract.

Test reporting immutability on the fitted object or diagnostic data directly,
not serialized closure/source-cache bytes. Keep genuine save/restore behavior
tests and frozen-file byte checks separate. Never replace exact rank, support,
status or RNG-state checks with approximate comparisons merely to quiet CI.

Keep numerical regressions tied to mathematical contracts. Use independently
demonstrable overflow to test overflow guards, alongside finite controls. Do
not require a mathematically finite result to fail merely because the current
accumulator does. Preserve the original public regression while permitting
valid accuracy improvements; do not encode an implementation defect as the
desired answer. Integrate cases into their owning suites and shared helpers
instead of accumulating copied fixtures and one-off exception tests.

Concentration/status regressions distinguish zero mass from positive mass whose
squares are subnormal or zero, and test safe common rescaling against an
independent reference. Include ordinary-path compatibility and actual
prediction/candidate reuse. Validation tests prove callback nonexecution for
invalid public signals and preserve callback failure visibility for valid
inputs. These are interaction contracts even when their individual source lines
already have coverage; see [R Packages' testing guidance](https://r-pkgs.org/testing-design.html).

## Engine and cross-package parity

Verify both within-package engine parity (C++ versus R, C++ versus Mata in the
companion) and cross-package parity for shared methods. Compare the same numeric
inputs, row order, targets, fitted geometry, signals, weight controls and
covariance definition. Exact and asymptotic weights are distinct methods;
frozen-rank and geometry-refitting deletions also have different contracts.
Engine equivalence applies to implementations of the same method and state.

Force each engine and establish that it actually ran. A compiled parity job
must fail if its required native entry point is unavailable; a portable-only
run may skip native checks explicitly, but must report that limit. Exercise
automatic fallback separately so a nominal native comparison cannot run the
fallback twice. Retain independent analytic or explicit-loop references.

For an affected primitive, cover ordinary inputs and consequential boundaries:
exact/near ties, large common offsets, tiny or constant coordinates, clamped
scales, signed weights, zero variance, missingness and overflow as applicable.
Follow changes through scalar/batch calls, reordered/duplicate targets, tuning
candidates, covariance and saved-fit prediction when those consume the primitive.

Use the same exported numeric data for cross-language fixtures. Equal seed
values in R and Stata do not establish identical random draws. Native portability
fixtures that depend on individual binary digits should use hexadecimal numeric
literals or exact binary constructions: decimal conversion can use different
intermediate precision across R builds. Preserve the intended input bits and
exact discrete assertions rather than relaxing tolerances.
Native portability jobs must run the relevant regressions on their supported architectures;
compiler contraction and parallel reduction changes need targeted coverage.
Record compiler, architecture, math-library and thread settings when they help
explain a discrepancy. The arithmetic contract is in
[the compiled-code policy](10_rcpp_and_compiled_code_policy.md).

## Simulation and replication

Public demo/data, subset and simulation entry points scope explicit `seed` or
`master_seed` values: restore the caller's RNG kind and saved state (or seed
absence) on success and error. `NULL` consumes the current stream. A simulation
replicate seeds before data generation and learner fitting; nested data
calls receive `seed = NULL` so learner/fold draws follow data draws unchanged.
Low-level DGPs implement that stream and do not own a separate restoration scope.
Capture before user callbacks as well. CATE seeds tuning callbacks separately,
then restarts the data/learner stream to preserve its established draws.
R's saved state excludes the Box-Muller
normal cache; the scope does not promise to restore that hidden state.

RNG restoration, repeatability, future-plan agreement and optional-dependency
agreement are separate contracts. Repetition requires the same seed, RNG kind,
dependency route, methods, tuning callbacks and execution configuration. CATE
and IV runners use `future_lapply(..., future.seed = TRUE)` when `future.apply`
is installed, including under a sequential future plan. The replicate bodies
then receive L'Ecuyer-CMRG streams, and their explicit `set.seed()` inherits that
kind. Package-controlled draws agree across sequential and multisession future
plans for a fixed RNG configuration; external learners and callbacks retain
their own reproducibility obligations.

Without `future.apply`, the `lapply` route inherits the caller's RNG kind.
The same integer replicate seeds can therefore yield different data draws
across dependency routes. A caller L'Ecuyer-CMRG kind aligns oracle draws with
fixed deterministic cell orders, but does not establish general dependency
invariance: future mapping advances the caller stream differently, so
stochastic cell-level tuning can select different orders in later cells.
Keep these established streams. Record RNG kind, optional-dependency
availability, package versions and the future plan with a simulation's seeds.
Bounded backend regressions use oracle methods, feasible fixed orders or small
stochastic tuning callbacks, a few replications and at most two workers. Restore
the caller's plan, options, RNG kind and seed on failure as well as success,
and shut down only resources created by the test.

Check seeded determinism, default resolution, return schemas and summarizers.
Small simulation tests establish computation rather than general coverage or
manuscript evidence. Test selection and commands belong in the
[test map](../../tests/README.md), with release gates in
[RELEASING.md](../../RELEASING.md).

## Installed-package and independent-reference coverage

Exercise registered S3 methods through public generics on fitted objects.
A test must not skip an entire method family because an unexported method
symbol is absent from the installed namespace. Similarly, public-interface
tests should not depend on test-only aliases or unrelated internal helpers.
Use explicit test-local reference calculations when independence is needed.

Prioritize score formulas, deterministic demo-data contracts, result reporting,
prediction, candidate views and compiled/fallback equivalence over a coverage
percentage. Oracle examples should agree with independently calculated truth;
a pair of optimized implementations agreeing is not an independent oracle.
Remaining defensive branches warrant tests when an observed defect, supported
input or changed contract makes them consequential. Unused helper coverage
and benchmark timing are not substitutes for user-facing regression tests.

Coverage must exercise the instrumented installed package. Copy current source
inputs before calling covr, because its cleaning/instrumentation modifies that
tree. A focused run uses installed `testthat::test_package()`, not a development
loader that can replace the instrumented namespace or DLL. Reject empty test
selections and propagate test failures. Record positive R hits and positive
handwritten C++ hits when the selection promises native execution, alongside
the actual installed library/DLL identity and source hashes. Full coverage
requires both kinds of execution; a deliberate pure R selection need not claim
native coverage. Report percentages as diagnostics rather than correctness
thresholds. The [covr API](https://covr.r-lib.org/reference/package_coverage.html)
defines the instrumentation lifecycle; commands belong in
[CONTRIBUTING.md](../../CONTRIBUTING.md#developer-analysis-tools).

Developer-tool regressions use isolated temporary fixtures. Cover real nonzero
tool diagnostics, missing tools, empty selections and scope exclusions, with
positive controls proving eligible work actually ran. Keep a small real
instrumentation/failure control in addition to argument and path tests; mocked
command construction alone does not establish a working coverage or sanitizer
pipeline. Optional-tool absence may skip infrastructure integration tests;
an explicitly requested developer command must fail actionably instead.

Native boundary regressions include null/deserialized contexts and wrong context
identity where an external-pointer interface is retained. Run a potentially
terminating case in an isolated R process and require both normal process
completion and the intended caught condition; an in-process error assertion
cannot diagnose a process abort. Include actual serial and OpenMP configurations,
not only requested thread counts. Deterministic native calls also need controls
for both an existing RNG state and an absent `.Random.seed`.
