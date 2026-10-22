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
substrings unless exact wording is part of the contract.

Keep the flat `tests/testthat/` layout and assign test blocks to their observable
contract: primitive arithmetic, estimator behavior, input validation, result
methods, retained-state reuse or execution infrastructure. Numerical consumer
regressions may intentionally span these layers. Shared fixtures and independent
references belong in `helper-*.R`; single-suite fixtures stay beside their tests.
Loading helpers must not consume or reseed the caller's RNG stream. Fixture
generators use explicit scoped seeds, and test-only state guards stay independent
of the production restoration code they help test. Keep future-plan/options
cleanup with the backend tests that own those resources.

When splitting suites, preserve reference calculations, assertion strength,
tolerances and native-engine requirements. Update the suite map and complete-family
filters together, then run each moved suite independently. Cross-estimator suites
such as subset shapes and simulation backends remain cohesive; include them in
the relevant family selections even when their filenames have a shared prefix.

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
