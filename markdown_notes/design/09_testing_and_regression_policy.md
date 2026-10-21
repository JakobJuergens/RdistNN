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

## Numerical assertions

Use explicit absolute/relative tolerances appropriate to quantity and scale.
Ranks, stable tie order, support, target identities, statuses and missingness
require exact equality after documented representation mappings. Where operation
order/rounding is controlled, require exact numeric equality, as in near-tie
distance tests. A small distance error cannot excuse a changed neighbor or
status. Investigate discrepancies before changing tolerances or fixtures.

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
