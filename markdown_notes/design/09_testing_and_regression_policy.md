# Testing and Regression Policy

This note standardizes how testing should evolve with the package.

## 1. General rule

Every nontrivial user-visible behavior change should come with tests.

That includes:

- new features;
- bug fixes;
- new validation rules;
- improved error handling;
- changed defaults.

## 2. Minimum test set for new exported behavior

For a new user-facing code path, add tests for:

1. a valid baseline case;
2. at least one edge case;
3. at least one invalid-input case;
4. the intended class / return shape.

## 3. Validation-path testing

When validation is added or improved, tests should exercise the exact branch.

Examples:

- missing required columns;
- unknown `return` entries;
- `s >= n` with jackknife enabled;
- malformed distance outputs;
- malformed batch weight outputs;
- invalid permutations for ranks;
- `k > s`.

Improved messages without tests are easy to regress.

## 4. Regression tests for bugs

Whenever a bug is fixed, add a test that would have failed before the fix.

The test should capture the bug at the most stable level available:

- public API if possible;
- internal helper only if necessary.

## 5. Determinism policy

Tests involving random draws must be deterministic.

Set seeds explicitly and locally.
Do not rely on ambient RNG state.

If helper generators are reused across test files, centralize them.

## 6. Numerical tests

For floating-point results:

- use explicit absolute/relative tolerances appropriate to the quantity and scale;
- test the property that matters;
- require exact equality for discrete outputs such as ranks, stable tie order,
  support, target identities, statuses and missing-value masks after documented
  representation mappings;
- use exact numeric equality where the numerical contract specifies identical
  operations and rounding, as in the distance near-tie regressions.

A small distance difference is consequential if it changes a neighbor or its
weight. Check ranks and downstream decisions separately from numeric distance
agreement. A tolerance must never conceal a changed ordering, support or status.
Investigate an unexpected discrepancy before changing a tolerance or fixture.

## 7. Error and warning tests

Use:

- `expect_error()` for hard failures;
- `expect_warning()` for soft diagnostics.

Match on meaningful substrings, not on full messages unless wording stability is
itself important.

## 8. Test the contract, not the implementation accident

Prefer tests that verify:

- output shape;
- invariants;
- equivalence across interfaces;
- known theoretical identities;
- stable semantics.

Avoid tests that only assert incidental internal sequencing.

## 9. Fast path vs fallback

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
jobs must run the relevant regressions on their supported architectures;
compiler contraction and parallel reduction changes need targeted coverage.
Record compiler, architecture, math-library and thread settings when they help
explain a discrepancy. The arithmetic contract is in
[the compiled-code policy](10_rcpp_and_compiled_code_policy.md).

## 10. Simulation / replication code

Paper-facing simulation helpers should still be tested for:

- deterministic output given a seed;
- correct default resolution;
- valid return schema;
- stable summarizer behavior.

## 11. Review checklist

Before merging, ask:

- does this test fail for the bug or misuse we care about?
- is the test deterministic?
- is the assertion specific enough to catch regression?
- is it narrow enough not to fail for unrelated refactors?

## 12. Standard

Code without tests may be exploratory.
Code merged into the package should not remain exploratory.


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
