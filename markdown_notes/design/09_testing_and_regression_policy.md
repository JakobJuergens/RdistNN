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

- compare with tolerance;
- test the property that matters;
- avoid brittle equality unless exact arithmetic is intended.

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

When there is both a compiled and an R fallback path, tests should ensure:

- both paths return compatible results where feasible;
- the fallback remains correct if compiled code is unavailable.

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
