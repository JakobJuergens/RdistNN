# Input validation contract

Validate required arguments, types, lengths/shapes, missingness, finiteness,
value constraints and cross-argument consistency before heavy computation.
Demo generators validate target inputs before consuming RNG state. Explicit
seeds are scoped as described in the [testing contract](09_testing_and_regression_policy.md).
Prefer this order:

1. Mutually exclusive or jointly required arguments.
2. Container types and required columns/fields.
3. Dimensions and lengths.
4. Missing/nonfinite values.
5. Semantic constraints.
6. Safe normalization or coercion.

Exports must catch obvious user errors rather than relying on native code.
Internal helpers may assume validated inputs only when they cannot be called
externally and the call path makes those assumptions explicit. Reused helpers
should remain robust at their own boundary.

## Columns and shapes

Use shared single-/multi-column extractors. Check column existence before
extraction, then validate numeric type and shape. Missing-column messages name
the column, selecting argument and available columns. Matrix errors report
expected and actual dimensions; vector errors report required/actual length.
Distinguish wrong type, missing values, nonfinite values and shape errors.

Check semantic constraints explicitly: `1 <= k <= s <= n`, `s < n` for
active delete-one inference, binary indicators and rank permutations of `1:n`.
Apply the constraints of the selected estimator rather than a generic coercion.
Accept integer-like scalars or vector-to-matrix conversion only when meaning
is preserved. Reject coercions that would change interpretation.

Plot target selectors are integer-valued numeric indices, with duplicates
removed in first-occurrence order. Reject fractional, logical, complex and
character selectors before indexing. Candidate reconstruction validates retained
controls, candidate/target identities and cached rank permutations; incomplete
or inconsistent saved state requires refitting.

## Messages and tests

Do not use `stopifnot()` on user-facing paths. Use explicit argument-specific
errors explaining what was expected, what was received and how to proceed;
see [the message policy](05_error_messages_and_warnings.md).

Each new public validation rule needs a valid case, an invalid case exercising
the branch and relevant boundary cases. Reuse small named validators for
columns, vector lengths, matrix shapes, binary values and permutations.
