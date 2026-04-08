# Input Validation Contract

This note standardizes how input validation should work across `RdistNN`.

## 1. General rule

Every exported user-facing function must validate:

- required arguments;
- argument type;
- length / shape;
- missingness / finiteness;
- value constraints;
- cross-argument consistency.

Validation should happen before heavy work begins.

## 2. Validation order

Use the following order whenever possible:

1. check mutually exclusive / jointly required arguments;
2. check container type (`data.frame`, matrix, vector, spec object, etc.);
3. check presence of referenced columns or named fields;
4. check shape and length;
5. check finiteness / NA constraints;
6. check semantic constraints (`k <= s`, `s < n` for delete-1 jackknife, etc.);
7. normalize or coerce only after the input is known to be admissible.

This ordering usually produces the clearest user experience.

## 3. User-facing functions vs internal helpers

### Exported functions
Exported functions must do explicit validation.
Do not rely on lower-level code to catch obvious user mistakes.

### Internal helpers
Internal helpers may assume validated inputs only if:

- they are never called directly from external interfaces; and
- the assumptions are obvious from the call path.

If a helper is reused in multiple places, prefer making it robust rather than
relying on fragile upstream assumptions.

## 4. No `stopifnot()` in user-facing paths

Do not use `stopifnot()` in exported functions or in helpers that are part of a
user-facing path.

Use explicit checks with explicit messages instead.

Reason:
`stopifnot()` is concise for internal debugging, but its failure messages are too
opaque for package users.

## 5. Column-based validation

Whenever a function accepts column names through `data = ...` interfaces:

- first check that the named column exists;
- then check that the extracted vector has the right length and type.

Missing-column errors should name:

- the missing column;
- the argument that referenced it; and
- the available columns in `data`.

Column extraction should be centralized in internal helpers where possible.

## 6. Shape validation conventions

### Vectors
When a vector is expected, validate:

- numeric vs other type;
- exact required length;
- finiteness if required;
- special structure if required (binary, permutation, probabilities, etc.).

### Matrices
When a matrix is expected, validate:

- matrix-ness;
- numeric type;
- exact dimensions;
- finiteness.

When rejecting a matrix, report both the expected shape and the actual shape.

## 7. Missingness and finiteness

The package should distinguish between:

- missing values (`NA`);
- non-finite values (`Inf`, `-Inf`, `NaN`);
- shape mismatch;
- type mismatch.

Do not collapse all of these into one generic error.

## 8. Semantic guardrails

Semantic constraints should be checked explicitly and early.

Examples:

- `k <= s`
- `1 <= s <= n`
- `s < n` when delete-1 jackknife is requested
- treatment indicators are binary
- rank vectors are permutations of `1:n`

These conditions are part of the API contract, not just implementation details.

## 9. Coercion policy

Coercion is allowed only when it is safe and unsurprising.

Acceptable examples:

- integer-like numeric scalar to integer;
- vector to matrix via a well-defined helper;
- data-frame column extraction to validated vector.

Avoid coercions that silently change meaning.

If coercion changes interpretation, error instead.

## 10. Validation helpers

Prefer small helpers for repeated validation patterns.

Examples of appropriate helpers:

- column extractor + validator
- vector-of-length-n validator
- matrix-of-shape-n-by-m validator
- binary indicator validator
- permutation validator

Helper names should make their contract obvious.

## 11. Failure mode standard

A validation failure should make it obvious:

- what was expected;
- what was received;
- which argument is at fault;
- what the user can do next, when relevant.

## 12. Testing requirement

For every new validation rule in exported behavior, add tests for:

- one valid case;
- one invalid case;
- edge cases near the boundary when relevant.

Validation that is not tested is not stable enough to rely on.
