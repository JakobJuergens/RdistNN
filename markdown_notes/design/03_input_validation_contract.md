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

NP, TDNN and CATE validate supplied outcome/pseudo-outcome missingness and
finiteness before fitted geometry setup or evaluation, including consumed CATE
nuisance columns with a finite-value contract. Unused nuisance inputs on a
supplied-pseudo branch remain ignored. A rejected signal
must not invoke distance callbacks. Keep the direct or column-selecting
argument in the condition. This public-boundary rule does not change score
primitives that deliberately validate contributions at evaluation time.

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
controls, candidate/target identities, the selected tie policy, cached rank
permutations and per-target distance groups; incomplete or inconsistent saved
state requires refitting.

## Messages and tests

Do not use `stopifnot()` on user-facing paths. Use explicit argument-specific
errors explaining what was expected, what was received and how to proceed;
see [the message policy](05_error_messages_and_warnings.md).

The shared argument helpers signal `rdistnn_input_error` with bounded argument,
expected and actual facts. At retained-state owners, reuse these validators
through `.dnn_validate_retained()` so malformed cached arrays become
`rdistnn_state_error` with refit guidance. Translate only owned input-validation
conditions there; do not relabel arbitrary errors, interruptions or resource
failures as bad input. Preserve the bounded argument, expected and actual facts
when translating an owned input condition into a retained-state condition. Missing optional extras have their own identity and
exact retention option, distinct from incomplete required scientific state.

Retained field names are exact schema keys. Check required names before using
`$`, or extract with `[[name, exact = TRUE]]`; a similarly prefixed field must
not replace a missing current field. This applies
[base R extraction semantics](https://stat.ethz.ch/R-manual/R-patched/library/base/html/Extract.html),
where `$` permits partial matching. Current fitted distance specifications retain
their `rdistnn_distance` identity; fit-time structural extension inputs do not
authorize upgrading an older retained specification during reporting.

The retained-matrix accessor owns optional `W`/`dx` shape and finite-value
checks. Missing optional `dx` remains unavailable, while a malformed present
matrix raises a state error with refit guidance. Raw numerical helper arguments
still use ordinary input-validation errors.

Each new public validation rule needs a valid case, an invalid case exercising
the branch and relevant boundary cases. Reuse small named validators for
columns, vector lengths, matrix shapes, binary values and permutations.
