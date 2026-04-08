# Error Messages and Warnings

This note standardizes user-facing diagnostics across `RdistNN`.

## 1. General rule

Every error and warning should be written for a user who:

- knows the function they called;
- does not know the internal implementation;
- wants to fix the problem quickly.

Messages are part of the API.

## 2. Standard error-message structure

A good error message should answer three questions:

1. What was expected?
2. What was received?
3. What should the user do instead?

Not every message needs all three in full prose, but all three should be
recoverable from the text.

## 3. Prefer specific over generic

Avoid messages like:

- "invalid input"
- "bad state"
- "must be numeric"
- "dimension mismatch"

Prefer messages like:

- "`k` must satisfy k <= s (got k=5, s=3)."
- "Delete-1 jackknife requires s < n (got s=100, n=100). Use `se = FALSE` or reduce `s`."
- "Column 'm1_oracle' (specified by `m1_col`) not found in `data`. Available columns: ..."

## 4. Use argument names verbatim

Whenever possible, refer to the user-facing argument name in backticks.

Examples:

- `` `X0` must have the same number of columns as `X`. ``
- `` `ehat` must have length equal to `y` (n=250); got 249. ``

This makes debugging much easier.

## 5. Include actual received values when informative

For scalar constraints, include the actual value.

Examples:

- `k=7`
- `s=300`
- `n=250`

For shape mismatches, include actual dimensions.

Examples:

- expected `100 x 3`, got `100 x 2`
- expected length `n=50`, got length `48`

## 6. Separate categories of failure

Do not conflate:

- wrong type;
- wrong length;
- wrong shape;
- missing column;
- NA values;
- non-finite values;
- semantic constraint violation.

Each category should get its own message.

## 7. Warnings vs errors

### Use an error when:
- the function cannot proceed safely;
- proceeding would silently change meaning;
- the result would likely be misleading.

### Use a warning when:
- the function can proceed safely;
- the input includes ignored or unused extras;
- the behavior is still well-defined and unsurprising.

## 8. Early warnings

Warnings about ignored inputs should happen at the call site, not late in the
pipeline.

Example:
unknown `return` fields should be warned about before computation begins.

## 9. `call. = FALSE`

Use `call. = FALSE` in user-facing `stop()` / `warning()` calls unless the call
itself adds genuine debugging value.

Default preference:
clean message first.

## 10. No passive-aggressive diagnostics

Avoid messages that sound like internal frustration.

Do not write:

- "this should never happen"
- "unexpected object"
- "internal bug" unless it truly is an internal bug

If a message is for developers, route it through a dedicated internal error path.

## 11. Internal invariant failures

For true internal impossible states, use a distinct internal helper such as
`.abort_bad_state()` with wording that signals:

- this is an internal inconsistency;
- the user may have found a bug;
- the package state is not trustworthy.

These should be rare.

## 12. Test messages that matter

When an error message is part of the intended UX, test it with `expect_error()`
using a meaningful substring.

Do not overfit to exact punctuation unless the exact wording matters.
