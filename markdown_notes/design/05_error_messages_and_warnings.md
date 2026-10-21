# Error messages and warnings

Write messages for a user who knows the call but not its implementation.
Name arguments verbatim, usually in backticks. Explain the expected constraint,
received value/shape and remedy where useful:

- "`k` must satisfy k <= s (got k=5, s=3)."
- "Delete-1 jackknife requires s < n (got s=100, n=100). Use `se = FALSE` or reduce `s`."
- "Column 'm1_oracle' (specified by `m1_col`) not found in `data`. Available columns: ..."

Avoid generic "invalid input" or "dimension mismatch" messages. Distinguish
wrong type, length/shape, absent columns, missingness, nonfinite values and
semantic violations. Give actual dimensions or lengths when they help.

## Errors and warnings

Error when computation cannot proceed safely or would change meaning or
mislead. Warn when computation remains well-defined, such as ignored optional
extras. Issue ignored-input warnings before expensive work. Use
`call. = FALSE` in user-facing `stop()`/`warning()` unless the call adds useful
context.

Reserve internal-inconsistency wording for genuine invariant failures, through
a distinct internal error path. Explain that the state is inconsistent and may
indicate a bug. Avoid frustrated or speculative messages such as "this should
never happen" for ordinary input errors.

Test consequential messages with `expect_error()`/`expect_warning()` and
meaningful substrings; require exact punctuation only when wording itself is
the contract. Validation order belongs in
[the input contract](03_input_validation_contract.md).
