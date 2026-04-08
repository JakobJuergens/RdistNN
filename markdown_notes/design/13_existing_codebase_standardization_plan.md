# Existing Codebase Standardization Plan

This note translates the design rules into concrete instructions for the current
`RdistNN` code base.

The goal is not a rewrite. The goal is a staged cleanup that standardizes the
existing code around a stable interface and validation contract.

## Phase 0: Treat this as a consistency pass, not a feature pass

Before adding more API surface, standardize the current package along the
dimensions below:

- input boundary
- validation style
- error messages
- comments
- tests
- return behavior
- namespace discipline
- compiled/fallback parity
- generated-file hygiene

Do not mix this work with unrelated feature expansion if it can be avoided.

---

## Phase 1: Standardize the input boundary

### 1.1 Establish a single repo-wide rule
Adopt the following rule everywhere:

> User-facing wrappers may accept either explicit matrices/vectors or
> `data.frame`-style inputs, including tibbles. After validation and interface
> resolution, all internal computation should proceed on canonical numeric
> matrices and vectors.

### 1.2 Audit exported functions
For each exported function, classify it as one of:

- high-level wrapper
- lower-level public primitive
- simulation / replication helper

Then enforce the intended input contract.

#### High-level wrappers
Examples:
- `np_dnn()`
- `cate_dnn()`
- `plqr_dnn()`

Instruction:
- allow matrix/vector interface
- allow `data.frame`/tibble interface where intended
- immediately resolve to canonical matrix/vector forms

#### Lower-level public functions
Examples:
- weight builders
- rank / score / distance evaluators

Instruction:
- prefer explicit matrices and vectors
- do not add `data=` interfaces unless there is a strong reason

### 1.3 Add a matrix-extraction helper for multi-column data-frame paths
You already have a single-column helper pattern.
Add a companion helper for extracting several columns as a numeric matrix.

Suggested contract:
- verify all requested columns exist
- preserve requested column order
- return numeric matrix
- error clearly if any column is missing or not coercible to numeric

Use this helper anywhere `x_cols`, `D_cols`, or similar multi-column inputs are
resolved from `data`.

---

## Phase 2: Standardize validation style

### 2.1 Remove `stopifnot()` from user-facing paths
Instruction:
replace `stopifnot()` in:
- exported functions
- helpers directly reachable from exported functions
- simulation helpers with user-visible misuse modes

with explicit `if (...) stop(...)` checks.

### 2.2 Standardize validation order
In every exported function, adopt this order:

1. mutually exclusive / jointly required arguments
2. input container type
3. named-column existence
4. length / shape
5. NA / non-finite checks
6. semantic constraints
7. normalization / coercion
8. computation

### 2.3 Centralize repeated validators
Create or standardize a small helper family.

Recommended helpers:
- single-column extractor / validator
- multi-column extractor / validator
- vector-of-length-n validator
- matrix-of-shape validator
- binary 0/1 validator
- rank permutation validator
- jackknife feasibility guard
- scalar integer validator where relevant

### 2.4 Separate NA and non-finite checks everywhere
Make it consistent across all relevant validation routines.

Rule:
- `NA` gets an `NA`-specific message
- `Inf`, `-Inf`, `NaN` get a non-finite-specific message
- shape/type mismatch gets its own message

---

## Phase 3: Standardize error and warning behavior

### 3.1 Adopt one message style
Every user-facing diagnostic should say:

- what was expected;
- what was received;
- what to do instead, if helpful.

### 3.2 Use argument names in messages
Use the exact user-facing argument names in backticks.

### 3.3 Include actual values where useful
For scalar constraints:
- include `k`, `s`, `n`, etc.

For shape mismatches:
- include actual dimensions.

### 3.4 Make `call. = FALSE` the default for user-facing diagnostics
Do this consistently across exported code and nearby validation helpers.

### 3.5 Standardize warning policy for `return=`
Choose one repo-wide rule and apply it across all wrappers.

Recommended default:
- warn and ignore unknown `return` values while the package remains experimental.

Then ensure:
- `np_dnn()`
- `cate_dnn()`
- `plqr_dnn()`
all behave the same way.

---

## Phase 4: Standardize wrapper structure

Every high-level wrapper should have the same broad structure.

### Block 1: interface / option validation
- `method`
- `jackknife`
- `return`
- mutually exclusive arguments
- high-level guardrails such as `s < n` when delete-1 jackknife is requested

### Block 2: dispatch and input resolution
- `data.frame` path
- matrix path
- column extraction
- target resolution

### Block 3: canonicalization
- convert to `x_mat`, `x0_mat`, validated vectors
- no more container-type branching after this point

### Block 4: main computation
- weights
- scores / pseudo-outcomes
- estimate
- optional SE path

### Block 5: object assembly
- stable main fields
- optional fields added according to `return`
- class assignment

Instruction:
Audit `np_dnn()`, `cate_dnn()`, and `plqr_dnn()` against this template and
refactor only enough to make the structure visibly consistent.

---

## Phase 5: Standardize namespace and exported surface

### 5.1 Audit exports
For each exported function, ask:
- is this truly public?
- is it documented as public?
- is it tested as public?
- is it an extension point or a user-facing primitive?

### 5.2 Keep helpers internal by default
Validation helpers, dispatch helpers, and wrapper-only utilities should remain
internal unless there is a strong external use case.

### 5.3 Prefer dependency restraint
Before adding a new dependency, document why base R or an existing dependency is
not sufficient.

---

## Phase 6: Standardize generated-file discipline

### 6.1 Never hand-edit generated files
This includes:
- `NAMESPACE`
- `RcppExports.cpp`
- `RcppExports.R`

### 6.2 Regenerate in the same pass
If a source change requires a generated-file update, commit both together.

### 6.3 Add CI drift checks later
Add checks for stale generated files once the current standardization pass has
settled.

---

## Phase 7: Standardize compiled vs fallback behavior

### 7.1 Keep the validation boundary in R
Do not push user-facing ambiguity into compiled code if it can be resolved in R.

### 7.2 Add fallback-equivalence tests
For any compiled fast path with an R fallback, add tests that compare:
- output shape
- numerical agreement within tolerance
- behavior on edge cases

### 7.3 Document indexing and shape assumptions
For nontrivial compiled routines, make wrapper-boundary assumptions explicit.

---

## Phase 8: Standardize return objects

### 8.1 Stable fields
Each user-facing result object should have a documented stable core.

### 8.2 Optional fields
Optional fields should only be added when requested.
Do not let optional returns silently alter the meaning of the core estimate.

### 8.3 Class consistency
Ensure all user-facing estimator results have:
- predictable class names;
- print methods where appropriate;
- summary / predict / confint support where promised.

---

## Phase 9: Standardize comments and documentation

### 9.1 Comments
Keep comments focused on:
- why a constraint exists;
- how code maps to the method;
- why a default or fallback exists.

### 9.2 Roxygen updates
Audit exported functions for consistency in documenting:
- accepted input formats
- matrix vs `data.frame` interfaces
- jackknife restrictions
- `return=` behavior
- tibbles accepted through the `data.frame` path

### 9.3 Documentation examples
Make sure examples and vignettes follow the preferred current API and remain
deterministic.

---

## Phase 10: Standardize tests around validation and interfaces

### 10.1 Add failure-path tests for every new validation branch
At minimum add tests for:
- missing column names in `data=` workflows
- unknown `return` warnings
- `s >= n` when delete-1 jackknife is requested
- malformed custom distance output (wrong type, wrong length, non-finite)
- malformed `rank_weight$eval_batch()` output
- invalid `ranks`
- invalid `k`, invalid `s`, and `k > s`

### 10.2 Add interface-equivalence tests
For each high-level wrapper that supports both interfaces:
- construct the same problem through matrix inputs and through `data=...`
- assert equivalent estimates
- assert equivalent optional outputs where appropriate

### 10.3 Keep deterministic helper generators centralized
If test helpers such as random matrix builders are useful across multiple files,
keep them in shared helper files and set seeds explicitly inside them.

---

## Phase 11: Standardize lifecycle and change management

### 11.1 Treat public API changes explicitly
Changes to exported names, argument meanings, defaults, result fields, or
method behavior count as public API changes.

### 11.2 Stage deprecations when possible
Prefer warning-based transitions before hard removal unless the old behavior is
actively harmful.

---

## Phase 12: Suggested implementation order

Recommended order for the actual cleanup:

1. Adopt the numbered design-note corpus and `AGENTS.md` guidance.
2. Standardize wrapper structure across `np_dnn()`, `cate_dnn()`, `plqr_dnn()`.
3. Add shared column/matrix extraction helpers.
4. Remove remaining `stopifnot()` from user-facing paths.
5. Standardize `return=` warning behavior.
6. Standardize error-message wording across the main public API.
7. Audit exports and generated-file discipline.
8. Add validation-path tests and interface-equivalence tests.
9. Add compiled/fallback equivalence tests where missing.
10. Do a final roxygen / vignette / pkgdown pass.

This ordering minimizes churn and keeps the package usable throughout the
cleanup.

---

## Final rule

The cleanup is done when the following statement becomes true across the package:

> A user may call the high-level estimators with either explicit matrices/vectors
> or `data.frame`/tibble inputs, receive early and actionable diagnostics on bad
> inputs, and rely on the fact that all deeper computation happens on canonical
> numeric matrices and vectors under a stable validation contract, with tests,
> generated files, compiled paths, and docs all aligned to the same package-wide
> design rules.
