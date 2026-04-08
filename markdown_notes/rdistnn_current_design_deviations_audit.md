# RdistNN audit note: current deviations from the stated design principles

This note summarizes concrete places where the current `RdistNN` repo still
deviates from the numbered design-note corpus and the root `AGENTS.md` rules.

It is based on the current `master` state inspected on 2026-04-07/08.

---

## 1. Executive summary

The package is already moving strongly in the direction of the new design-note
system, but a few mismatches remain:

1. **`stopifnot()` still appears on exported or user-facing execution paths.**
2. **The canonical internal “numeric matrix/vector” boundary is not fully enforced.**
3. **Column extraction is not yet standardized across wrappers.**
4. **Some validation helpers still use “bad state” errors for ordinary user input failures.**
5. **The repo still contains a stale legacy design-note layer alongside the numbered corpus.**

The highest-priority fixes are:

- merge the cleanup PR that removes the old unnumbered design-note duplicates;
- replace remaining `stopifnot()` calls on exported/user-facing paths;
- enforce numeric canonical matrices at wrapper boundaries;
- centralize multi-column extraction helpers and route all `data=` interfaces through them.

---

## 2. Relevant governing principles

The following current repo principles are the main standards used for this audit.

### Root-level agent rules
From `AGENTS.md`:

- user-facing high-level estimators may accept matrices/vectors or `data.frame` inputs, including tibbles;
- after validation and interface resolution, internal computation must proceed on canonical numeric matrices and vectors;
- do not add `stopifnot()` to exported functions or to helpers on user-facing paths;
- validation must happen early and produce actionable, argument-specific messages;
- do not introduce new interface conventions without updating the design notes.

### Design-note rules
From the numbered corpus:

- `03_input_validation_contract.md`: no `stopifnot()` in exported functions or helpers on user-facing paths; validate early and explicitly.
- `04_api_shapes_and_input_formats.md`: public API may be flexible, but internal computation should move to canonical matrix/vector forms.
- `05_error_messages_and_warnings.md`: user-facing diagnostics should be argument-specific and not be framed as opaque internal bad states.
- `13_existing_codebase_standardization_plan.md`: standardize wrapper structure, matrix extraction, return behavior, and the cleanup toward one consistent note system.

---

## 3. Concrete deviations

## 3.1 Remaining `stopifnot()` on exported or user-facing paths

This is the clearest direct deviation from the stated rules.

### A. Exported demo/helper functions still using `stopifnot()`

#### `R/demo-data-npr.R`
Functions:
- `np_demo_data()`
- `np_demo_subsets()`
- `np_y_demo_sets()`

Current issue:
- these functions use `stopifnot(...)` for input validation.

Why this conflicts with the design notes:
- `03_input_validation_contract.md` explicitly says not to use `stopifnot()` in exported functions or in helpers on user-facing paths.
- `AGENTS.md` repeats this as a non-negotiable rule.

Recommended fix:
- replace `stopifnot(...)` with explicit validation using the standardized helper layer or explicit `if (...) .abort_bad_arg(...)` / `stop(..., call. = FALSE)` checks.

#### `R/demo-data-plqr.R`
Functions:
- `.r_ald()`
- `plqr_demo_data()`
- `plqr_moment()`
- `plqr_demo_subsets()`
- `plqr_demo_sets()`

Current issue:
- these functions still use `stopifnot(...)`.

Why this conflicts:
- same reason as above.
- `plqr_demo_data()`, `plqr_moment()`, and `plqr_demo_subsets()` are exported, so this is especially clear-cut.

Recommended fix:
- same replacement pattern as in the validation cleanup already applied elsewhere in the package.

#### `R/sim-cate-data.R`
Function:
- `sim_cate_data()`

Current issue:
- still starts with `stopifnot(kx >= 2L, clip > 0, clip < 0.5)`.

Why this conflicts:
- exported function;
- exactly the kind of user-visible API boundary that should now give explicit, argument-specific messages.

Recommended fix:
- replace with explicit scalar validation helpers and messages, e.g.:
  - `kx` must be an integer and `>= 2`;
  - `clip` must lie in `(0, 0.5)`.

### B. Exported function with user-facing fallback still using `stopifnot()`

#### `R/jackknife-rank.R`
Functions:
- `jackknife_fast_rank()`
- `.jackknife_fast_rank_R()`

Current issue:
- `.jackknife_fast_rank_R()` begins with `stopifnot(is.matrix(ranks))`.

Why this conflicts:
- `jackknife_fast_rank()` is exported;
- the R fallback is part of the user-facing execution path;
- your note system says helpers on user-facing paths should not use `stopifnot()` either.

Recommended fix:
- replace with explicit matrix validation and argument-aware error messages for `ranks`.

---

## 3.2 Canonical internal numeric matrices/vectors are not fully enforced

The design notes clearly say that once the flexible public interface has been
resolved, internal computation should proceed on canonical **numeric** matrices
and vectors.

That is not fully true yet.

### A. `np_dnn()` does not fully enforce numeric canonical matrices

#### `R/user-dnn-npr.R`
Current behavior:
- `x_mat` is created from either `data[, x_cols, drop = FALSE]` or `X`;
- it is passed through `.validate_matrix_nonempty()`.

#### `R/utils-validate.R`
Current issue:
- `.validate_matrix_nonempty()` delegates to `.validate_matrix(..., numeric = FALSE)`,
  so it does **not** enforce numeric type.

Why this conflicts:
- `AGENTS.md` and `01_interface_standardization_decisions.md` say the internal
  boundary should be canonical numeric matrices and vectors.
- `04_api_shapes_and_input_formats.md` says internal computation should move to
  canonical matrix/vector forms after validation.

Consequence:
- non-numeric or unexpectedly coercible tabular inputs may enter deeper weight
  and distance code paths before failing.

Recommended fix:
- either:
  1. make `.validate_matrix_nonempty()` numeric by default for covariate/design inputs, or
  2. introduce a dedicated `.validate_numeric_matrix_nonempty()` / multi-column extractor helper
     and use it systematically in wrappers.

### B. `plqr_dnn()` has the same issue for `X`, `D`, and `mD`

#### `R/user-dnn-plqr.R`
Current behavior:
- `x_mat`, `d_mat`, and `mD_mat` all go through `.validate_matrix_nonempty()`,
  which currently does not enforce numeric type.

Why this conflicts:
- same internal-boundary rule as above.

Recommended fix:
- route all three through a numeric matrix validator.

### C. `cate_dnn()` only partially conforms

#### `R/user-dnn-cate.R`
Current behavior:
- many scalar/vector nuisance columns already go through `.dnn_col_vec()`, which is good;
- but `x_mat` still uses `data[, x_cols, drop = FALSE]` plus `.validate_matrix_nonempty()`.

Why this conflicts:
- the covariate matrix is still not guaranteed to be canonical numeric at the boundary.

Recommended fix:
- make `x_mat` extraction flow through the same stricter numeric-matrix helper that should also be used in `np_dnn()` and `plqr_dnn()`.

---

## 3.3 Column extraction is not yet standardized across wrappers

Your design notes call for shared extraction helpers, especially for column-based
`data=` interfaces. Right now this is only partially implemented.

### A. `cate_dnn()` is ahead of the others
#### `R/user-dnn-cate.R`
Good current pattern:
- `.dnn_col_vec()` is used for `y_col`, `t_col`, `e_col`, `m0_col`, `m1_col`, `m_col`.

This is aligned with the design-note direction.

### B. `np_dnn()` still uses raw slicing for the data-frame path
#### `R/user-dnn-npr.R`
Current behavior:
- `data[, x_cols, drop = FALSE]`
- `data[[y_col]]`

Why this conflicts:
- missing-column and bad-type diagnostics are not fully centralized;
- behavior may diverge from the `cate_dnn()` style.

### C. `plqr_dnn()` still uses raw slicing heavily
#### `R/user-dnn-plqr.R`
Current behavior:
- `data[, x_cols, drop = FALSE]`
- `data[, D_cols, drop = FALSE]`
- `data[, mD_cols, drop = FALSE]`
- `data[[q_col]]`
- `data[[y_col]]`

Why this conflicts:
- the design notes call for a shared extraction layer;
- the code base still has multiple wrapper-specific extraction idioms.

Recommended fix:
- add a multi-column extractor/validator helper, for example:
  - one helper for numeric vector extraction from `data`;
  - one helper for numeric matrix extraction from `data` with ordered columns;
- then refactor all wrappers to use them.

---

## 3.4 `abort_bad_state()` is still used for ordinary user-input validation

The error-handling notes distinguish between:

- ordinary user argument failures; and
- true internal impossible states.

At the moment, some validation helpers still route ordinary input problems
through `.abort_bad_state()`.

#### `R/utils-validate.R`
Examples:
- `.validate_ncol_relation()`
- `.validate_same_n()`
- parts of `.validate_subset_indices()`
- parts of `.validate_x0_inv_s()`

Why this conflicts:
- `05_error_messages_and_warnings.md` says internal impossible states should be rare
  and distinct;
- many of these failures are really user-facing argument errors, not corrupted internal state.

Why this matters:
- messages become less consistent with the new argument-aware style;
- this makes the package less aligned with its own new diagnostic contract.

Recommended fix:
- for validation helpers used on user inputs, prefer `.abort_bad_arg()` and argument-specific messages;
- reserve `.abort_bad_state()` for genuinely impossible internal states.

---

## 3.5 The repo still contains a stale legacy design-note layer

This is a repo-structure mismatch rather than a runtime-code mismatch, but it is real.

### Current state
The root `AGENTS.md` says the numbered note corpus is the source of truth.
However, the repo still contains:

- old unnumbered design-note files under `markdown_notes/design/`
- the old unnumbered `markdown_notes/existing_codebase_standardization_plan.md`
- a design-note `README.md` that still reflects the older layout

This means the repo currently exposes two overlapping note systems.

Why this conflicts:
- the new note system was meant to create a single ordered canonical corpus;
- the presence of both versions still creates ambiguity for humans and LLM agents.

Recommended fix:
- merge the cleanup PR that:
  - updates `markdown_notes/design/README.md`;
  - removes the unnumbered duplicate files;
  - removes the old unnumbered standardization-plan copy.

---

## 4. Things that already align well with the new design principles

This audit is about deviations, but it is worth recording where the repo is
already conforming well.

### Validation helper layer
`R/utils-validate.R` already provides a substantial shared validation layer:
- `.abort_bad_arg()`
- `.validate_scalar_int()`
- `.validate_scalar_prob()`
- `.validate_binary01()`
- permutation/subset validators
- `.warn_unknown_return()`

This is exactly the right direction.

### Wrapper-level unknown `return` policy
The high-level wrappers now warn on unknown `return` values via `.warn_unknown_return()`,
which is aligned with the current experimental-package policy.

### High-level wrapper structure
`np_dnn()`, `cate_dnn()`, and `plqr_dnn()` are already much closer than before to
the intended common wrapper structure:
- option validation;
- dispatch;
- target resolution;
- weight computation;
- estimate;
- optional jackknife;
- output assembly.

### Validation-path tests
At least for `np_dnn()` and `plqr_dnn()`, there are now meaningful validation
tests covering:
- bad `k` and `s`;
- bad `tau`;
- `s < n` jackknife requirement;
- bad interface combinations;
- unknown `return` warnings.

So the code base is not far from the design-note standard; the remaining issues
are mostly cleanup and consistency work.

---

## 5. Priority remediation plan

## Highest priority

### 1. Merge the design-note cleanup PR
Reason:
- removes ambiguity immediately;
- makes the note corpus truly canonical.

### 2. Remove remaining `stopifnot()` from exported/user-facing paths
Target files:
- `R/demo-data-npr.R`
- `R/demo-data-plqr.R`
- `R/sim-cate-data.R`
- `R/jackknife-rank.R`

Reason:
- direct conflict with stated rules;
- cleanest visible design-note mismatch.

## Second priority

### 3. Enforce numeric canonical matrices at wrapper boundaries
Target files:
- `R/user-dnn-npr.R`
- `R/user-dnn-plqr.R`
- `R/user-dnn-cate.R`
- `R/utils-validate.R`

Reason:
- this is one of the core architectural promises of the interface notes.

### 4. Add shared multi-column extraction helpers
Target files:
- `R/utils-validate.R`
- all three user wrappers

Reason:
- closes the remaining gap between `cate_dnn()` and the other wrappers;
- standardizes missing-column and bad-type diagnostics.

## Third priority

### 5. Reduce `.abort_bad_state()` use in ordinary validation helpers
Target files:
- `R/utils-validate.R`

Reason:
- improves consistency with the newer error-message design.

---

## 6. File-by-file TODO list

## `R/demo-data-npr.R`
- replace `stopifnot()` in exported/user-facing functions
- convert to argument-specific validation messages

## `R/demo-data-plqr.R`
- replace `stopifnot()` in exported/user-facing functions
- convert to argument-specific validation messages

## `R/sim-cate-data.R`
- replace exported-function `stopifnot()` with explicit validation helpers

## `R/jackknife-rank.R`
- remove `stopifnot()` from `.jackknife_fast_rank_R()`
- replace with explicit validation for `ranks`

## `R/utils-validate.R`
- add a dedicated numeric matrix extractor/validator for wrapper boundaries
- add a multi-column `data.frame` extractor
- convert ordinary user-input failures from `.abort_bad_state()` to `.abort_bad_arg()` where appropriate

## `R/user-dnn-npr.R`
- route `x_mat` and `y` extraction through stricter canonical helpers

## `R/user-dnn-plqr.R`
- route `x_mat`, `d_mat`, `mD_mat`, `q_vec`, `y_vec` through standardized helpers

## `R/user-dnn-cate.R`
- route `x_mat` through the same standardized numeric matrix helper used by the other wrappers

## repo structure
- retire old unnumbered design-note duplicates after cleanup PR merge

---

## 7. Bottom line

The repo is already broadly aligned with the new design-note system, but it has
not fully converged yet.

The most important current mismatches are:

- lingering `stopifnot()` on exported/user-facing paths;
- incomplete enforcement of canonical numeric matrices/vectors internally;
- wrapper-specific rather than standardized column extraction;
- stale legacy note duplication in the repo.

This is a **cleanup-and-standardization** problem, not a fundamental design problem.
The note system is already good; the remaining work is to make the code base and
repo layout fully obey it.
