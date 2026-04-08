# Agent Instructions for RdistNN

Read this file before modifying code. It is the canonical instruction file for
LLM agents and automated coding assistants working in this repository.

## Mandatory read order

Before making changes, read these files in order:

1. `markdown_notes/design/README.md`
2. `markdown_notes/design/interface_standardization_decisions.md`
3. `markdown_notes/design/input_validation_contract.md`
4. `markdown_notes/design/error_messages_and_warnings.md`
5. `markdown_notes/existing_codebase_standardization_plan.md`

For changes that affect comments, docs, or tests, also read:

- `markdown_notes/design/comments_and_function_docs.md`
- `markdown_notes/design/testing_and_regression_policy.md`
- `markdown_notes/design/api_shapes_and_input_formats.md`
- `markdown_notes/design/code_design_principles.md`

## Repository-wide non-negotiable rules

- User-facing high-level estimators may accept either:
  - explicit matrices and vectors; or
  - `data.frame` inputs, including tibbles.
- After validation and interface resolution, internal computation must proceed on
  canonical numeric matrices and vectors.
- Do not add `stopifnot()` to exported functions or to helpers on user-facing
  paths.
- Validation must happen early and produce actionable, argument-specific
  messages.
- Diagnostics should use clean user-facing `stop()` / `warning()` messages,
  generally with `call. = FALSE`.
- Unknown `return` values must be handled consistently across wrappers.
- New user-visible validation or behavior changes must include tests.
- Do not introduce new interface conventions without updating the design notes.

## Working rules for agents

When editing existing code:

1. Prefer standardizing toward the design notes over local convenience.
2. Preserve the public API unless there is a documented reason to change it.
3. Keep the matrix/data.frame dual-interface policy at the outer boundary only.
4. Reuse or extend shared validators instead of duplicating validation logic.
5. Keep error messages specific about what was expected, what was received, and
   what the user should do instead when helpful.
6. If a code path changes user-visible behavior, update tests and docs in the
   same pass.

## In case of conflict

If code, comments, and design notes conflict, treat the design notes as the
current source of truth unless there is a stronger repository-level instruction
added after this file.
