# RdistNN Design Notes

This bundle contains a proposed design-note system for the `RdistNN` package,
together with a concrete standardization plan for the existing code base.

## Suggested repo placement

Recommended location inside the repo:

- `markdown_notes/code_design_principles.md`
- `markdown_notes/input_validation_contract.md`
- `markdown_notes/api_shapes_and_input_formats.md`
- `markdown_notes/error_messages_and_warnings.md`
- `markdown_notes/comments_and_function_docs.md`
- `markdown_notes/testing_and_regression_policy.md`
- `markdown_notes/interface_standardization_decisions.md`
- `markdown_notes/existing_codebase_standardization_plan.md`

## Read order

1. `interface_standardization_decisions.md`
2. `code_design_principles.md`
3. `input_validation_contract.md`
4. `api_shapes_and_input_formats.md`
5. `error_messages_and_warnings.md`
6. `testing_and_regression_policy.md`
7. `existing_codebase_standardization_plan.md`
8. `comments_and_function_docs.md`

## Purpose of each note

- `code_design_principles.md`:
  repo-level design philosophy and layering rules
- `input_validation_contract.md`:
  package-wide validation rules and ordering
- `api_shapes_and_input_formats.md`:
  standard for public interfaces and canonical internal forms
- `error_messages_and_warnings.md`:
  style guide for diagnostics
- `comments_and_function_docs.md`:
  rules for comments and roxygen
- `testing_and_regression_policy.md`:
  test requirements for features, fixes, and validation behavior
- `interface_standardization_decisions.md`:
  concrete package-level decisions, especially on matrices/vectors vs data frames/tibbles
- `existing_codebase_standardization_plan.md`:
  actionable instructions for standardizing the current code base

## Core package-wide decision

User-facing estimators may accept either:

- explicit matrices and vectors; or
- `data.frame`-style inputs, including tibbles.

After validation and interface resolution, all internal computation should proceed
on canonical numeric matrices and vectors.
