# Comments and Function Documentation

This note standardizes how code should be documented.

## 1. General rule

Comments should explain:

- why the code exists;
- what invariant it is enforcing;
- why a design choice was made;
- what subtlety would otherwise be easy to miss.

Comments should not restate obvious syntax.

## 2. Prefer structure over noise

Use comments to mark conceptual blocks, not every line.

Good pattern:

- short block header;
- 3–10 lines of focused code;
- next block header.

Bad pattern:

- line-by-line commentary of trivial assignments.

## 3. Comment the boundary between theory and code

This package is theory-linked.
Whenever code implements a mathematically meaningful construction, comments should
make the mapping legible.

Examples:

- what score is being computed;
- what jackknife path is being approximated;
- what weighting rule corresponds to the paper notation.

Do not assume future readers will remember the paper by heart.

## 4. Explain why a validation rule exists

When a constraint is not obvious from generic programming logic, explain it.

Example:
why delete-1 jackknife requires `s < n`.

Short rationale comments are valuable when they connect an implementation guardrail
to the underlying method.

## 5. Roxygen standard

Exported functions should have roxygen that covers:

- what the function does;
- the main input formats;
- the most important constraints;
- what is returned;
- key method choices when relevant.

Roxygen should prioritize usability over exhaustiveness.

## 6. Examples

Examples should be:

- small;
- deterministic;
- representative of the main usage;
- aligned with the current API.

Examples are part of the interface contract and should be maintained when APIs
change.

## 7. Internal helpers

Internal helpers do not always need full roxygen, but they should usually have:

- a one-line purpose comment;
- enough naming clarity that their contract is obvious.

If a helper has nontrivial behavior, add a short block comment above it.

## 8. Keep comments current

Outdated comments are worse than missing comments.

When behavior changes, update the nearby comment in the same commit.

## 9. Avoid speculative comments

Do not leave vague notes like:

- "maybe optimize later"
- "not sure if needed"
- "temporary fix"

If something is unresolved, put it in a design note or issue, not in production
code comments.

## 10. Style

Prefer comments that are:

- calm;
- direct;
- concrete;
- short enough to scan.

A good comment should make the next read faster, not slower.
