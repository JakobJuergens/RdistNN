# Namespace and Dependency Policy

This note standardizes how `RdistNN` should manage its exported surface,
internal helpers, and package dependencies.

## 1. General rule

The package should expose a clean and deliberate public API.

Exports are part of the long-run maintenance burden. A function should be
exported only if at least one of the following is true:

- it is a primary user-facing entry point;
- it is a documented extension point;
- it is an S3 method that must be visible;
- it is a lower-level public primitive with a clear external use case.

Everything else should remain internal.

## 2. Export discipline

### Export when:
- the function has a stable user story;
- the function is documented in roxygen;
- the function has tests at the public API level;
- the function is expected to remain meaningful outside the current file layout.

### Do not export when:
- the function exists only to support one wrapper;
- the function is a validation helper;
- the function is transitional or convenience-only;
- the function is likely to be renamed or redesigned soon;
- the function is an implementation detail of a compiled or fallback path.

## 3. Extension points vs internal helpers

The package intentionally exposes extensibility through objects and registries,
such as distances, rank weights, and scores.

This is different from exposing internal plumbing.

Rule:
- extension points may be exported and documented;
- plumbing helpers should remain internal unless there is a strong external need.

## 4. Namespace stability

Avoid unnecessary export churn.

Adding an export is easy. Removing or renaming one later is costly.
Before exporting a new function, ask:

> Is this a real public interface, or just a currently useful helper?

If the answer is unclear, keep it internal.

## 5. Generated namespace files

`NAMESPACE` is generated and must not be edited by hand.

Rules:
- update roxygen comments in source files;
- regenerate `NAMESPACE` through the package tooling;
- do not make direct manual edits to generated namespace content.

Generated-file drift should be treated as a package maintenance issue.

## 6. Dependency discipline

Prefer a small and stable dependency footprint.

New Imports should be added only when they provide a clear benefit that outweighs:

- installation burden;
- CI burden;
- cross-platform complexity;
- maintenance cost.

Base R should remain the default unless a dependency clearly improves the package.

## 7. Imports vs Suggests

### Use Imports when:
- the package functionality truly depends on the dependency;
- the dependency is required for normal user-facing operation.

### Use Suggests when:
- the dependency is only needed for examples, tests, vignettes, benchmarking,
  optional workflows, or comparisons.

Suggested packages must be treated as optional everywhere.

## 8. Optional dependency behavior

If a feature uses a Suggested package, the code should:

- check for availability explicitly when needed;
- fail gracefully with a clear message;
- not make the entire package unusable.

Examples and vignettes should guard optional dependencies appropriately.

## 9. Documentation alignment

Every exported function should be documented and discoverable in one of:

- the pkgdown reference index;
- the vignettes;
- the README when it is a primary entry point.

An exported function with no clear documentation role should usually not be exported.

## 10. Review checklist

Before adding a new export or dependency, ask:

- Is this function truly public?
- Is it documented and tested as public?
- Does the dependency justify its long-run cost?
- Could this remain internal for now?
- Does this change keep the package architecture clearer, not noisier?
