# Lifecycle and Deprecation Policy

This note standardizes how `RdistNN` should evolve while it remains an
experimental research package.

## 1. General rule

The package may evolve, but it should not evolve chaotically.

Even in an experimental package, changes to public behavior should be deliberate,
documented, and staged where possible.

## 2. Lifecycle posture

`RdistNN` retains its package-wide experimental lifecycle. Supported workflows
include CATE (`cate_dnn()`), nonparametric regression (`np_dnn()`, `np_tdnn()`),
conditional LATE (`clate_dnn()`), ordinary supplied-signal ratios
(`affine_ratio_dnn()`), and their supporting weights, tuning, diagnostics,
reporting, plots and frozen-model reuse.

Core CLATE support covers supplied-prediction binary-IV signals and the common
paired-ratio calculation under the documented pointwise inference conditions.
The signal constructor `iv_augmented_signals()` is part of that workflow.
Supporting computational options does not expand the theorem's statistical
scope. The optional `ratio_confidence_set()` inversion API remains experimental
and does not establish weak-instrument-uniform coverage. Public arguments,
result schemas and numerical statuses follow the deprecation discipline below.

All PLQR functionality is explicitly experimental: estimator and score
constructors, demo helpers, dedicated post-estimation methods, specialized
plots, and PLQR dispatch through shared generics. Local IV-GMM is also
experimental. These workflows remain available but fall outside the
supported workflow focus. Mark their status in help, articles, and
the pkgdown reference; do not change computation or emit new runtime warnings
solely to communicate lifecycle status.

This means:

- APIs may still improve and consolidate;
- internal structure may continue to change;
- some public interfaces may not yet be fully settled.

It does not mean:

- arbitrary breakage is acceptable;
- undocumented changes are acceptable;
- user-facing churn should be casual.

## 3. Public API changes

Any change to the following counts as a public API change:

- exported function name;
- argument name;
- argument meaning;
- default value;
- return-object fields;
- object classes;
- method behavior;
- warning or error behavior when it is part of the intended UX.

Public API changes should be treated as design decisions, not just code edits.

## 4. Preferred order of change

When a public interface needs improvement, prefer this sequence:

1. add the improved behavior;
2. preserve old behavior temporarily when feasible;
3. warn clearly on deprecated usage;
4. update docs and examples;
5. remove compatibility only after the change has been visible and tested.

## 5. Deprecation warnings

When deprecating behavior, warnings should say:

- what is deprecated;
- what to use instead;
- when or under what condition the old behavior is expected to disappear.

Do not emit vague deprecation warnings.

## 6. Soft vs hard deprecation

### Soft deprecation
Use when:
- the old interface still works;
- the package is transitioning toward a better API;
- compatibility is still valuable.

Behavior:
- keep functionality;
- warn clearly;
- update examples and docs toward the new path.

### Hard deprecation
Use when:
- the old behavior is actively harmful or misleading;
- maintaining compatibility creates serious inconsistency;
- there is no reasonable migration ambiguity.

Behavior:
- error clearly;
- explain what replaces the old behavior.

## 7. Return-object stability

Result objects should evolve conservatively.

Rules:
- do not remove core fields lightly;
- do not silently rename fields;
- do not change class structure without documenting the change;
- keep methods aligned with any object-structure updates.

Result-object changes count as API changes.

## 8. Experimental does not mean undocumented

If behavior changes in a user-visible way, update in the same pass:

- roxygen docs;
- examples if affected;
- vignettes if affected;
- tests if affected.

The package should never rely on users inferring the new behavior from code alone.

## 9. Backward-compatibility shims

Short-lived compatibility shims are acceptable when they reduce user pain.
They should not become permanent clutter.

If a shim is added, document:
- why it exists;
- what future cleanup it enables;
- what the intended end state is.

## 10. Review checklist

Before changing a public interface, ask:

- Is this really an improvement in consistency or usability?
- Does it break existing examples or return contracts?
- Should this be staged through a warning first?
- Have docs and tests been updated together with the code?
- Does this reduce long-run interface confusion?
