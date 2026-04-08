# Examples, Vignettes, and Pkgdown Policy

This note standardizes how user-facing documentation code should be written and
maintained across `RdistNN`.

## 1. General rule

Examples, vignettes, and the pkgdown site are part of the package interface.

They should reflect the preferred way to use the package now, not historical or
incidental usage patterns.

## 2. Examples

Examples in roxygen should be:

- small;
- deterministic;
- fast enough for routine checking;
- representative of the intended user workflow;
- aligned with the current API.

Examples should teach the preferred interface, not merely exercise code.

## 3. Determinism

All documentation code should be deterministic given a seed.

This includes:
- examples;
- vignette chunks;
- simulation illustrations;
- comparison workflows.

Do not rely on ambient RNG state.

## 4. Lightweight examples

Examples should avoid unnecessary runtime and dependency burden.

Prefer:
- small matrices;
- small demo datasets;
- short code paths;
- minimal optional dependencies.

If a richer workflow is valuable, put it in a vignette rather than a roxygen example.

## 5. Vignettes

Vignettes should have a clear role.

Recommended categories:
- getting started;
- conceptual / theory details;
- extension patterns;
- replication-specific workflows.

Each vignette should teach one coherent slice of the package rather than mix too
many purposes.

## 6. Preferred API in docs

Documentation should consistently privilege the preferred current interface.

If both matrix and `data.frame` interfaces are supported:
- explain both;
- present the preferred user path clearly;
- keep the examples aligned with the design notes.

## 7. Optional dependencies in docs

If examples or vignettes rely on Suggested packages:

- guard the code appropriately;
- fail gracefully when the package is unavailable;
- avoid making the entire document unusable.

Suggested-package workflows should be clearly optional.

## 8. Paper replication content

Replication-facing articles are valuable, but they should remain clearly labeled
as paper or study specific.

Rule:
- keep replication workflows separate from general user onboarding;
- do not let paper-specific material define the package’s primary public identity.

## 9. Pkgdown structure

The pkgdown site should mirror the conceptual architecture of the package.

Reference groups and article sections should remain aligned with:
- core primitives;
- high-level wrappers;
- extension interfaces;
- replication tools.

If the architecture changes, the site structure should change with it.

## 10. Documentation maintenance rule

When user-visible behavior changes, update in the same pass:
- examples;
- affected vignette text or code;
- pkgdown grouping if the conceptual surface changed.

Outdated docs are treated as a bug.

## 11. Review checklist

Before merging a documentation-affecting change, ask:

- Is the example deterministic and lightweight?
- Does it reflect the preferred current API?
- Is the vignette still conceptually focused?
- Are optional dependencies handled gracefully?
- Does the pkgdown organization still match the package architecture?
