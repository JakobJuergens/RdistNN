# Code Design Principles for RdistNN

This note records the core design principles that should guide package growth.
It is normative, not merely descriptive: new code should follow these rules
unless there is a documented reason to do otherwise.

## 1. Package identity

`RdistNN` is a modular toolkit for distance-localized estimation built around:

- distance specifications;
- rank-to-weight mappings;
- score / kernel specifications; and
- user-level estimators built on top of those abstractions.

The package is not just a collection of wrappers. Its long-run value comes from
maintaining a clean separation between:

1. reusable core primitives;
2. high-level estimation workflows; and
3. replication / simulation infrastructure.

New features should reinforce that structure rather than blur it.

## 2. Layering

The package should remain organized in three conceptual layers.

### Core layer
This includes:

- distances;
- order / rank helpers;
- rank-weight maps;
- score objects;
- jackknife primitives;
- shared validation and shape helpers.

This layer should be:

- small in conceptual surface area;
- reusable across estimators;
- strict about invariants; and
- as stateless as possible.

### User layer
This includes:

- `np_dnn()`
- `cate_dnn()`
- `plqr_dnn()`
- other future user-facing estimators

This layer should:

- validate early;
- accept friendly inputs;
- produce informative objects;
- fail with actionable messages; and
- hide internal machinery unless explicitly requested.

### Replication / simulation layer
This includes:

- DGP generators;
- paper simulation harnesses;
- result summarizers.

This layer may be more paper-specific, but it should still respect the same
standards for validation, determinism, and diagnostics.

## 3. Design for extension, not for cleverness

Prefer explicit interfaces over implicit magic.

Good extension points:

- `new_distance()`
- `new_rank_weight()`
- `new_score()`
- registries and `get_*()` accessors

Bad extension points:

- hidden assumptions about object shape;
- code paths that depend on undocumented fields;
- silent coercions that change semantics.

A new feature should have an obvious place in the existing abstraction tree.
If it does not, reconsider the design before writing more code.

## 4. Validate at the boundary

Validation should happen as close as possible to the user-facing call site.

The user should fail:

- before entering heavy computation;
- before entering LOO / jackknife loops;
- before compiled code is called;
- before fallback dispatch branches become complicated.

Do not postpone validation to deep internal helpers when the problem can be
detected at the outer API boundary.

## 5. Prefer robust fallbacks over brittle fast paths

Compiled and optimized paths are valuable, but correctness and usability are
more important than speed.

The package should prefer:

- correct and well-tested fallbacks;
- deterministic behavior;
- clear diagnostics when fast paths cannot be used.

Fast paths are an optimization, not the contract.

## 6. Determinism is part of correctness

Anything involving:

- rank order;
- ties;
- randomization;
- sampling;
- Monte Carlo infrastructure

must be deterministic when given the same seed and the same inputs.

Non-determinism is a bug unless it is an intentional and documented part of the
API.

## 7. Error handling is part of API design

Errors, warnings, and returned object shapes are part of the user interface.

A function is not complete when it only works on the happy path.
It is complete when:

- invalid inputs are rejected cleanly;
- ambiguous requests are clarified or blocked;
- optional behavior is documented through predictable return structure.

## 8. Keep internals small and named

When a validation or transformation pattern appears twice, extract it into a
small internal helper with a clear name.

Examples of good internal helpers:

- shape validators;
- column extractors;
- target resolvers;
- scalar guard checks.

Helpers should be:

- short;
- single-purpose;
- easy to test indirectly or directly.

## 9. Minimize hidden conventions

Do not rely on conventions that are only “known by the current author”.

Every important assumption should live in one of:

- the function signature;
- a helper with a descriptive name;
- roxygen docs;
- a design note like this one.

## 10. Backward compatibility posture

The package is still experimental, so APIs may evolve.
However, avoid gratuitous churn.

When changing public behavior:

- preserve object fields where reasonable;
- keep legacy aliases only if they add real value;
- update docs and vignettes in the same pass;
- add regression tests for the old failure mode or ambiguity.

## 11. Standard for merging nontrivial code

A nontrivial code change is not done until all of the following hold:

- the interface location is appropriate;
- validation is explicit;
- messages are user-facing and actionable;
- tests cover both happy path and failure path;
- docs or vignettes are updated if user-visible behavior changed.

## 12. Guiding question

Before adding code, ask:

> Is this making the package more modular, more predictable, and easier to use
> correctly?

If the answer is no, redesign first.
