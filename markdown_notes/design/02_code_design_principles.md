# Code design principles

RdistNN separates reusable primitives, user estimators and study infrastructure.
New code should fit those layers:

| Layer | Responsibilities |
| --- | --- |
| Core | Distances/ranks, rank weights, scores, jackknife primitives and shared validation. Keep explicit numeric contracts and little state. |
| User | Resolve friendly inputs, validate before expensive work, compute estimates and assemble informative classed results. |
| Simulation/replication | DGPs, orchestration and summaries with explicit seeds, schemas and study-specific assumptions. |

Use documented constructors and registries (`new_distance()`,
`new_rank_weight()`, `new_score()`, `get_*()`) as extension points. Internal
fields and silent coercion are not extension APIs. Give each function one
implementation owner and extract repeated validation/transformation into small,
single-purpose helpers. Important assumptions belong in signatures, named
helpers, documentation or the relevant design contract.

Validate at the public boundary before jackknife loops, native dispatch or
heavy computation. Keep fallbacks correct and testable; optimized code
implements the same selected method. Inputs, row order and seed must determine
ranks, ties, sampling and simulation results unless randomization is explicitly
part of the API.

Errors, warnings, output shapes and edge behavior are part of the interface.
Avoid gratuitous API churn even during the experimental lifecycle; use
[the lifecycle policy](07_lifecycle_and_deprecation_policy.md) for changes.
Nontrivial code needs appropriate placement, explicit validation, actionable
messages, deterministic valid/invalid/edge tests, and current user documentation.
