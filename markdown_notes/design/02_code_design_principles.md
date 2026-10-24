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

## Captured configuration and retained work

Evaluate owned configuration arguments before a factory returns closures.
Validation or insertion into a result may already evaluate an argument;
otherwise use `force()` explicitly. Rebinding a caller's variable before first
use must not change the created metric, score or model. Validate properties
that need training data at setup. Arbitrary external state in a user callback
remains the callback author's responsibility. See
[Advanced R, function factories](https://adv-r.hadley.nz/function-factories.html).

Retain only context used by later computation. Scaled L1/L2 contexts hold
per-column inverse scales; evaluators operate on the supplied training matrix
and do not retain an unused scaled copy. Remove orphaned helpers when their
last consumer disappears. Use profiling and realistic allocation measurements
before making performance claims; see
[Advanced R, measuring performance](https://adv-r.hadley.nz/perf-measure.html)
and [doing less work](https://adv-r.hadley.nz/perf-improve.html).

## Maintaining this design reference

Record adopted design principles and best practices in the existing note that
owns the subject, with an authoritative source, the applicable boundary and
the regression contract when relevant. Distinguish a source's recommendation
from the package's chosen application of it. Keep notes synchronized when an
implementation changes; an intended policy alone does not prove completion.
Put session plans, execution logs and historical receipts in private working
records rather than turning these reference notes into an audit diary.
