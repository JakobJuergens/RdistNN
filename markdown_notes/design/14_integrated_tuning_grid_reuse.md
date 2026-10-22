# Integrated tuning-grid reuse for DNN wrappers

This contract describes the implemented sensitivity workflow. Grid evaluation
is a mode of the existing estimators, not a separate estimator API. The shared
[result architecture](13_estimator_and_result_architecture.md) owns generic
reporting and plotting; the [ratio contract](17_clate_ratio_inference.md)
adds paired component covariance and target-specific numerical statuses.

## Public behavior and default anchor

`np_dnn()`, `cate_dnn()`, `plqr_dnn()`, `affine_ratio_dnn()` and `clate_dnn()`
accept `tuning_mode = c("default", "grid")`. A supplied tuning specification
resolves the default order and candidate grid. Without explicit `s` or tuning,
the estimator uses its documented paper-default rule. Explicit `s` stays
scalar; combining it with grid mode is an error.

Primary estimates, standard errors, metadata and requested extras always refer
to the default candidate. Grid mode adds `fit$tuning_grid` without changing
that primary result. It is sensitivity analysis, not automatic selection or
post-selection inference. The paper-rule parameter and candidate column are
named `lambda` / `lambda_grid`.

## Stored schemas

`tuning_grid` contains four fields: `candidates`, `estimates`, `diagnostics`
and `meta`.

| Candidate field | Meaning |
| --- | --- |
| `candidate_id` | Unique stable row identifier within the fit |
| `computed_key` | Computation identity; duplicate orders share a result |
| `is_default` | Exactly one anchor row |
| `feasible` | Whether the candidate satisfies the requested calculation |
| `status`, `reason` | Evaluation status and explanatory reason |
| `lambda`, `s` | DNN tuning label and resolved integer order |

IDs are deterministic and disambiguated when labels collide. Consumers should
select stored IDs rather than recreate their formatting. For repeated resolved
orders, plot labels include candidate IDs where needed to distinguish them.
Reconstruction validates retained controls, candidate/target identities and
cached permutations before computing a view; refit incomplete or inconsistent
state. If a custom grid
omits the resolved default, the resolver inserts an anchor row. Equal DNN
orders have equal computation keys, while distinct user rows remain present.

The long estimates table has candidate identity/labels, `target_id`,
`component`, `estimate`, `se`, `status` and `reason`. Scalar NP/CATE/ratio
components use `NA`; PLQR uses moment-component names. Missing SEs remain
`NA`. Confidence levels belong to reporting consumers, so fit tables do not
store confidence bounds. Ratio grids additionally retain `ratio_status`,
`covariance_status` and `wald_available`; candidate success is distinct from
an undefined ratio at an individual target.

The diagnostics table contains candidate-level weight information: `ess`,
`max_weight`, `top_n`, `top5_weight_share` and `degenerate`.
The concentration summary uses the first five rank positions,
independently of estimator `k`. Target-dependent support or score diagnostics
belong to estimator-specific diagnostic adapters. Metadata records the method,
order resolution, sample dimensions and candidate computation counts.

Do not store one model or weight matrix per candidate. Retain frozen model
state once and materialize candidate weights transiently or on request.

## Feasibility, failure and warnings

Validate all orders before heavy computation. Every candidate needs `k <= s`
and `s <= n`; active delete-one covariance additionally needs `s < n`.
An invalid default is an error. Infeasible non-default candidates retain rows
with `NA` estimates/SEs, status and reason, with a summary warning identifying
the skipped candidates. A failed non-default calculation is recorded without
aborting other candidates; a failed default propagates its error.

Duplicate computation keys reuse calculations but preserve user-facing rows.
The scalar near-uniform-weight warning applies to the default candidate;
non-default concentration is reported in the diagnostic table rather than
emitting repeated warnings during a sensitivity sweep.

## Shared engine and frozen state

`R/user-dnn-tuning.R` resolves orders and candidate identities.
`R/user-dnn-grid.R` provides `.dnn_grid_eval()`, driven by a candidate table,
cached ranks, a rank-weight provider and estimate/SE callbacks. Precomputed
default results anchor the grid. The engine deduplicates by `computed_key`,
evaluates feasible candidates and maps results back to every requested row.
Estimator callbacks own score construction, solving and uncertainty; the
engine does not assume every estimator is a scalar weighted mean.
NP, CATE and PLQR callbacks are owned by their respective `user-dnn-*-grid.R`
modules; paired-ratio callbacks remain in `user-dnn-affine-ratio-grid.R`.
The shared engine module owns infeasible-candidate warning text and the
NP/CATE/PLQR descriptive metadata builder. Estimator adapters retain their
existing schemas, private entry points and statistical policies.

Grid point estimates use the same observations, targets, nuisance predictions
and cached full-sample geometry. Changing `s` changes their rank weights;
it does not refit nuisance learners. The fast jackknife reuses this geometry
and rank tables, so its candidate loop requires no new ranking or RNG work.
Randomized preprocessing belongs before that loop.

The NP/CATE/PLQR `jackknife = "loo"` reference path remains an
explicit-deletion calculation: it recomputes weights and reranks retained
samples for each candidate. Grid acceleration therefore primarily concerns
`jackknife = "fast"`. With custom data-dependent or stateful metrics, the
reference loop can differ from a frozen-full-rank deletion and can consume RNG
if the metric itself does. Do not silently substitute the fast path for a
requested loop. Paired affine-ratio/CLATE estimators deliberately freeze the
full-sample geometry under both `"fast"` and `"loo"`, as specified in the
[ratio contract](17_clate_ratio_inference.md).

For table-capable rank-weight specifications, `spec$table(n)` returns a
length-n rank-to-weight vector, normalized consistently with the specification.
`.rank_weight_table()` validates this capability. Batch evaluation indexes the
table with each target's rank permutation. Arbitrary user-supplied rank rules
can lack this representation; preserve their supported evaluation fallback.

Candidate concentration summaries use actual averaged weights at every target.
The one-row-per-candidate display takes the minimum ESS and maximum absolute
weight/top-five absolute-mass share across targets; the near-uniform flag is
true if any target reaches the existing threshold. With duplicates, weight
profiles can differ by target even at one fixed candidate order.

## Rank weights and delete-one identity

For exact DNN with fixed `(n,k,s)`, the weight at rank `r` is

```
w_rank[r] = (1/k) * sum_{j=1}^k C(r-1,j-1) * C(n-r,s-j) / C(n,s)
```

Out-of-range binomials contribute zero; the unaveraged table's support ends at
rank `n-s+k`. X-aware construction averages each exact duplicate-covariate
block, so actual nonzero observation support can extend to the end of a block
crossing that edge. Reuse the shared block and deletion helpers for every candidate.
The table depends on ranks and parameters rather than observed signal values.
Asymptotic DNN also has a rank table. The low-level TDNN weight layer has
signed tables; its availability there does not authorize TDNN for CATE/CLATE.

Every retained rank vector after one deletion is a permutation of `1..n-1`.
Consequently, normalization by the sum of a rank table is the same across
all deletions and targets and can be performed once for that sample size.

Without duplicate-covariate blocks, full rank `r_i` and the size-`n-1` table
`w` give a frozen signal's deleted average

```
L_i = sum_{r < r_i} signal_(r) * w[r]
    + sum_{r > r_i} signal_(r) * w[r-1]
```

Prefix and shifted suffix sums calculate all deletions in a linear pass once
the rank order is known, per target and signal component. The compiled
implementation is `src/jackknife_rank_full.cpp`; the R reference and explicit
loop remain independent comparison paths. Never use unguarded zero indices
inside `ifelse()` for boundary ranks: R drops index zero and can silently
misalign a vector. Pad cumulative arrays or guard indices themselves.

With duplicate blocks, reuse the shared original-signal block deletion path:
the removed observation's block gets its actual survivor average from the
`n-1` table, and outside blocks use unshifted or shifted block averages.
Retain prefix/suffix and within-block correction states through their merge;
do not substitute rounded block-mean signals. The
[weighted accumulation reference](19_weighted_accumulation.md) owns this
arithmetic and its linear-work versus interleaved-group fallback bounds.

## Candidate views and plots

`candidate_fit(fit, id_or_index)` and `candidate_fits(fit, candidates)` return
scalar-style views with the parent's scientific class. They update the chosen
order, estimate, uncertainty, metadata and any requested `W`, `ess` or
`neighbors`. Neighbors are the first `k` observations, not the full DNN
support. Ratio views also rebuild components, covariance and optional
paired deletion/pseudo-value arrays. NP views update retained weight controls,
allowing joint target covariance
and prediction to use the selected order. Covariance is reconstructed on
request from common deletions. Retain candidate-independent nuisance and
geometry state; do not copy stale candidate-dependent quantities.

The common display grammar is `candidates = "default"`, `"all"` or candidate
IDs for grid-aware estimate/weight/diagnostic consumers. `plot_tuning()` uses
stored grid results directly; scalar-fit sensitivity uses retained model
state and the same internal grid engine, not repeated public estimator calls.
Stored grids reuse their stored uncertainty. Rowwise ratio Wald availability
must survive tidy extraction and plotting. Dashboards include tuning panels
when appropriate, under the same optional-package policy as scalar displays.

Scalar tuning plots re-evaluate all candidates under one
requested SE policy; stored grid fits reuse their recorded uncertainty.
An explicit scalar order's one-row bookkeeping payload is not a tuning sweep:
use the paper grid unless a tuning specification was retained or supplied.
Optional return fields cannot change these candidate orders.

## Simulation integration and verification

CATE simulation methods emit tidy rows per method and candidate. Benchmarks
that do not depend on `s` run once per replication and have missing candidate
labels/orders; do not duplicate them into the raw results for every candidate.
Summaries group by candidate identity with missing benchmark labels handled
explicitly. Reuse data, folds and nuisance predictions across candidates.
The fast-path candidate loop is RNG-neutral; custom metrics used by the reference
explicit-deletion path are subject to the caveat above.

The acceptance contracts are candidate/scalar equality, a preserved default
anchor, duplicate-order deduplication, explicit infeasible/failed rows, frozen
geometry and predictions, correct normalization at both sample sizes,
prefix/loop and R/compiled agreement, complete candidate views, serialization,
and consistent tidy-data/plot selection. Verify missing covariance and mixed
ratio statuses as well as ordinary finite output. The
[test-suite map](../../tests/README.md) identifies the corresponding suites.
