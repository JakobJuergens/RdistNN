# Integrated tuning-grid reuse for DNN wrappers

*Created: 2026-07-05. Revised: 2026-07-05, after a code-level review of the
tuning, weight, jackknife, plotting, and simulation layers. This revision
resolves the previously open design questions, corrects the runtime analysis,
replaces the proposed rank-cache helper with reuse of existing internals, and
expands the note into a full implementation guide covering the estimators, the
simulation layer, result-object ergonomics, and the plotting family. It is the
authoritative specification for the grid-reuse work; deviations during
implementation must be recorded back here.*

---

## 1. Purpose

This note records the package-side design for evaluating a grid of DNN kernel
orders without rerunning the entire estimator workflow for each value of
`kappa` or `s`.

The immediate motivation is the corrected CATE simulation workflow for the
conditional Z-estimation paper. After correcting the paper-facing simulations
to use `k = 1`, the old grid of `kappa` values is less informative than it
looked under the mistaken `k = 5` setup. The manuscript simulation code now
needs a larger grid of `kappa` values to show the bias-variance and
localization tradeoff as `s` varies. Running that grid by launching a complete
simulation run per `kappa` is wasteful, because many expensive pieces do not
depend on `s`.

The same issue is not merely a simulation inconvenience. Empirical users will
naturally want sensitivity checks over plausible `kappa` values. The package
should make that workflow ergonomic and fast enough to use routinely.

---

## 2. Central Decisions

> **Decision 1 — integrated mode, not a parallel API.** Grid evaluation is a
> mode of the existing estimators and simulation functions. Do not create new
> exported functions whose names encode grid-specific signatures, such as
> `cate_dnn_grid()`, `np_dnn_kappa_grid()`, or `sim_cate_kappa_grid_run()`.
>
> **Decision 2 — algebraic groundwork first.** Before any grid work, make
> rank-to-weight tables first-class in the `rank_weight` contract and replace
> the O(n²) fast-jackknife loop with an O(n) prefix-sum evaluation
> (Section 4). This is a pure internal speedup of the scalar path, it is where
> most of the current runtime actually goes, and it makes grid-mode standard
> errors nearly free per candidate.
>
> **Decision 3 — engine over callbacks, not over cross-products.** The shared
> grid engine's contract is "candidate table + weight provider + estimator
> callback". It must not assume the estimate is `crossprod(signal, W)`,
> because `gmm_dnn()` (a stated future target) solves a nonlinear program per
> target.

Decision 1 is consistent with the existing design notes: user-facing
estimators keep a stable grammar, tuning lives in the tuning abstraction, and
optional behavior is additive rather than forcing users onto a parallel API
surface. Decisions 2 and 3 are new in this revision and are justified in
Sections 4 and 8.

---

## 3. Current State

The package already has most of the pieces needed. The earlier draft of this
note understated how much already exists; implementers should read this
section against the code before writing anything.

### 3.1 Tuning layer (`R/api-tuning.R`)

`tuning_dnn_paper()` stores both a default `kappa` (used to choose the primary
`s`) and a `kappa_grid` (used by `tuning_eval(..., mode = "grid")`).
`tuning_eval()` returns candidate metadata with `kappa`, `s`, `feasible`, and
`reason`, and its grid validator (`.tuning_validate_grid_output()`) already
computes per-candidate feasibility from `k <= s` and, when a jackknife is
requested, `s < n`.

### 3.2 Candidate resolution (`R/user-dnn-core.R`, `.dnn_resolve_s()`)

`.dnn_resolve_s()` already evaluates **both** modes on every tuning-based fit:
it calls `tuning_eval(mode = "default")` and `tuning_eval(mode = "grid")` and
stores both in the tuning payload (`payload$default`, `payload$grid`). Its
manual-`s` branch already constructs a one-row candidate grid with the same
columns. The "one schema for scalar and grid mode" goal is therefore already
met at the resolver level. What is missing is only:

- identifier columns on the candidate table (Section 6.5), and
- wrappers that *evaluate* the non-default candidates.

### 3.3 Estimator wrappers

`np_dnn()` (`R/user-dnn-npr.R`), `cate_dnn()` (`R/user-dnn-cate.R`),
`plqr_dnn()` (`R/user-dnn-plqr.R`), and eventually `gmm_dnn()`
(`R/user-dnn-gmm.R`) share the pattern:

1. validate the public interface;
2. resolve `data=` or matrix inputs into canonical matrices and vectors;
3. resolve a scalar `s` through `.dnn_resolve_s()`;
4. resolve targets via `.dnn_resolve_targets()`;
5. compute DNN weights via `.dnn_weights_with_extras()`;
6. combine weights with a signal or moment
   (`.dnn_estimate_from_signal()` for NP/CATE/PLQR; a GMM solve for GMM);
7. optionally compute delete-1 jackknife SEs via
   `.dnn_jackknife_se_loo_matrix()`;
8. attach return extras and hidden `.model` state.

This structure is good and is reused, not replaced. Two facts matter for the
grid design:

- **The wrappers already cache ranks.** With `jackknife = "fast"`, they
  request `"ranks"` from `weights_dnn()` once and pass them to the jackknife
  as `precomputed_ranks`. The grid path extends this existing caching; it does
  not introduce caching.
- **`gmm_dnn()` is not a signal-times-weights estimator.** It runs two-step
  GMM with `stats::optim()` per target; `W` enters the weighted moment
  `gbar`. Any engine abstraction that assumes a signal matrix will not extend
  to it.

### 3.4 Weight and rank layer (`R/api-weights.R`, `R/wrank-*.R`)

`weights_matrix()` is ranker-agnostic, and `.resolve_ranking()` is the
canonical distances → orders → ranks path with the deterministic tie-breaker
`.order_with_ties()` (distance, then index). For fixed `X`, `X0`, distance,
and tie rule, the rank matrix does not depend on `s`; only the rank-to-weight
map does.

What the current layer does *not* exploit is that the DNN rank-to-weight map
is a pure table lookup (Section 4.1):

- `rank_weight_dnn_exact()$eval_batch` and `rank_weight_dnn_asymp()$eval_batch`
  (`R/wrank-dnn-wrapper.R`) loop over target columns and call the
  `*_from_ranks()` builders per column, recomputing the identical table `m`
  times per evaluation;
- the fast jackknife calls `rank_weight_eval_batch()` once per deletion,
  recomputing the identical size-`(n-1)` table `n` times.

### 3.5 Fast jackknife (`R/user-dnn-core.R`, `.dnn_jackknife_se_loo_matrix()`)

The `"fast"` path already uses precomputed full-sample ranks and the
closed-form delete-1 rank update

```r
r_i            <- matrix(ranks_full[i, ], nrow = n - 1L, ncol = m, byrow = TRUE)
rank_matrix_i  <- ranks_full[keep, , drop = FALSE]
rank_matrix_i  <- rank_matrix_i - (rank_matrix_i > r_i)
```

That is the right direction, but it is scalar in `s`, loops over deletions in
R with per-iteration matrix allocations, and rebuilds rank-weight tables
inside the loop. **This loop — not nuisance fitting, not distance evaluation —
is the dominant per-fit cost at realistic simulation sizes** (n = 2000: an
O(n²m) R-level loop versus a handful of `glm`/`lm` fits in the low
milliseconds).

### 3.6 Simulation layer (`R/sim-cate-run.R`, `R/sim-cate-benchmarks.R`)

Registered methods implement
`.method_XXX(dgp, x_cols, k, s, alpha, ...) -> named numeric vector`, consumed
positionally through an `na_row` template in `sim_cate_one_rep()`.
`sim_cate_run()` resolves one scalar `s` per cell (from `s_fn` or
`tuning_eval(mode = "default")`) and prepends cell identity columns. The
external kappa-grid workaround in the paper repository loops over whole
`sim_cate_run()` calls.

### 3.7 Plotting (`R/plot-tuning.R`)

`plot_tuning()` answers exactly the grid question ("how sensitive is the
estimate to `s`?") and currently does so by **refitting the estimator at each
candidate through the public wrappers**, including a full delete-one jackknife
per refit. This work is concurrent with this note; Section 9 defines how the
two must converge.

---

## 4. Algebraic Foundations

Four facts about the DNN weight construction determine the efficient design.
They are recorded here in full because they change the runtime story, the
implementation plan, and the diagnostics design. Nontrivial conclusions are
recorded per the repository rule that derivations must land in code, tests, or
design notes.

### 4.1 Fact: DNN weights are pure rank-table lookups

For fixed `(n, k, s)`, the exact DNN weight of an observation depends only on
its rank. `dnn_weights_exact_full_from_ranks()` (`R/wrank-dnn_exact.R`)
computes a length-`n` vector `w_rank` from the hypergeometric law

```
w_rank[j] = (1/k) * sum_{i=1}^{k} C(j-1, i-1) * C(n-j, s-i) / C(n, s)
```

(with out-of-range binomials contributing zero; exact support `j <= n - s + k`)
and then maps `w = w_rank[ranks]`. **The table does not depend on the data.**
The same holds for the asymptotic DNN rule
(`dnn_weights_asymp_full_from_ranks()`) and for the TDNN tables
(`tdnn_rank_table_exact()`, `tdnn_rank_table_asymp()`; signed entries, columns
still sum to 1).

Consequences:

- a batch weight evaluation is `W <- matrix(tab[R], nrow(R), ncol(R))` — one
  O(n) table build per unique `(n, k, s)`, then pure indexing;
- in the delete-1 jackknife, every deletion has size `n - 1`, so **one** table
  serves all `n` deletions and all `m` targets;
- the `table(n)` builders already exist implicitly:
  `dnn_weights_exact_full_from_ranks(seq_len(n), k, s)$w_rank` *is* the table
  (identity permutation), so Phase 0 adds no new weight mathematics.

### 4.2 Fact: normalization constants are deletion-invariant

`rank_weight_dnn_asymp()` (and defensively also the exact spec) normalizes the
realized weights by their sum. Every deletion of one observation realizes each
rank `1..n-1` **exactly once** (ranks are permutations, and the delete-1
update is a bijection onto `1..n-1`). Hence the realized sum equals
`sum(tab)` for every deletion and every target, and the table can be
pre-normalized once. The same argument covers TDNN's `normalize = TRUE`. This
is what licenses folding normalization into the table on both the point and
jackknife paths.

### 4.3 Fact: the fast delete-1 jackknife collapses to prefix sums

Fix one target column with full-sample ranks `r_1, ..., r_n` (a permutation of
`1..n`) and let `w[.]` be the rank-weight table for sample size `n - 1`
(pre-normalized per 4.2). Under deletion of observation `i`, kept observation
`j` has rank `r_j - 1{r_j > r_i}`, so the leave-one-out estimate for a signal
vector `sig` is

```
est_i = sum_{j != i} sig_j * w[ r_j - 1{r_j > r_i} ].
```

Re-index the sum by rank, writing `sig_(r)` for the signal of the observation
whose full-sample rank is `r`:

```
est_i = sum_{r < r_i} sig_(r) * w[r]  +  sum_{r > r_i} sig_(r) * w[r - 1].
```

Both sums are cumulative sums over the rank-ordered signal. Define

```
A[t] = sum_{r <= t} sig_(r) * w[r]          (t = 1..n; w[n] treated as 0)
B[t] = sum_{r >= t} sig_(r) * w[r - 1]      (t = 1..n; w[0] treated as 0)
```

Then all `n` leave-one-out estimates are

```
est_i = A[r_i - 1] + B[r_i + 1]
```

with the convention `A[0] = 0` and `B[n + 1] = 0`. This is **O(n) per
(target, s)** after one rank sort, instead of the current O(n²) per-deletion
loop. For a q-column signal (PLQR moment components), run the two passes per
component: O(nq) per (target, s). The rank sort (`order(ranks[, j])`) is per
target and is reused across all `s` values, so over a grid of `U` unique
feasible `s` values the jackknife costs `m` sorts plus `O(n m q U)` cumulative
sums.

**Verification record.** The identity was checked numerically against the
current per-deletion loop algebra (the exact update from Section 3.5, with the
exact hypergeometric table) at `n = 400`, `k = 1`, `s ∈ {3, 17, 60}`, on
random signals and random rank permutations: maximum absolute deviation was at
double-precision noise (~1e-16 to ~1e-17), and an independent tiny case
(`n = 6`) was verified element by element. The implementation must add a
package-level equivalence test against the existing loop (kept as the
reference path), per the fallback-equivalence rule in `AGENTS.md`.

**Known implementation pitfall (hit during verification).** A vectorized
implementation of the boundary cases must not index with `r_i - 1L` directly
inside `ifelse()`: when some `r_i = 1`, the index `0` produces a *shorter*
vector (R drops zero indices) and silently misaligns every subsequent element.
Guard the index itself:

```r
loo <- ifelse(r > 1L, A[pmax(r - 1L, 1L)], 0) +
       ifelse(r < n,  B[pmin(r + 1L, n)],  0)
```

or use explicit padding (`A2 <- c(0, A); B2 <- c(B, 0)` and index shifted).
The equivalence test must include targets where some observation has rank 1
and rank n (always true for permutations, but the test should not mask the
boundary by accident).

**Applicability conditions**, all of which hold on the current fast path:

- the weight rule is a pure rank-table lookup for the given sample size
  (exact DNN, asymptotic DNN, TDNN);
- ranks are deterministic-tie-broken permutations (guaranteed by
  `.resolve_ranking()`);
- `s <= n - 1`, which is exactly the existing `s < n` jackknife feasibility
  condition;
- any normalization constant is deletion-invariant (Fact 4.2).

Custom user-supplied `rank_weight` specs are arbitrary functions of the rank
vector and need not be table-representable. For those, the existing
per-deletion loop remains the correct fallback (Section 8.1).

### 4.4 Fact: candidate diagnostics are table functionals

For a pure rank-table weight, each column of `W` is a permutation of the table
entries. Therefore, per candidate:

- `colSums(W^2)` is constant across targets and equals `sum(tab^2)`;
- effective sample size `ess = 1 / sum(tab^2)` is **one scalar per
  candidate**, identical at every target;
- `max(abs(W[, j]))` is constant across targets and equals `max(abs(tab))`;
- top-weight shares are functions of `sort(tab)` only.

Consequences: grid-mode weight diagnostics are O(n) per candidate, computed
from the table without materializing any candidate's `W`; and the degenerate-
weights check (`.dnn_warn_degenerate_weights()`'s `ess >= 0.9 n` criterion)
reduces to one scalar comparison per candidate. Data-dependent diagnostics
(for example CATE overlap measures that combine `W` with treatment status) are
the exception and may still need `W`; the engine therefore materializes each
candidate's `W` transiently (Section 8.3) but never stores it.

These four facts move most of the value of "grid reuse" into the scalar path
itself. After Phase 0, the marginal cost of an additional unique candidate `s`
is: one O(n) table for the point estimate, one O(n) table for the jackknife,
one O(nm) indexing pass, and one O(nmq) prefix-sum pass.

---

## 5. What Should Not Be Built

Do not add one-off exported functions such as:

- `cate_dnn_grid()`;
- `np_dnn_grid()`;
- `plqr_dnn_grid()`;
- `cate_dnn_kappa_grid()`;
- `sim_cate_kappa_grid_run()`;
- `weights_dnn_kappa_grid()`.

Those names solve the immediate problem by duplicating the public surface, and
they create long-run problems: every estimator would need a scalar and a grid
function; documentation and examples would split across two APIs; users would
have to learn which arguments belong to which call; simulation code would
drift away from the estimator API empirical users see; and future wrappers
(a general conditional Z-estimation front end) would inherit the naming debt.
The package already has a tuning abstraction. The grid is a tuning behavior,
not a separate estimator.

The one deliberate exception to "no new exported functions" is the pair of
result accessors `candidate_fit()` / `candidate_fits()` (Section 6.6): the ban
targets parallel *estimation* APIs, not accessors over an already-computed
fit.

In addition, the following **internal** machinery must not be built:

- **No `.dnn_compute_rank_cache()` helper.** The earlier draft proposed a new
  helper that would "use the same distance resolution, batch-distance
  evaluation, and deterministic tie-breaking as `weights_matrix()`". That is a
  specification for duplicating `.resolve_ranking()`. The rank cache the grid
  path needs is simply the `ranks` matrix (and optionally `ord`) that the
  wrappers already obtain from their single `weights_dnn()` call. Reuse that
  call; never add a second distance-touching code path.
- **No replacement of `.dnn_resolve_s()`.** It already returns the candidate
  grid with feasibility metadata in both the tuning and manual-`s` branches.
  Extend its candidate table with identifier columns (Section 6.5); do not
  introduce a second resolver.
- **No second grid engine inside `plot_tuning()`** (Section 9).
- **No per-candidate weight storage by default**, and no repurposing of
  `return = "W"` (Section 13).
- **No stored list of per-candidate fit objects.** A grid fit must not embed
  one `np_result`/`cate_result` per candidate: each copy would duplicate
  `.model` state (including the `n x p` training matrix), `x0`, and metadata,
  and it would make the primary fields ambiguous. Per-candidate fit objects
  exist, but as on-demand views via `candidate_fit()` (Section 6.6), never as
  stored payload.
- **No parallelism inside the candidate loop.** Rep-level parallelism owns the
  cores in simulations (the same reasoning as `num.threads = 1L` for the GRF
  benchmarks); the candidate loop is cheap after Phase 0 and must stay
  single-threaded and allocation-light.

---

## 6. Public API Specification

### 6.1 The mode argument

One shared control on the high-level DNN wrappers:

```r
tuning_mode = c("default", "grid")
```

This name is **decided** (previously open): it is explicit, short, parallel to
the existing `tuning` argument, and mirrors `tuning_eval(mode =)`. The mode
lives **outside** the tuning object, deliberately: the same tuning spec (for
example the paper spec with its `kappa_grid`) must be usable for a scalar fit
and for a sensitivity fit without rebuilding the spec, and simulations must be
able to flip the mode without touching the spec.

`tuning_mode = "default"` reproduces today's behavior exactly. All wrappers
adding the argument must place it directly after `tuning` in the signature and
resolve it with `match.arg()`.

### 6.2 Interaction rules

- `tuning_mode = "grid"`, `tuning` supplied: resolve candidates from the
  tuning object; evaluate all feasible candidates.
- `tuning_mode = "grid"`, no `tuning`: use `tuning_dnn_paper()` — and
  therefore its default `kappa_grid` — exactly as the scalar default path
  already falls back to the paper rule inside `.dnn_resolve_s()`. This is the
  intended ergonomic path for a quick sensitivity check and requires no new
  logic.
- `tuning_mode = "grid"` together with explicit scalar `s`: **error**, in the
  style of the existing `s`/`tuning` exclusivity error (`.abort_bad_arg()`
  with an actionable message: grid mode requires a tuning object; omit `s`).
- `tuning_mode = "default"` behaves exactly as today in every combination.

### 6.3 Explicit `s` stays scalar

The public wrappers keep scalar `s`. Users who want a grid use `tuning`.
A raw vector `s = c(15, 23, 31, 47)` lacks the `kappa` labels, formula
metadata, feasibility metadata, and default-candidate semantics of the tuning
API. If vector-valued `s` is ever added, it must be translated into the
internal candidate schema of Section 6.5 (as a manual grid with `kappa = NA`)
and must not produce a separate return shape. Deferred; see Section 17.

### 6.4 Return shape

Primary fields remain anchored at the **default candidate** in both modes:

- `fit$tau` (`cate_dnn`), `fit$mu` (`np_dnn`), `fit$moment` (`plqr_dnn`);
- `fit$se`;
- `fit$meta` (including `meta$s` = default candidate's `s`);
- scalar-style extras: `fit$W`, `fit$ranks`, `fit$ord`, `fit$dx`, `fit$ess`,
  `fit$neighbors`, `fit$psi`, `fit$nuisance`, `fit$y`;
- hidden `.model` state, unchanged.

Grid mode adds exactly one additive payload, `fit$tuning_grid` (name
**decided**), a list with components:

**`tuning_grid$candidates`** — one row per user candidate (Section 6.5).

**`tuning_grid$estimates`** — tidy long table, one row per
(candidate × target × component):

| column         | type      | semantics                                              |
|----------------|-----------|--------------------------------------------------------|
| `candidate_id` | character | key into `candidates`                                  |
| `kappa`        | numeric   | `NA` for manual grids                                  |
| `s`            | integer   | candidate kernel order                                 |
| `target_id`    | integer   | row index into `fit$x0` (1..m)                         |
| `component`    | character | `NA` for NP/CATE; moment-component name for PLQR       |
| `estimate`     | numeric   | point estimate                                         |
| `se`           | numeric   | fast-jackknife SE; `NA` if `se = FALSE`                |
| `status`       | character | `"ok"`, `"infeasible"`, `"failed"`                     |

Deliberately **no** CI columns: the wrappers take no coverage level, so
interval construction belongs to consumers (`plot_tuning(level =)`, summary
methods, the simulation layer with its `alpha`). Recording `se` is sufficient
and avoids baking a default level into fit objects.

**`tuning_grid$diagnostics`** — one row per candidate (not per target; see
Fact 4.4):

| column        | type      | semantics                                         |
|---------------|-----------|---------------------------------------------------|
| `candidate_id`| character |                                                   |
| `s`           | integer   |                                                   |
| `ess`         | numeric   | `1 / sum(tab^2)`; identical across targets        |
| `max_weight`  | numeric   | `max(abs(tab))`                                   |
| `top_share`   | numeric   | share of total absolute weight in the top `k` ranks (extend with further table functionals as needed) |
| `degenerate`  | logical   | `ess >= 0.9 * n` (the scalar-path warning criterion) |

Per-target, data-dependent diagnostics (CATE overlap, PLQR conditioning) are
added in Phase 4 as additional long tables; they must not block Phases 0–2.

**`tuning_grid$meta`** — compact: tuning `name`, `n`, `kx`, `k`, `method`,
`distance` name, `jackknife`, number of feasible/skipped candidates.

Do **not** store per-candidate weight matrices in the payload (Section 13).
The payload is tidy tables rather than a list of per-candidate fit objects;
Section 6.6 records that decision and specifies the accessor that provides the
list-of-fits view on demand.

### 6.5 Candidate schema (engine-facing)

The candidate table is where generality is cheap now and expensive later. The
DNN wrappers are tuned by one `s` derived from `kappa`, but the package
already ships TDNN wrappers with two scales `(s1, s2)`, and a `(k, s)` grid is
a plausible future request. Therefore the **engine** interprets only generic
columns and treats estimator-specific tuning columns as opaque labels:

| column         | type      | interpreted by engine | semantics                                             |
|----------------|-----------|:---------------------:|--------------------------------------------------------|
| `candidate_id` | character | yes                   | unique, stable row id (rules below)                    |
| `computed_key` | character | yes                   | deduplication key: equal keys share one computed result |
| `is_default`   | logical   | yes                   | exactly one `TRUE` row                                 |
| `feasible`     | logical   | yes                   | from the tuning validator                              |
| `reason`       | character | yes                   | `"ok"` or the infeasibility reason                     |
| `status`       | character | yes                   | filled by the engine: `"ok"`, `"infeasible"`, `"failed"` |
| `kappa`        | numeric   | no                    | DNN wrappers; `NA` for manual grids                    |
| `s`            | integer   | no*                   | DNN wrappers (*used by the DNN weight provider, not by the engine core) |

Candidate-id rules (must be deterministic, locale-safe, filename-safe):

- tuning grids with `kappa`: `"kappa_" %+% gsub("-", "m", gsub("\\.", "p",
  format(kappa, trim = TRUE, scientific = FALSE)))` — e.g. `kappa_0p75`,
  `kappa_2`;
- manual / `NA`-kappa rows: `sprintf("s_%d", s)`;
- collisions (identical formatted labels): append `"__2"`, `"__3"`, … in row
  order. Ids must be unique within one fit but need not be globally unique.

`computed_key` for the DNN wrappers is `sprintf("s_%d", s)`. A future TDNN
grid would use `sprintf("s1_%d_s2_%d", s1, s2)` with its own tuning columns;
the engine core needs no change for that.

For scalar calls the table has exactly one row (`is_default = TRUE`), so all
downstream code has a single path regardless of mode.

### 6.6 Per-candidate fit views: `candidate_fit()`

A natural alternative shape for a grid fit is a list with one entry per
candidate — one full result object per `s` (or per `kappa`, depending on the
preferred labeling). That shape is **rejected as the stored representation**
but **provided as an on-demand view**. Rejected as storage because:

- storing G result objects duplicates `.model` state (including the `n x p`
  training matrix), `x0`, and metadata G times, for no informational gain
  over the tidy tables plus the shared caches;
- it makes the primary fields ambiguous (`fit$tau` versus which entry?) and
  breaks the anchored-default backward-compatibility rule of Section 6.4;
- it invites per-candidate loops where vectorized payload consumption is the
  right tool (simulations, comparison plots).

What downstream code actually needs from the list-of-fits perspective is the
ability to hand *one candidate* to existing single-fit machinery —
`predict()`, `diagnostics()`, any single-fit plot. That is an accessor
problem, not a storage problem:

```r
candidate_fit(fit, candidate)              # -> np_result / cate_result / ...
candidate_fits(fit, candidates = "all")    # -> named list of the above
```

Specification:

- `candidate` is a `candidate_id` string or an integer row index into
  `fit$tuning_grid$candidates`;
- the returned object has the **same class and field layout as a scalar fit**
  at that candidate's `s`: primary estimate and `se` filled from the payload;
  `meta$s` and tuning metadata updated; `meta$candidate_id` records the
  provenance;
- candidate-dependent extras that the parent fit carried (`W`, `ess`,
  `neighbors`) are **recomputed on demand** from the cached ranks and the
  candidate's rank-weight table — O(nm), nothing stored;
- candidate-independent components (`x0`, `psi`/`y`/`nuisance`,
  `ranks`/`ord`/`dx` extras, `.model` state) are carried over unchanged;
  sharing `.model` is what makes the view cheap and makes `predict()` work;
- requesting the default candidate returns an object equal to the scalar-call
  result — a tested parity anchor (Section 15);
- calling `candidate_fit()` on a fit without a `tuning_grid` payload errors
  with the actionable hint to refit with `tuning_mode = "grid"`; requesting an
  infeasible or failed candidate errors with the recorded `reason`;
- `candidate_fits()` is a thin convenience wrapper (`lapply` over the singular
  accessor, list named by `candidate_id`) — it exists because the
  named-list-keyed-by-candidate shape is what users reach for in exploratory
  work, and it must never become a separately maintained code path.

`candidate_fit()` ships in Phase 1 alongside the payload: it makes every
existing S3 method and single-fit plot grid-capable on day one, before the
plotting family gains native multi-candidate displays in Phase 4.

### 6.7 Diagnostics S3 behavior

`diagnostics(fit)` keeps returning default-candidate diagnostics for both
scalar and grid fits (current behavior preserved). Phase 4 adds
`diagnostics(fit, candidates = c("default", "all"))` on the existing S3
methods; `"all"` requires a grid fit and errors otherwise with an actionable
message. No `diagnostics_grid()` function, ever. Until Phase 4, users read
`fit$tuning_grid$diagnostics` directly; that table is documented and stable
from Phase 1 on.

---

## 7. Warning and Error Policy in Grid Mode

Validation and feasibility (mostly existing, restated as the binding
contract):

- validate `k` and **all** candidate `s` values before any heavy computation;
- `k <= s` is a feasibility condition per candidate;
- `s < n` is a feasibility condition per candidate when delete-1 jackknife is
  requested (`se = TRUE` and `jackknife != "none"`);
- an infeasible **default** candidate is an error, exactly as today (same
  message text as the scalar path);
- infeasible **non-default** candidates are skipped: `status = "infeasible"`,
  `reason` preserved, estimate/SE rows emitted with `NA` values and the same
  status, and **one** summary warning per fit listing the skipped candidate
  ids and reasons;
- if no feasible candidate remains, error.

This policy deliberately treats "feasible for the point estimate but not for
the jackknife" as infeasible when `se = TRUE`. A candidate with `s >= n` is
useless for the paper's coverage-versus-`kappa` analysis, and a finer status
taxonomy ("estimated, no SE") can be added later without breaking the schema
(new `status` value) if a real use case appears.

Degenerate-weight warnings get an explicit grid-mode rule, because a
sensitivity grid *deliberately* includes small `s` values:

- the near-uniform-weights warning (`.dnn_warn_degenerate_weights()`) fires
  for the **default candidate only** — at most once per fit, exactly as
  today;
- non-default candidates never warn; the payload records per-candidate `ess`
  and `degenerate` flags (Section 6.4), where the information is useful for
  plotting the localization tradeoff instead of drowning the console.

Without this rule an eight-candidate grid could emit up to eight identical
warnings per fit and thousands per simulation run.

Error-message style follows design note 05 throughout: argument-specific,
expected/actual where helpful, `call. = FALSE`.

---

## 8. Internal Design

The implementation factors the current scalar path into reusable internal
helpers. Public wrappers stay readable and do not each implement their own
grid loop. All helpers below live in `R/user-dnn-core.R` unless noted.

### 8.1 Rank-weight tables in the `rank_weight` contract

Add one **optional** field to `rank_weight` specs:

```r
spec$table <- function(n) numeric(n)   # rank -> weight lookup for sample size n
```

Semantics and requirements:

- `table(n)` returns the weight assigned to rank `r` at position `r`, for a
  sample of size `n`, **already normalized** the way the spec's `eval()`
  would normalize realized weights (licensed by Fact 4.2);
- providers: `rank_weight_dnn_exact()`, `rank_weight_dnn_asymp()`, and both
  TDNN constructors (`R/wrank-dnn-wrapper.R`, `R/wrank-tdnn-wrapper.R`). The
  implementations are one-liners over existing code, e.g.
  `function(n) dnn_weights_exact_full_from_ranks(seq_len(n), k, s)$w_rank`
  (then defensively renormalized exactly as `eval_one` does). No new weight
  mathematics;
- `rank_weight_eval_batch()` (in `R/api-ranks.R`) uses the table when present:

  ```r
  tab <- spec$table(nrow(R))
  W   <- matrix(tab[R], nrow(R), ncol(R))
  ```

  falling back to the existing per-column loop when `spec$table` is absent
  (custom user specs). This single change removes the per-column table
  rebuild from **every** existing call site, scalar paths included;
- an internal accessor centralizes the capability check and output
  validation:

  ```r
  .rank_weight_table(spec, n)   # -> numeric(n) or NULL if unsupported
  ```

  Validation: numeric, length `n`, all finite; positive sum before
  normalization. Invalid tables are an internal-state error
  (`.abort_bad_state()`), since only package-provided specs supply them in
  Phase 0.

Feasibility interplay: building the jackknife table requires `s <= n - 1`,
which is the existing `s < n` condition; building the point-estimate table
requires `s <= n`, which the resolver validates.

### 8.2 Candidate resolver extension

Extend `.dnn_resolve_s()` — do not replace it — so the candidate tables it
already builds (both the tuning branch and the manual-`s` branch) carry the
full schema of Section 6.5. Concretely, its return value stays
`list(s, tuning, tuning_used)` with `tuning$grid$candidates` gaining
`candidate_id`, `computed_key`, `is_default`, `status` (initialized `"ok"` /
`"infeasible"` from `feasible`). Marking `is_default`: the row whose `s`
equals the default resolution and, among ties, the one whose `kappa` equals
the tuning object's default `kappa` when that column is available; if no grid
row matches the default `s` (possible with hand-written `eval_grid`
functions), **prepend** a synthetic default row flagged
`is_default = TRUE, kappa = NA` so grid mode always contains the anchor
candidate.

Scalar calls keep returning the same scalar `s` and a one-row candidate table.
No behavior change for any existing caller.

### 8.3 The grid engine

One internal helper drives all wrappers:

```r
.dnn_grid_eval(
  candidates,     # resolved candidate table (Section 6.5 schema)
  ranks,          # n x m integer rank matrix, cached from the default weights call
  weight_spec_fn, # function(candidate_row) -> rank_weight spec (e.g. exact/asymp with that s)
  estimate_fn,    # function(W, candidate_row) -> data.frame of estimate rows
                  #   (columns: target_id, component, estimate)
  se_fn = NULL,   # function(candidate_row, tab_jk) -> data.frame of SE rows
                  #   (columns: target_id, component, se); NULL when se = FALSE
  ctx = NULL
)
```

Engine algorithm:

```
results <- empty list keyed by computed_key
for key in unique(candidates$computed_key[candidates$feasible]):
    row  <- first candidate row with this key
    spec <- weight_spec_fn(row)
    tab  <- .rank_weight_table(spec, n)          # point-estimate table (size n)
    W    <- matrix(tab[ranks], n, m)             # transient; never stored
    est  <- estimate_fn(W, row)
    if se requested:
        tab_jk <- .rank_weight_table(spec, n - 1)
        se     <- se_fn(row, tab_jk)             # prefix-sum jackknife inside
    diag <- table functionals of tab (Section 4.4)
    results[[key]] <- list(est, se, diag)
map results back over ALL candidate rows (duplicates share results;
infeasible rows get NA rows with status = "infeasible")
assemble tuning_grid payload (candidates, estimates, diagnostics, meta)
```

Duplicate handling: computation is keyed on `computed_key`; result rows are
repeated per user candidate row so the reported grid matches the user's grid.
No warning for duplicates unless **all** candidates collapse to one unique
key (then warn once that the grid is effectively a single value).

Failure handling: a candidate whose evaluation throws is caught
(`tryCatch`), recorded with `status = "failed"` and the condition message in
`reason`, and does not abort the remaining candidates — unless the failing
candidate is the default, which re-throws (the scalar contract).

If `spec$table` is unavailable (not the case for any package spec, but the
engine must not crash on future extensions), the engine falls back to
`rank_weight_eval_batch()` for `W` and to the loop jackknife for SEs.

### 8.4 Point estimates and the default candidate

The wrapper flow in grid mode:

1. resolve inputs, pseudo-outcome/signal, and targets once, exactly as today;
2. call `weights_dnn()` **once for the default candidate** with
   `return = unique(c(return, "ranks"))` — the same call the scalar path
   already makes for the fast jackknife. This yields `W_default`, the cached
   `ranks`, and any user-requested extras, and it keeps the default candidate
   on the identical code path as a scalar fit (same warnings, same extras);
3. compute the default estimate and SE exactly as the scalar path does;
4. hand `ranks` plus the candidate table to `.dnn_grid_eval()` for the
   remaining candidates (the engine skips the default row's computation and
   splices in the already-computed default results, so the default is never
   evaluated twice);
5. attach `fit$tuning_grid`.

There is no new rank-cache abstraction anywhere in this flow: the cache is
the `ranks` matrix the wrapper already holds.

### 8.5 Prefix-sum fast jackknife

New internal helper implementing Fact 4.3:

```r
.dnn_loo_prefix(sig, ranks, tab_jk)
# sig     : n x q numeric signal matrix (q = 1 for NP/CATE)
# ranks   : n x m integer full-sample rank matrix
# tab_jk  : numeric(n - 1) pre-normalized rank-weight table for size n - 1
# returns : n x m matrix (q = 1) or n x q x m array (q > 1) of LOO estimates,
#           matching the dimensions .dnn_jackknife_se_from_loo() consumes
```

Per target `j` (hoist the sort out of any per-`s` loop when evaluating a
grid — the ordering depends only on `ranks[, j]`):

```r
ord <- order(ranks[, j])                 # observation index by rank position
for (c in seq_len(q)) {
  sr    <- sig[ord, c]                   # sig_(r)
  A     <- cumsum(sr * c(tab_jk, 0))     # A[t]; tab_jk[n] treated as 0
  B     <- rev(cumsum(rev(sr * c(0, tab_jk))))   # B[t]; tab_jk[0] treated as 0
  r     <- ranks[, j]
  loo   <- ifelse(r > 1L, A[pmax(r - 1L, 1L)], 0) +
           ifelse(r < n,  B[pmin(r + 1L, n)],  0)
  ...store loo...
}
```

(The `pmax`/`pmin` guards are load-bearing; see the pitfall in Section 4.3.)

Integration into `.dnn_jackknife_se_loo_matrix()`:

- on the `"fast"` path, when the rank-weight spec supplies a table, compute
  `loo` via `.dnn_loo_prefix()` and feed the **unchanged**
  `.dnn_jackknife_se_from_loo(loo, q, n, m)`;
- keep the current per-deletion loop, unexported, as the reference
  implementation and as the fallback for non-table specs;
- the `"loo"` path is untouched (Section 8.6).

Grid usage: `se_fn` in the engine closes over `sig` and the per-target `ord`
vectors (computed once) and calls `.dnn_loo_prefix()` with the candidate's
`tab_jk`.

The earlier draft's "deletions outer, candidates inner" loop nesting is
**superseded**: on the table path there is no per-deletion loop at all.

### 8.6 `jackknife = "loo"` in grid mode

For `"loo"`, each deletion changes the reference sample before ranking, so the
full-sample rank matrix cannot be reused. Grid mode remains *correct* for
`"loo"` via the conservative recomputation path (the engine's `se_fn` simply
calls the existing machinery per candidate), and the documentation must state
that grid mode is intended primarily for `jackknife = "fast"`. Never silently
switch the user's jackknife method.

### 8.7 Wrapper integration

**`np_dnn()` / `cate_dnn()` (Phase 1).** Signal is a single column (`y_vec` /
`psi`). `estimate_fn` is a one-liner over `.dnn_estimate_from_signal()`;
`component` is `NA` in the estimates table. All input resolution, pseudo-
outcome construction (AIPW / Robinson / direct), and target resolution happen
once, before the engine; the engine must never re-touch nuisance columns.

**`plqr_dnn()` (Phase 3).** Signal is the `n x q` moment-component matrix.
`estimate_fn` returns long rows with `component` set to the moment-component
names; the scalar public field `fit$moment` remains a matrix for backward
compatibility. The long form is the payload's canonical shape; this is why
`plqr_dnn()` is sequenced after the schema has settled in Phase 1.

**`gmm_dnn()` (later, out of scope for Phases 0–3).** `estimate_fn` wraps the
per-target GMM solve with the candidate's `W`; SEs come from the GMM
asymptotic covariance, not the jackknife, so `se_fn` is estimator-specific
there. Grid mode reuses model setup, ranks, and weights, but each candidate
pays a full optimization; the documentation must say so plainly. The engine
contract (Decision 3) is what makes this integration possible without a
rewrite — nothing in the engine assumes `crossprod`.

### 8.8 Nuisance and pseudo-outcome reuse

Inside `cate_dnn()` the pseudo-outcome is computed once — direct
`pseudo_outcome` validated once; AIPW nuisance columns extracted, validated,
and combined once; Robinson likewise. Grid evaluation begins only after `psi`
exists as a canonical numeric vector. In the simulation layer the same rule
scales up (Section 10): one DGP draw, one full-sample nuisance fit, one
cross-fit, one AIPW score per replicate and method, all candidates evaluated
on the same score and rank cache.

---

## 9. Plotting and Downstream Display

The plotting family is uniform — `plot_*(object, ...)` S3 generics over
result objects (`plot_estimates()`, `plot_weights()`, `plot_overlap()`,
`plot_influence()`, `plot_calibration()`, `plot_conditioning()`,
`plot_tuning()`, `plot_dashboard()`, plus `autoplot()`/`plot()` methods).
Grid fits must be supported natively, in a way that is consistent across the
family rather than invented per function.

### 9.1 Two-tier design

- **Comparison-native plots** answer a question *about the grid* (how does
  the estimate, CI, or weight profile change with `s`?). They consume
  `fit$tuning_grid` directly — no refits, no per-candidate object
  construction. Members: `plot_tuning()`, `plot_estimates()`,
  `plot_weights()` (rank-profile panel), the diagnostics autoplot, and a
  dashboard tuning panel.
- **Single-fit plots** answer a question about *one* fit and gain
  per-candidate access through `candidate_fit()` (Section 6.6):
  `plot_overlap(candidate_fit(fit, id))` works from Phase 1 on with zero
  plotting-code changes. Native side-by-side comparison displays for these
  plots are added only when a concrete display design exists (deferred;
  Section 17).

### 9.2 Shared `candidates =` argument

Comparison-native plot functions gain one shared argument:

```r
candidates = "default"   # default: current behavior, default candidate only
candidates = "all"       # all evaluated candidates
candidates = c("kappa_0p75", "kappa_1p5")   # explicit candidate ids
```

Rules:

- the default value `"default"` preserves today's output for every existing
  call — vdiffr baselines for scalar fits must not change;
- non-default values require a grid fit; on a scalar fit they error with the
  actionable hint to refit with `tuning_mode = "grid"` (exception:
  `plot_tuning()`, whose whole purpose is the sweep — for scalar fits it
  computes the grid itself via the internal engine, Section 9.4);
- infeasible and failed candidates are dropped from the display; when any
  were dropped, the plot subtitle says how many and why (one line, not one
  warning per candidate — consistent with Section 7);
- unknown candidate ids error, listing the available ids.

### 9.3 Visual encoding conventions

Candidates are an **ordered discrete series** (ordered by `s`, ascending),
not an unordered factor:

- default mapping is color/linetype with legend labels of the form
  `"kappa = 0.75 (s = 23)"` (falling back to `"s = 23"` for `NA`-kappa
  candidates), legend title `"candidate"`, following the palette conventions
  in `plot-common.R`;
- when a plot already binds color to something else, candidates map to
  facets instead;
- the default candidate is visually anchored (drawn last / emphasized line
  width) so the sensitivity reading is always "how far do the alternatives
  move from the default".

### 9.4 Per-function behavior

| function              | grid-native behavior                                                            | data source          | phase    |
|-----------------------|----------------------------------------------------------------------------------|----------------------|----------|
| `plot_tuning()`       | estimate ± CI versus `s` (its existing purpose); refit loop **deleted**          | payload; internal engine on `.model` for scalar fits | 4 |
| `plot_estimates()`    | profile style: one line + ribbon per candidate; pointrange style: dodged points  | `tuning_grid$estimates` | 4 |
| `plot_weights()`      | new comparison panel: weight-versus-rank curves overlaid for all candidates — these are exactly the rank-weight tables (Fact 4.4), O(nU), no `W` needed; target-specific weight displays go per-candidate via `candidate_fit()` | tables / accessor | 4 |
| diagnostics autoplot  | with `candidates = "all"`: `ess`, max weight, top shares versus `s`              | `tuning_grid$diagnostics` | 4 |
| `plot_dashboard()`    | adds a tuning panel when a grid payload is present; other panels stay default-candidate | payload | 4 |
| `plot_overlap()`      | per-candidate via `candidate_fit()`; native comparison display deferred          | accessor             | deferred |
| `plot_influence()`, `plot_calibration()`, `plot_conditioning()` | per-candidate via `candidate_fit()`; native comparison deferred | accessor | deferred |
| `plot_tdnn_scales()`  | out of scope until TDNN grids are scheduled                                      | —                    | deferred |

### 9.5 `plot_tuning()` specifics and the no-second-engine rule

`plot_tuning()` currently refits at each candidate through the public
wrappers, including a full delete-one jackknife per refit. Once the engine
exists, the package must not carry two grid paths, one of them slow:

- `plot_tuning()` consumes `fit$tuning_grid` directly when present (grid
  fits);
- for scalar fits it calls the internal engine on the stored `.model` state
  (which already carries `x_mat`, the signal, `distance`, `method`,
  `normalize`, `ctx`) instead of refitting through the public wrappers;
- until that integration lands (Phase 4), `plot_tuning()` keeps its refit
  loop but must not grow engine-like caching of its own;
- the same applies to any future summary or dashboard panel that sweeps `s`.

CI construction in plots uses each function's `level` argument on top of the
payload's `se` column, consistent with the no-CI-columns decision
(Section 6.4).

---

## 10. Simulation Integration

The simulation layer gains no new public grid-specific entry point.
`sim_cate_run()` and `sim_cate_one_rep()` gain the same `tuning_mode`
semantics as the estimators.

### 10.1 Method contract change (internal, wholesale)

Grid mode forces the internal method contract to change; make the change
wholesale rather than bolting candidate handling onto the named-vector
plumbing. New contract for every registered `.method_XXX()`:

```r
.method_XXX(dgp, x_cols, k, s, tuning, tuning_mode, alpha, ...)
  -> data.frame, one or more rows, columns:
     candidate_id (chr), kappa (dbl), s (int),
     tau_hat, tau0, bias, sq_err, se, cover, ci_len, neff (dbl),
     status (chr: "ok" | "failed" | "infeasible")
```

- **DNN methods** in scalar mode call the estimator with scalar `s` exactly
  as today and emit one row (`candidate_id = NA`, `kappa` from the tuning
  metadata when available, `s` = the scalar value, `status = "ok"`). In grid
  mode they call e.g. `cate_dnn(..., tuning = tuning, tuning_mode = "grid",
  se = TRUE, jackknife = "fast")` and emit one row per candidate directly
  from `fit$tuning_grid` (point estimates and SEs from `estimates`, `neff`
  from `diagnostics$ess`), including `status` for skipped candidates.
- **Non-DNN benchmarks** ignore `tuning`/`tuning_mode`, run **once** per
  replicate, and emit **one** row with `candidate_id = NA`, `kappa = NA`,
  `s = NA`.
- A method-level failure is caught in `sim_cate_one_rep()` as today and
  yields one all-`NA` metrics row with `status = "failed"`.

`.sim_metric()` gains a thin vectorized sibling (or is applied rowwise) that
returns data.frame rows instead of a named vector; coverage is computed from
`tau_hat`, `se`, `tau0`, and `alpha` per row exactly as now.

### 10.2 Raw output schema

`sim_cate_run()` output gains columns and changes one convention. Full column
set:

```
method, design, n, kx, k, tuning_source,
candidate_id, kappa, s,
seed,
tau_hat, tau0, bias, sq_err, se, cover, ci_len, neff,
status
```

**Deliberate breaking change, recorded here:** `s` moves from a cell-level
constant into the method-emitted rows, and benchmark rows carry `s = NA`
(they do not depend on `s`; the current convention of stamping the cell's `s`
onto benchmark rows was an artifact of the assembly code). Consumers —
`summarise_sim_results()` in `R/sim-cate-summary.R` and the paper repository's
table code — must group by `(method, design, n, kx, candidate_id)` treating
`NA` candidate ids as their own group. This lands in the same pass as the
contract change, with tests.

### 10.3 Flow through `sim_cate_one_rep()` / `sim_cate_run()`

- `sim_cate_run(tuning =, tuning_mode = "grid")` validates the combination
  (`tuning_mode = "grid"` requires `tuning`; `s_fn` remains mutually
  exclusive with `tuning` as today), then passes both through to
  `sim_cate_one_rep()`, which passes them to the methods. Candidate
  resolution happens inside the estimator per replicate — its cost is
  trivial (Section 3.2) and this guarantees the simulation and estimator
  paths can never disagree about the candidate set;
- the scalar path (`s_fn`, or `tuning` with `tuning_mode = "default"`)
  resolves one `s` per cell exactly as today — behavior-identical;
- the per-cell `verbose` message prints the candidate summary in grid mode
  (e.g. `s in {15, 23, 31, 47}` with the count of infeasible candidates)
  instead of a single `s`;
- fault tolerance is unchanged: total rep failure produces all-`NA` rows for
  every method — in grid mode, for every (method × candidate) expected row,
  which requires resolving the candidate set once in `sim_cate_run()` purely
  to shape failure rows (an acceptable, read-only second use of
  `tuning_eval()`).

### 10.4 Benchmarks are not refit and not repeated

Between the two reporting options in the earlier draft (repeat benchmark rows
per candidate for plotting convenience, or emit once), **emit once is
decided**. Repeated benchmark rows silently inflate replicate counts in any
naive downstream aggregation — a footgun for the paper tables. Broadcasting
benchmark results across candidates for a comparison plot is one join in the
summary/plot layer, where it belongs. `dnn_crossfit` must not refit folds per
`kappa`; the fold assignment and nuisance fits happen before the candidate
loop by construction (they precede the single `cate_dnn()` call).

### 10.5 RNG neutrality (contract, not accident)

The candidate loop must not consume RNG. DNN estimation is deterministic
given the data, and `dnn_crossfit`'s fold `sample()` happens before
estimation, so the property holds by construction — but it is elevated to a
tested contract: **per-candidate results in grid mode must match an external
per-`kappa` run with the same seed exactly** (not just to tolerance — the
arithmetic is identical). This is also the migration check for retiring the
external kappa-grid runners in the conditional Z-estimation repository
(Section 18).

### 10.6 `summarise_sim_results()`

Update the grouping to `(method, design, n, kx, candidate_id, kappa, s)` with
`NA`-safe grouping, keep `n_failed` semantics (now derivable from `status`),
and add nothing else in Phase 2. Plot-facing convenience (candidate-versus-
benchmark broadcasting) is a summary-layer helper, added only when the paper
plotting code needs it.

---

## 11. Runtime Implications

Let `G` = number of candidates, `U` = unique feasible `s` values after
rounding, `n` = sample size, `m` = number of targets, `q` = signal components,
`p` = covariate dimension.

The earlier draft framed nuisance refitting and distance/rank recomputation as
the main waste of an external sweep. Code review shows the dominant per-fit
cost at paper sizes is the fast-jackknife loop itself. Corrected accounting:

| piece                          | today (scalar fit)             | after Phase 0 (scalar)   | grid mode, per extra unique `s`     |
|--------------------------------|--------------------------------|--------------------------|-------------------------------------|
| distances + sort               | O(nmp) + O(mn log n)           | unchanged                | **0** (ranks cached)                |
| point-estimate weights         | O(nm) + m table builds         | O(nm) + **1** table      | O(nm) + 1 table                     |
| point estimate                 | O(nmq)                         | unchanged                | O(nmq)                              |
| fast jackknife                 | **O(n²m)** + nm table builds   | **O(nmq)** + 1 table + m sorts | O(nmq) (sorts reused)         |
| nuisances / pseudo-outcome     | once                           | once                     | **0**                               |
| DGP draw / benchmarks (sims)   | once per rep                   | once per rep             | **0**                               |

Consequently:

- Phase 0 is the largest single speedup available and applies to every scalar
  fit, independent of grid work;
- grid-mode standard errors on the fast path are nearly free per candidate;
  the earlier guidance that "grid SEs are not free" survives only for
  `jackknife = "loo"` (full recomputation per deletion) and for `gmm_dnn()`
  (a full solve per candidate);
- the performance target remains: never repeat work whose inputs are
  unchanged — but after Phase 0 the remaining per-candidate work is genuinely
  small, and no further micro-optimization of the candidate loop is
  warranted (see the no-parallelism rule in Section 5).

---

## 12. Duplicate and Infeasible Candidates (Summary)

Handled by the engine via `computed_key` and `status` (Sections 6.5, 7, 8.3):
user candidate rows are always preserved because the labels matter for
reporting; computation happens once per unique feasible key; duplicated rows
share results; infeasible non-default rows are skipped, recorded, and warned
about once in aggregate; an infeasible default errors. Duplicates warn only
when the whole grid collapses to a single unique value.

---

## 13. Interaction With Return Extras

The `return` argument stays additive and is never the trigger for expensive
grid evaluation:

- `tuning_mode = "default"` with `return = "tuning"`: returns default tuning
  metadata and candidate grid *metadata*, as today — metadata, not
  evaluation;
- `tuning_mode = "grid"`: evaluates feasible candidates and attaches the
  payload; `return = "tuning"` remains independent and unchanged;
- `return = "W"`, `"ranks"`, `"ord"`, `"dx"`, `"ess"`, `"neighbors"` in grid
  mode: all refer to the **default candidate**, unchanged.

If full per-candidate weights are ever exposed, that happens through an
explicit new extra (working name `return = "grid_weights"`), only after
memory behavior is measured on realistic simulation sizes, and without
changing the meaning of `return = "W"`. Deferred (Section 17).

---

## 14. Implementation Phases

Phases 0–2 are the paper-critical path. Each phase lands with its tests and
documentation in the same pass, per the repository working rules.

### Phase 0 — rank-weight tables and prefix-sum jackknife

Touches: `R/wrank-dnn-wrapper.R`, `R/wrank-tdnn-wrapper.R`, `R/api-ranks.R`
(`rank_weight_eval_batch()`), `R/user-dnn-core.R`
(`.rank_weight_table()`, `.dnn_loo_prefix()`,
`.dnn_jackknife_se_loo_matrix()` fast path), tests.

Deliverables: `table(n)` on the four package `rank_weight` constructors;
table-based `rank_weight_eval_batch()` with per-column fallback; prefix-sum
fast jackknife with the loop retained as reference/fallback.

Acceptance: equivalence tests of Section 15 pass; no public API change; no
user-visible behavior change beyond speed (numerical differences at floating-
point noise level are acceptable and covered by test tolerances).

### Phase 1 — grid mode in `np_dnn()` and `cate_dnn()`

Touches: `R/user-dnn-core.R` (`.dnn_resolve_s()` extension, `.dnn_grid_eval()`),
`R/user-dnn-npr.R`, `R/user-dnn-cate.R`, a new `R/user-dnn-candidates.R` for
the accessors, docs, tests.

Deliverables: `tuning_mode` argument; candidate schema; `tuning_grid` payload
with `candidates` / `estimates` / `diagnostics` / `meta`; grid-mode warning
and error policy; `candidate_fit()` and `candidate_fits()` accessors
(Section 6.6), which make every existing S3 method and single-fit plot
grid-capable immediately.

Acceptance: scalar calls byte-identical in behavior; grid default candidate
identical to scalar fit; per-candidate parity with separate scalar calls
(Section 15). `plqr_dnn()` deliberately excluded — its long-form output is
the most schema-sensitive piece and must not block the paper.

### Phase 2 — CATE simulation integration

Touches: `R/sim-cate-benchmarks.R` (method contract), `R/sim-cate-run.R`,
`R/sim-cate-summary.R`, tests.

Deliverables: data.frame method contract; `tuning_mode` on `sim_cate_run()` /
`sim_cate_one_rep()`; candidate-level raw rows with the Section 10.2 schema;
updated `summarise_sim_results()`.

Acceptance: RNG-neutrality parity test versus external per-`kappa` runs;
benchmark single-row rule; nuisance-reuse test with a counted mock. This is
the step that matters most for the conditional Z-estimation paper.

### Phase 3 — `plqr_dnn()` grid mode

Same engine, long-form callback with `component` rows, same payload
conventions, same tests adapted to `q > 1`.

### Phase 4 — diagnostics, summaries, and the plotting family

Touches: `R/user-dnn-*-diagnostics.R`, `R/plot-tuning.R`,
`R/plot-estimates.R`, `R/plot-weights.R`, `R/plot-dashboard.R`,
`R/plot-common.R` (candidate series encoding helpers), tests (including
vdiffr baselines for the new displays).

Deliverables, per Section 9:

- `diagnostics(fit, candidates = c("default", "all"))`;
- the shared `candidates =` argument and encoding conventions on the
  comparison-native plots;
- `plot_tuning()` consumes `tuning_grid` (grid fits) or the internal engine
  on `.model` state (scalar fits) — the refit loop is deleted;
- `plot_estimates()` candidate overlays; `plot_weights()` rank-profile
  comparison panel (tables only); dashboard tuning panel;
- summary methods surface estimate, SE, `ess`, max weight, top-weight
  shares, and (per-target, Phase-4-only) CATE overlap and PLQR conditioning
  diagnostics.

Acceptance: scalar-fit vdiffr baselines unchanged with default arguments;
no plot function computes estimator fits itself anymore.

### Phase 5 — weight-layer grid support (optional, unscheduled)

Consider extending `weights_dnn()` to accept a tuning object. Do not add
`weights_dnn_grid()`. The return-shape question (list of matrices vs. 3-d
array) must not block anything above.

---

## 15. Testing Requirements

Tests live at the public wrapper level, with shared-helper tests where the
public surface cannot reach a behavior. All equivalence tests compare against
retained reference implementations, not stored constants.

**Phase 0 (equivalence and contract):**

- `table(n)` output equals `eval()` output on the identity permutation, for
  exact, asymp (`normalize` both ways), and both TDNN specs, over a range of
  `(n, k, s)` including `k = s` and `s = n`;
- table-based `rank_weight_eval_batch()` equals the per-column path on random
  rank matrices (multiple targets);
- prefix-sum fast jackknife equals the per-deletion loop to numerical
  tolerance: `q = 1` and `q > 1`, `m = 1` and `m > 1`, small and large `s`,
  and covariates with heavy ties (exercising the deterministic tie-breaker);
- boundary regression: datasets where some observation has rank 1 / rank n at
  some target (the Section 4.3 pitfall);
- custom `rank_weight` specs without `table` still work end-to-end via the
  loop path (fallback contract).

**Phase 1 (grid mode):**

- scalar calls unchanged for `np_dnn()` and `cate_dnn()` (estimates, SEs,
  extras, warnings, errors);
- `tuning_mode = "grid"` default-candidate fields equal the scalar default
  call exactly;
- each candidate's estimate and fast-jackknife SE equals a separate scalar
  call with the same `s` exactly;
- duplicate `s` candidates computed once (assert via a counting shim on the
  table builder) and reported per candidate row;
- infeasible non-default candidates skipped and recorded; exactly one summary
  warning; infeasible default still errors with the scalar message;
- degenerate-weight warning fires only for the default candidate; payload
  `ess` present for all candidates;
- `return = "W"` refers to the default candidate in grid mode;
- `tuning_mode = "grid"` + explicit `s` errors; grid mode without `tuning`
  uses the paper default grid;
- candidate-id determinism (formatting, `NA`-kappa fallback, collision
  suffixes);
- `candidate_fit()` parity: the view for candidate `s` equals a scalar call
  with that `s` in all public fields (estimate, `se`, `meta$s`, recomputed
  `W`/`ess` when the parent carried them); the default-candidate view equals
  the grid fit's own primary fields;
- `candidate_fit()` errors: scalar fit without payload, unknown id,
  infeasible/failed candidate (message carries the recorded `reason`);
- `candidate_fits()` returns a list named by `candidate_id` and is
  element-wise identical to repeated `candidate_fit()` calls;
- `predict()` and `diagnostics()` work on a candidate view and agree with the
  same calls on the equivalent scalar fit.

**Phase 2 (simulation):**

- grid mode draws one dataset per replicate and returns stable candidate-level
  rows with the Section 10.2 schema;
- grid-mode candidate rows match external per-`kappa` runs seed-for-seed
  (exact equality — RNG-neutrality contract);
- benchmarks appear once per replicate with `candidate_id = NA` and are not
  refit per candidate;
- cross-fitted methods fit nuisances once per replicate: verify with a
  counted nuisance-fitting mock (the reuse contract, not wall-clock time);
- failure paths: whole-rep failure and per-method failure produce correctly
  shaped `NA`/`status` rows in both modes;
- `summarise_sim_results()` groups candidates correctly and keeps benchmark
  groups separate.

**Phase 4 (plotting and diagnostics):**

- scalar-fit vdiffr baselines unchanged for every plot function with default
  arguments;
- `candidates = "all"` / id vectors: correct series ordering by `s`, correct
  labels, default candidate anchored; errors on scalar fits (except
  `plot_tuning()`), unknown ids error listing available ids;
- dropped infeasible/failed candidates appear in the subtitle, not as
  warnings;
- `plot_tuning()` output from a grid fit's payload matches the output from
  the engine-on-`.model` path for an equivalent scalar fit;
- no plot function calls the public estimators (assert via a counting shim
  in tests, so refit loops cannot silently return).

---

## 16. Documentation Requirements

User-facing documentation must say:

- the ordinary estimator call still evaluates the default `s`;
- `tuning = tuning_dnn_paper(kappa_grid = ...)` defines candidate values, and
  `tuning_mode = "grid"` evaluates the feasible ones efficiently, reusing
  ranks, tables, and pseudo-outcomes;
- on the fast-jackknife path, per-candidate SEs are cheap;
  `jackknife = "loo"` and (later) `gmm_dnn()` remain expensive per candidate,
  and grid mode does not change that;
- full per-candidate weights are not returned by default;
- the `tuning_grid` payload schema (documented once, referenced from all
  wrappers);
- `candidate_fit()` / `candidate_fits()` as *the* way to treat a grid fit as
  a list of per-candidate fits, with the shared-`.model` and
  recomputed-extras semantics stated;
- the `candidates =` plotting convention (documented once in the plotting
  reference, cross-referenced from each comparison-native plot).

The CATE simulation documentation must state that grid mode is the preferred
way to run paper sensitivity checks once implemented, and must document the
raw-schema change of Section 10.2. External looping over `kappa` remains
valid for backwards compatibility but is no longer the recommended route.

`_pkgdown.yml` gains reference entries for `candidate_fit()` and
`candidate_fits()` in Phase 1; apart from those accessors, no new functions
are exported through Phase 3, and the wrappers' existing entries cover the new
argument.

---

## 17. Resolved and Remaining Design Questions

Resolved in this revision:

- **Argument name:** `tuning_mode = c("default", "grid")`, mirroring
  `tuning_eval(mode =)`; mode lives outside the tuning object (Section 6.1).
- **Grid payload name:** `fit$tuning_grid`, with the Section 6.4 schemas.
- **CI columns in the payload:** none; consumers construct intervals from
  `se` with their own level.
- **Candidate schema:** generic engine columns plus opaque estimator-specific
  tuning columns, keyed by `computed_key` (Section 6.5).
- **Rank cache:** no new helper; reuse the wrappers' existing
  `weights_dnn(..., return = c(..., "ranks"))` call (Section 8.4).
- **Jackknife strategy:** prefix-sum evaluation on the table path; loop kept
  as reference/fallback (Sections 4.3, 8.5).
- **Warning policy:** default-candidate-only degenerate warning; one summary
  warning for skipped candidates (Section 7).
- **Benchmark reporting in simulations:** one row per replicate,
  `candidate_id = NA`; broadcasting happens in the summary/plot layer
  (Section 10.4).
- **Simulation method contract:** data.frame rows, wholesale (Section 10.1),
  with the recorded raw-schema break (Section 10.2).
- **Grid-fit object shape:** tidy tables as stored payload; the
  list-of-fits-per-candidate shape is provided on demand by
  `candidate_fit()` / `candidate_fits()` and never stored (Section 6.6).
- **Plotting convention:** two-tier design; comparison-native plots gain a
  shared `candidates =` argument with the Section 9.3 encoding conventions;
  single-fit plots go per-candidate through the accessor (Section 9).

Still open, deliberately deferred:

- **Per-candidate weight storage** (`return = "grid_weights"`): decide only
  after memory behavior is measured on realistic simulation sizes (Phase 5
  at the earliest).
- **`diagnostics(fit, candidates =)`**: Phase 4; the payload already carries
  candidate diagnostics before then.
- **Vector-valued `s` in the public wrappers:** revisit only on a concrete
  need; must reuse the internal candidate schema if added (Section 6.3).
- **Finer `status` taxonomy** ("estimated, no SE"): only if a real use case
  appears (Section 7).
- **TDNN grid mode:** the candidate schema and engine support it by design
  (`computed_key`, opaque tuning columns); scheduling is a separate decision.
- **Native comparison displays for `plot_overlap()`, `plot_influence()`,
  `plot_calibration()`, `plot_conditioning()`:** per-candidate use works via
  `candidate_fit()` from Phase 1; a faceted or overlaid comparison design is
  added per function only when a concrete display need appears (Section 9.4).
- **`plot_tdnn_scales()` grid support:** out of scope until TDNN grids are
  scheduled.

---

## 18. Implication for the Conditional Z Paper Simulations

The current manuscript-side workaround loops over `kappa` externally (the
PowerShell kappa-grid runners in the paper repository). That remains
acceptable for pilots and for validating the expanded grid, but it is not the
right long-run route.

After Phases 0–2 land in RdistNN, the paper simulation runner should move to:

- one simulation cell per design, `n`, `kx`, method set, and replicate block;
- one `tuning_dnn_paper(kappa_grid = ...)` object passed into the DNN methods
  with `tuning_mode = "grid"`;
- candidate-level rows emitted by the package;
- no repeated nuisance fitting or rank ordering across `kappa` values.

The RNG-neutrality contract (Section 10.5) gives the migration check:
grid-mode candidate rows must reproduce the external per-`kappa` runs
seed-for-seed before the external runners are retired. The integrated route
also makes the paper output more internally coherent: candidate rows share
the same dataset draws and nuisance estimates by construction, which is the
right comparison when the goal is to isolate the effect of varying `s`.

---

## 19. Final Decision

The package treats tuning-grid evaluation as an integrated mode of the
existing DNN estimators and simulation functions, built on two internal
primitives that ship first on their own: rank-to-weight tables as part of the
`rank_weight` contract, and a prefix-sum evaluation of the fast delete-1
jackknife (verified equivalent to the current loop; regression-tested against
it).

The implementation preserves the scalar default path, reuses the canonical
matrix/vector internals and the ranks the wrappers already cache, computes
pseudo-outcomes once, evaluates unique candidate `s` values through
rank-weight tables keyed by a generic candidate schema, and returns
candidate-level results as an additive `tuning_grid` payload. The estimate
step is an estimator callback over weights, so the same engine serves
`np_dnn()`, `cate_dnn()`, `plqr_dnn()`, and eventually `gmm_dnn()` without
assuming a linear cross-product.

This gives users and paper simulations the performance benefits of grid reuse
without fragmenting the public API — and makes the scalar estimators
substantially faster along the way.
