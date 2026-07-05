# Integrated tuning-grid reuse for DNN wrappers

*Created: 2026-07-05.*

## Purpose

This note records the package-side design for evaluating a grid of DNN kernel
orders without rerunning the entire estimator workflow for each value of
`kappa` or `s`.

The immediate motivation is the corrected CATE simulation workflow for the
conditional Z-estimation paper. After correcting the paper-facing simulations to
use `k = 1`, the old grid of `kappa` values is less informative than it looked
under the mistaken `k = 5` setup. The manuscript simulation code now needs a
larger grid of `kappa` values to show the bias-variance and localization
tradeoff as `s` varies. Running that grid by launching a complete simulation
run per `kappa` is wasteful, because many expensive pieces do not depend on
`s`.

The same issue is not merely a simulation inconvenience. Empirical users will
naturally want sensitivity checks over plausible `kappa` values. The package
should make that workflow ergonomic and fast enough to use routinely.

The central design decision is:

> Integrate grid evaluation into the existing functions and tuning machinery.
> Do not create new exported functions whose names encode grid-specific
> signatures, such as `cate_dnn_grid()`, `np_dnn_kappa_grid()`, or
> `sim_cate_kappa_grid_run()`.

This is consistent with the existing package design notes: user-facing
estimators should have a stable grammar, tuning should live in the tuning
abstraction, and optional behavior should be additive rather than forcing users
onto a parallel API surface.

---

## Current State

The package already has most of the pieces needed for an efficient integrated
grid path.

### Tuning layer

`tuning_dnn_paper()` already stores both:

- a default `kappa`, used to choose the primary `s`;
- a `kappa_grid`, used by `tuning_eval(..., mode = "grid")`.

`tuning_eval()` already returns candidate metadata with:

- `kappa`,
- `s`,
- `feasible`,
- `reason`.

The missing part is not candidate generation. The missing part is that the main
estimator wrappers use the grid only as metadata and still evaluate the
estimator at one resolved scalar `s`.

### Estimator wrappers

The main wrappers:

- `np_dnn()`,
- `cate_dnn()`,
- `plqr_dnn()`,
- and, eventually, `gmm_dnn()`,

all follow the same broad pattern:

1. validate the public interface;
2. resolve `data=` or matrix inputs into canonical matrices and vectors;
3. resolve a scalar `s` through `.dnn_resolve_s()`;
4. resolve targets;
5. compute DNN weights;
6. combine weights with a signal or moment;
7. optionally compute delete-1 jackknife standard errors;
8. attach optional return extras and hidden model state.

This structure is good. The grid path should reuse it rather than adding a
separate user-facing wrapper family.

### Weight and rank layer

`weights_matrix()` is already ranker-agnostic. It can consume:

- a distance spec or name;
- precomputed ranks;
- precomputed orders;
- precomputed distances;
- custom rank providers.

The DNN-specific weights are rank weights. For fixed `X`, `X0`, distance, and
tie-breaking rule, the rank matrix does not depend on `s`. Only the
rank-to-weight map depends on `s`.

This is exactly the separation we need:

- compute distances, orders, and ranks once;
- evaluate `rank_weight_dnn_exact(k, s)` or `rank_weight_dnn_asymp(k, s)` once
  per candidate `s`.

### Fast jackknife path

The `"fast"` jackknife path already uses precomputed full-sample ranks. That is
the right direction, but it is currently scalar in `s`. If users run a grid by
calling the full estimator repeatedly, the package recomputes the same
full-sample ranks and the same delete-1 rank updates once per candidate.

The integrated grid path should instead:

- compute full-sample ranks once;
- for each deletion, compute the leave-one-out rank update once;
- evaluate all candidate rank-weight maps on that updated rank matrix.

The jackknife still depends on `s`, so the grid is not free. But the rank work
and the upstream estimator work should not be repeated.

---

## What Should Not Be Built

Do not add one-off exported functions such as:

- `cate_dnn_grid()`;
- `np_dnn_grid()`;
- `plqr_dnn_grid()`;
- `cate_dnn_kappa_grid()`;
- `sim_cate_kappa_grid_run()`;
- `weights_dnn_kappa_grid()`.

Those names solve the immediate problem by duplicating the public surface. They
would also create several long-run problems:

- every estimator would need a scalar function and a grid function;
- documentation and examples would split across two APIs;
- users would have to learn which arguments belong in the scalar call and which
  belong in the grid call;
- simulation code would drift away from the estimator API used by empirical
  users;
- future wrappers, such as a general conditional Z-estimation front end, would
  inherit unnecessary naming debt.

The package already has a tuning abstraction. The grid is a tuning behavior,
not a separate estimator.

---

## Public API Direction

The preferred public design is to extend the existing wrappers with one shared
control for whether the tuning grid is only recorded or actively evaluated.

Working name:

```r
tuning_mode = c("default", "grid")
```

The exact name can be revisited before implementation, but the behavior should
be common across all high-level DNN wrappers.

### Default behavior

The default remains scalar:

```r
fit <- cate_dnn(
  data = dat,
  x_cols = x_cols,
  y_col = "Y",
  t_col = "D",
  m0_col = "m0_hat",
  m1_col = "m1_hat",
  e_col = "e_hat",
  x0 = x0
)
```

This call resolves the default `s` and returns the same class and primary fields
as today.

### Grid behavior

Users who want the candidate grid should call the same function:

```r
fit <- cate_dnn(
  data = dat,
  x_cols = x_cols,
  y_col = "Y",
  t_col = "D",
  m0_col = "m0_hat",
  m1_col = "m1_hat",
  e_col = "e_hat",
  x0 = x0,
  tuning = tuning_dnn_paper(
    kappa = 1,
    kappa_grid = c(0.5, 0.75, 1, 1.25, 1.5, 2, 2.5, 3)
  ),
  tuning_mode = "grid"
)
```

No estimator-specific `kappa_grid` argument is needed. The grid belongs inside
the tuning object. This keeps `cate_dnn()`, `np_dnn()`, `plqr_dnn()`, and future
wrappers on the same grammar.

### Explicit `s`

For the first implementation pass, keep explicit `s` scalar in the public
wrappers. Users who want a grid should use `tuning`.

Supporting `s = c(15, 23, 31, 47)` can be considered later, but it should not be
the initial public path. A raw vector of `s` values lacks the `kappa` labels,
formula metadata, feasibility metadata, and default-candidate semantics already
available through the tuning API.

If vector-valued `s` is added later, it should be treated as a manual tuning
grid and should follow the same internal candidate schema as `tuning` grids. It
should not produce a separate return shape.

### Return shape

Existing primary fields should remain anchored at the default candidate:

- `fit$tau` for `cate_dnn()`;
- `fit$mu` for `np_dnn()`;
- `fit$moment` for `plqr_dnn()`;
- `fit$se`;
- `fit$meta`;
- scalar-style extras such as `fit$W`, `fit$ranks`, `fit$dx`, and `fit$ess`.

When `tuning_mode = "grid"`, add a structured grid payload, for example:

```r
fit$tuning_grid
```

The field name can be finalized during implementation, but the payload should
be additive. It should not replace the scalar primary fields.

Recommended contents:

- `candidates`: candidate-level metadata from `tuning_eval(..., mode = "grid")`,
  augmented with candidate ids, default flags, duplicate-`s` flags, and
  evaluation status;
- `estimates`: tidy long table with candidate, target, estimate, SE, CI, and
  diagnostics;
- `diagnostics`: candidate-level DNN diagnostics when available;
- optional compact internal metadata needed by `summary()` or `diagnostics()`.

Do not store every candidate weight matrix by default. For realistic grids this
can be large, especially for many targets. The grid path should compute and
store diagnostics by default; full per-candidate weights should be opt-in and
carefully named if they are ever exposed.

### Diagnostics behavior

`diagnostics(fit)` should continue to work for scalar fits. For grid fits, it
should be able to return either:

- the default-candidate diagnostics, preserving current behavior; or
- all candidate diagnostics when explicitly requested.

This can be done later with a small argument on the existing S3 method, for
example `candidates = c("default", "all")`. Avoid a separate
`diagnostics_grid()` function.

---

## Internal Design

The implementation should factor the current scalar path into reusable internal
helpers. The public wrappers should stay readable and should not each implement
their own grid loop.

### Candidate resolver

Replace or extend `.dnn_resolve_s()` with an internal helper that can return a
full candidate table while preserving scalar behavior.

Working helper:

```r
.dnn_resolve_s_candidates(
  s,
  tuning,
  n,
  p,
  k,
  se,
  jackknife,
  ctx,
  tuning_mode = c("default", "grid")
)
```

Recommended output:

```r
list(
  default_s = 31L,
  tuning = tuning_payload,
  tuning_used = TRUE,
  candidates = data.frame(
    candidate_id = c("kappa_0p5", "kappa_0p75", "kappa_1"),
    kappa = c(0.5, 0.75, 1),
    s = c(15L, 23L, 31L),
    is_default = c(FALSE, FALSE, TRUE),
    feasible = c(TRUE, TRUE, TRUE),
    reason = c("ok", "ok", "ok"),
    duplicate_s = c(FALSE, FALSE, FALSE),
    stringsAsFactors = FALSE
  )
)
```

For scalar calls, `candidates` should contain one row. That lets downstream code
use one path regardless of mode.

Validation rules:

- validate `k` and all candidate `s` values before heavy computation;
- keep `k <= s` as an explicit feasibility condition;
- keep `s < n` as an explicit feasibility condition when delete-1 jackknife is
  requested;
- error if the default candidate is infeasible, as today;
- for grid mode, mark non-default infeasible candidates and skip evaluation
  rather than aborting the whole fit, unless no feasible candidate remains.

Skipping infeasible non-default candidates is important for practical grids. A
large `kappa` value may be useful for large `n` but infeasible for a small pilot
cell. The result should record that fact instead of failing the whole workflow.

### Rank cache

Add an internal helper that computes the rank information once:

```r
.dnn_compute_rank_cache(
  x_mat,
  x0_mat,
  distance,
  ctx,
  return = c("ranks", "ord", "dx")
)
```

This helper should use the same distance resolution, batch-distance evaluation,
and deterministic tie-breaking as `weights_matrix()`. It should return:

- `ranks`;
- optionally `ord`;
- optionally `dx`;
- distance name and any useful metadata.

The cache is independent of `s`. It is valid for all candidates that share:

- the same `X`;
- the same `X0` or `x0`;
- the same distance specification;
- the same distance setup context;
- the same deterministic tie rule.

This cache should be treated as internal infrastructure. Public users can
already pass precomputed ranks through `weights_matrix()` when they need that
level of control.

### Weight evaluation from ranks

Add a scalar internal helper:

```r
.dnn_weights_from_rank_cache(
  rank_cache,
  method,
  k,
  s,
  normalize,
  ctx
)
```

This should:

1. construct the appropriate rank-weight spec;
2. call `rank_weight_eval_batch()` on cached ranks;
3. apply the same degenerate-weight warning policy as the current scalar path;
4. return `W` plus simple metadata.

Then add a vectorized helper:

```r
.dnn_evaluate_weight_grid(
  rank_cache,
  candidates,
  method,
  k,
  normalize,
  ctx
)
```

This helper should compute each unique feasible `s` once. If two `kappa`
values round to the same `s`, reuse the computed weight matrix and repeat the
result rows so the reported grid still matches the user's candidate grid.

### Estimate evaluation

The signal-combination part is independent of how the signal was created.
Create one shared internal path that takes:

- a signal matrix with `n` rows;
- a rank cache;
- candidate metadata;
- estimator identity metadata;
- a function for assembling estimator-specific result rows.

For scalar CATE and NP this is simply:

```r
estimate <- crossprod(signal, W)
```

For PLQR, the signal is a matrix of moment components, so the result should be
stored in long form:

- `candidate_id`,
- `kappa`,
- `s`,
- `target_id`,
- `component`,
- `estimate`,
- `se`,
- diagnostics.

That long form should be used inside the grid payload even if the scalar public
field remains a matrix for backward compatibility.

### Fast jackknife over a grid

For `jackknife = "fast"`, the efficient loop should be organized around
deletions first and candidate `s` values second.

Current repeated-call behavior is effectively:

```r
for (s in s_grid) {
  compute full ranks
  for (i in seq_len(n)) {
    compute delete-i ranks
    compute W_i(s)
    compute estimate_i(s)
  }
}
```

The integrated behavior should be:

```r
compute full ranks
for (i in seq_len(n)) {
  compute delete-i ranks once
  for (s in unique_feasible_s_grid) {
    compute W_i(s)
    compute estimate_i(s)
  }
}
```

This avoids repeated full-rank work and repeated delete-1 rank updates. It does
not avoid the fact that the leave-one-out estimate changes with `s`; the
rank-weight evaluation and cross-products still need to happen per unique `s`.

The helper should return a table or array that can be mapped back to candidate
rows. For scalar CATE and NP this can be a matrix with dimensions:

- candidate by target.

For PLQR this naturally becomes:

- component by candidate by target,

or a long table immediately after computation.

### Leave-one-out jackknife over a grid

For `jackknife = "loo"`, the package cannot reuse the full-sample rank matrix
in the same way, because each deletion changes the reference sample before the
distance ranking is computed.

Grid mode should still be correct for `"loo"`, but it will remain much more
expensive. The first implementation can support `"loo"` by reusing the
outer candidate schema while doing the conservative recomputation path. It is
acceptable for the documentation to say that grid mode is intended primarily
for `jackknife = "fast"`.

Do not silently change the user's jackknife method from `"loo"` to `"fast"`.

### Nuisance and pseudo-outcome reuse

This is especially important for `cate_dnn()` and the simulation harness.

Inside `cate_dnn()`, the pseudo-outcome should be computed once:

- direct `pseudo_outcome` path: validate once;
- AIPW nuisance-column path: extract and validate nuisances once, then compute
  AIPW once;
- Robinson path: extract and validate nuisances once, then compute Robinson
  pseudo-outcome once.

Grid evaluation should begin only after `psi` is available as a canonical
numeric vector. The grid path should never recompute nuisance extraction or
pseudo-outcomes per candidate.

In `sim_cate_run()`, the improvement is larger:

- draw the DGP once per replicate;
- fit full-sample nuisances once per replicate and method;
- fit cross-fitted nuisances once per replicate and method;
- compute the AIPW score once per replicate and method;
- evaluate all candidate `s` values using the same score and rank cache.

This is the major reason to integrate the grid into the package rather than
running an external PowerShell or R wrapper that launches separate simulations
per `kappa`.

---

## Simulation Integration

The simulation layer should not gain a new public grid-specific entry point.
Instead, extend `sim_cate_run()` and `sim_cate_one_rep()` with the same tuning
semantics as the estimators.

Recommended direction:

- keep `sim_cate_run()` as the public simulation entry point;
- keep `tuning` as the way to provide the paper grid;
- add the same mode concept used by the estimator wrappers;
- return candidate-level rows when grid mode is active.

The raw simulation output should include at least:

- `method`;
- `design`;
- `n`;
- `kx`;
- `k`;
- `candidate_id`;
- `kappa`;
- `s`;
- `tuning_source`;
- `seed`;
- `tau_hat`;
- `tau0`;
- `bias`;
- `sq_err`;
- `se`;
- `cover`;
- `ci_len`;
- `neff`;
- `status` or `failed` metadata.

For non-DNN benchmarks such as `ols_interact`, `local_constant`, `causal_forest`,
or `llf`, there are two reasonable reporting choices:

1. run them once per replicate and record one non-grid row;
2. run them once per replicate and repeat the same result across candidate rows
   for plotting convenience.

The first choice is cleaner statistically. The second choice can simplify
paper-facing tables that compare every DNN candidate against the same
benchmarks. The important rule is that these benchmarks must not be refit once
per DNN candidate.

For DNN methods, the simulation method implementations should be refactored so
the expensive steps happen before the candidate loop. For example,
`dnn_crossfit` should not refit folds once per `kappa`.

---

## Runtime Implications

Let:

- `G` be the number of candidate `s` values;
- `U` be the number of unique feasible `s` values after rounding;
- `n` be sample size;
- `m` be number of targets;
- `p` be covariate dimension.

The current external sweep repeats almost everything `G` times:

- DGP draw in simulations;
- nuisance fitting;
- pseudo-outcome construction;
- distance evaluation;
- rank ordering;
- weight construction;
- point estimate;
- jackknife SE;
- global benchmarks.

The integrated path should reduce the repeated work to the pieces that truly
depend on `s`:

- DNN rank-weight evaluation;
- point-estimate cross-products;
- leave-one-out estimate arrays for the jackknife;
- candidate-specific diagnostics.

Distance evaluation and rank ordering should be approximately one-time per
fit. Nuisance fitting and pseudo-outcome construction should also be one-time
per fit or per simulation method. This matters because the expanded paper grid
has more candidate values than before, and future users may try even larger
grids.

The grid will still increase runtime. It should not be advertised as free. A
reasonable expectation is:

- point estimates over a grid should be much cheaper than `G` full refits;
- fast-jackknife grid SEs should be materially cheaper than `G` full
  jackknife refits, but still grow roughly with `U` in the rank-weight and
  cross-product pieces;
- cross-fitted CATE simulations should see the largest practical improvement,
  because nuisance fitting is no longer repeated per candidate.

The performance target is not zero marginal cost. The target is to avoid
repeating work whose inputs are unchanged.

---

## Duplicate and Infeasible Candidates

The package must handle duplicate `s` values cleanly.

Duplicate `s` values can occur when different `kappa` values round to the same
integer, especially in small pilots. The result should preserve the user's
candidate rows because the labels matter for reporting, but the computation
should use unique feasible `s` values internally.

Recommended behavior:

- keep all candidate rows in `fit$tuning_grid$candidates`;
- add `duplicate_s = TRUE/FALSE`;
- add `computed_s_id` or equivalent to map duplicate candidate rows to a
  unique computed value;
- repeat estimates and diagnostics for duplicated candidate rows;
- do not emit a warning for duplicates unless all candidates collapse to one
  unique value.

Infeasible candidates should also be represented explicitly.

Recommended behavior in grid mode:

- evaluate feasible candidates;
- skip infeasible non-default candidates;
- record `status = "infeasible"` and the feasibility reason;
- warn once if any candidates were skipped;
- error only if the default candidate is infeasible or no feasible candidates
  are available.

This is more useful than aborting the whole fit because one large `kappa`
violates `s < n` for a small pilot sample.

---

## Interaction With Return Extras

The existing `return` argument should remain additive and should not be used as
the primary trigger for expensive grid evaluation.

In particular, avoid designs where a user accidentally computes a large grid
only because they requested `return = "tuning"`. Returning tuning metadata is
not the same as evaluating all candidate fits.

Recommended behavior:

- `tuning_mode = "default"` and `return = "tuning"`: return default tuning
  metadata and candidate grid metadata, as today;
- `tuning_mode = "grid"`: evaluate feasible candidates and attach the grid
  payload;
- `return = "W"` in grid mode: keep returning the default-candidate `W` unless
  a later design explicitly adds an opt-in way to store all candidate weights.

If full per-candidate weights are added later, prefer an explicit extra such as
`return = "grid_weights"` rather than changing the meaning of `return = "W"`.
That addition should come only after memory behavior is tested.

---

## Implementation Phases

### Phase 1: shared internal grid engine

Add internal helpers for:

- resolving scalar and grid candidates into one schema;
- computing and validating a rank cache once;
- evaluating DNN weights from cached ranks for one `s`;
- evaluating estimates and diagnostics across candidate rows;
- computing fast jackknife SEs across unique feasible `s` values.

Do this behind the existing scalar wrappers first, with scalar behavior
unchanged.

### Phase 2: integrate `np_dnn()`, `cate_dnn()`, and `plqr_dnn()`

Add the shared mode argument to the existing wrappers.

Each wrapper should:

1. validate and canonicalize inputs once;
2. build the estimator signal once;
3. resolve targets once;
4. call the shared grid engine;
5. populate scalar primary fields from the default candidate;
6. attach the grid payload only when grid mode is requested.

The three wrappers should use the same naming, validation, and output
conventions.

### Phase 3: integrate CATE simulations

Refactor `sim_cate_one_rep()` and `sim_cate_run()` so grid mode evaluates all
candidate DNN `s` values on the same simulated dataset, nuisance fits, and
pseudo-outcomes.

This is the implementation step that will matter most for the conditional
Z-estimation paper simulation runs.

### Phase 4: extend diagnostics and summaries

Teach `diagnostics()` and summary methods how to expose the grid payload.

The first version should focus on high-value columns:

- estimate;
- SE;
- CI length;
- `n_eff`;
- maximum weight;
- top-weight shares;
- overlap diagnostics for CATE;
- local moment conditioning diagnostics for PLQR when available.

### Phase 5: consider public weight-layer grid support

After estimator-grid semantics settle, consider extending `weights_dnn()` to
accept a vector of `s` values or a tuning object. Do not add
`weights_dnn_grid()`.

The weight-layer return shape needs more care because a natural output is either
a list of weight matrices or a three-dimensional array. That decision should
not block estimator and simulation integration.

---

## Testing Requirements

The implementation should include tests at the public wrapper level and at the
shared helper level where necessary.

Minimum tests:

- scalar calls are unchanged for `np_dnn()`, `cate_dnn()`, and `plqr_dnn()`;
- `tuning_mode = "grid"` returns the same default estimate as a scalar default
  call;
- each candidate estimate matches a separate scalar call with the same `s`;
- duplicate `s` candidates are computed once and reported for each candidate
  row;
- infeasible non-default candidates are skipped and recorded;
- infeasible default candidate still errors;
- fast-jackknife grid SEs match repeated scalar fast-jackknife calls;
- `return = "W"` still refers to the default candidate;
- diagnostics for grid fits include candidate identifiers and `s`;
- `sim_cate_run()` grid mode draws one dataset per replicate and returns stable
  candidate-level rows;
- cross-fitted simulation methods do not refit nuisances per candidate.

For runtime-sensitive behavior, add at least one small regression test with a
mock or counted nuisance-fitting function if that can be done without making
the public API artificial. The test should verify the intended reuse contract,
not exact wall-clock time.

---

## Documentation Requirements

User-facing documentation should say:

- the ordinary estimator call still evaluates the default `s`;
- `tuning = tuning_dnn_paper(kappa_grid = ...)` defines candidate values;
- the shared grid mode evaluates feasible candidates efficiently;
- the grid reuses ranks and pseudo-outcomes where possible;
- jackknife SEs still depend on `s`, so grid SEs are not free;
- full per-candidate weights are not returned by default.

The CATE simulation documentation should also state that grid mode is the
preferred way to run paper sensitivity checks once it is implemented. External
looping over `kappa` should remain valid for backwards compatibility, but it
should no longer be the recommended route.

---

## Open Design Questions

These should be resolved before implementation starts.

### Argument name

The working name is `tuning_mode`. Alternatives include:

- `tuning_eval = c("default", "grid")`;
- `evaluate_tuning = c("default", "grid")`;
- `sensitivity = FALSE/TRUE`.

`tuning_mode` is currently preferred because it is explicit, short, and
parallel to the existing `tuning` argument.

### Grid payload name

The working name is `tuning_grid`. Alternatives include:

- `sensitivity`;
- `tuning_results`;
- `grid`;
- `candidates`.

`tuning_grid` is currently preferred because it is specific and avoids implying
that the object is a generic sensitivity analysis over all estimator inputs.

### Diagnostics method argument

A later diagnostics extension needs a way to choose default-candidate
diagnostics or all candidate diagnostics. A possible argument is:

```r
diagnostics(fit, candidates = c("default", "all"))
```

This is not necessary for the first estimator-grid implementation if the grid
payload already contains candidate diagnostics.

### Weight storage

The package should not store all candidate weight matrices by default. The open
question is whether an explicit return extra for all candidate weights is worth
supporting. This should be decided only after checking memory behavior on
realistic simulation sizes.

---

## Implication for the Conditional Z Paper Simulations

The current manuscript-side workaround is to loop over `kappa` externally. That
is acceptable for pilots and for validating the expanded grid, but it is not the
right long-run simulation route.

After the integrated grid path lands in RdistNN, the paper simulation runner
should move toward:

- one simulation cell per design, `n`, `kx`, method set, and replicate block;
- one `tuning_dnn_paper(kappa_grid = ...)` object passed into the DNN methods;
- candidate-level rows emitted by the package;
- no repeated nuisance fitting or rank ordering across `kappa` values.

This will also make the paper output more internally coherent. Candidate rows
will share the same dataset draws and nuisance estimates by construction,
which is the right comparison when the goal is to isolate the effect of varying
`s`.

---

## Final Decision

The package should treat tuning-grid evaluation as an integrated mode of the
existing DNN estimators and simulation functions.

The implementation should preserve the scalar default path, reuse canonical
matrix/vector internals, compute ranks and pseudo-outcomes once, evaluate unique
candidate `s` values through the rank-weight layer, and return candidate-level
results as an additive payload.

This gives users and paper simulations the performance benefits of grid reuse
without fragmenting the public API.
