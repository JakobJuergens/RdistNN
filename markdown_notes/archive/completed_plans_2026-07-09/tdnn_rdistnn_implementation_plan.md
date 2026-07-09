> **Status (2026-07-09): implemented and archived.** The TDNN rank-weight
> layer, `np_tdnn()`, and the fast jackknife described here are all
> implemented. The 2026-07-08 namespace pruning later made the raw
> `tdnn_rank_table_*` and `dnn_weights_tdnn_*_full_from_ranks` helpers
> internal; the public surface is `rank_weight_tdnn_*()` / `weights_tdnn*()`.

# TDNN for `RdistNN`: weight layer and fast jackknife plan

Date: 2026-03-30

## Bottom line

The clean way to add the Demirkaya two-scale DNN (TDNN) scheme to `RdistNN` is to add it first at the `rank_weight` layer, not at the distance or ranking layer.

In package terms, TDNN is just a signed linear combination of two single-scale DNN rank-weight tables evaluated on the same rank matrix:

- `W_tdnn = w1^* W_dnn(s1) + w2^* W_dnn(s2)`
- with
  - `w1^* = 1 / (1 - (s1 / s2)^(-2 / d_loc))`
  - `w2^* = 1 - w1^*`

Here `d_loc` is the localization dimension from the Demirkaya bias expansion. In the simplest Euclidean case this is just `ncol(X)`, but I would expose it explicitly at the low level.

For the **exact** `k = 1` TDNN weights, the rank-order sequence is

- `a_j^TDNN(n, s1, s2) = w1^* * choose(n - j, s1 - 1) / choose(n, s1) + w2^* * choose(n - j, s2 - 1) / choose(n, s2)`

so the exact support cutoff is

- `j <= n - s1 + 1`

not `n - s2 + 1`.

The reason is simple: the embedded `s1` piece can still assign mass to ranks between `n - s2 + 2` and `n - s1 + 1`. A smaller `s1`-subset can exclude many closer points even when the enclosing `s2`-subset cannot.

The main implementation conclusion for variance estimation is:

- phase 1: reuse the existing fast delete-1 rank-update path once TDNN rank weights exist;
- phase 2: if we want a truly strong rank-based jackknife implementation, add a new large-support / full-rank-sequence jackknife engine.

For the **exact** TDNN rule, there is a genuine compact-support representation with cutoff `n - s1 + 1`, so an analogue of the current finite-support jackknife machinery is possible. The asymptotic TDNN rule is different: like the current asymptotic DNN approximation, it is effectively full-rank.

## Relevant existing package pieces

The main files I would build on are:

- `R/wrank-dnn-wrapper.R`
- `R/wrank-dnn_exact.R`
- `R/wrank-dnn_approx.R`
- `R/api-weights.R`
- `R/user-dnn-core.R`
- `R/jackknife-rank.R`
- `src/jackknife_rank.cpp`

The main theory references on our side are:

- `../JK_GenUStat/Paper/Parts/Main_TDNN.tex`
- `../JK_GenUStat/Paper/Parts/Supp_TDNN_Notation.tex`
- `../JK_GenUStat/Paper/Parts/Supp_TDNN_HajekDominance.tex`

The useful structural takeaway from the jackknife repo is the same one we want computationally:

- TDNN should be treated as an inheritance object built from two single-scale DNN pieces;
- the right primitive is the induced observation-level signed weight/localizer, not a new ranking rule.

## What I would implement first

### 1. Add paper-faithful TDNN rank weights

I would start with the Demirkaya construction itself, meaning `k = 1`.

That matters because the paper's bias-cancellation weights use the `s^(-2/d)` expansion for the single-scale DNN, and that is the object the TDNN theory and the jackknife repo are discussing.

So my initial low-level API would look something like:

- `rank_weight_tdnn_exact(s1, s2, dim)`
- `rank_weight_tdnn_asymp(s1, s2, dim, normalize = TRUE)`
- `weights_tdnn_exact(...)`
- `weights_tdnn_asymp(...)`
- `weights_tdnn(...)`

I would not initially expose a general `k` parameter for TDNN unless we are prepared to separately justify the corresponding bias-correction formula. If we later want a generalized multiscale bagged-kNN layer, that should be a separate extension.

I would also keep the exact and asymptotic TDNN layers visibly separate:

- exact TDNN: compact support ending at `n - s1 + 1`;
- asymptotic TDNN: effectively full rank support.

### 2. Make the weight layer table-based instead of recomputing by column

This is important for TDNN, but it also improves the current DNN implementation.

Right now the DNN rank-weight specs conceptually work column-by-column from a full rank vector. But for a fixed sample size `n`, the DNN weight sequence in **rank order** depends only on:

- `n`
- `s`
- and the exact/asymptotic method

It does **not** depend on the target point once the ranks are known.

So the right internal primitive is not "compute weights from one full rank permutation" but instead:

- `dnn_rank_table_exact(n, k, s)`
- `dnn_rank_table_asymp(n, k, s, normalize = TRUE)`

returning the rank-order weight table `w_rank[1:n]`.

Then a batch rank matrix `R` can be evaluated by simple indexing:

- `W[, t] = w_rank[R[, t]]`

TDNN then becomes trivial:

- `tdnn_rank_table_exact(n, s1, s2, dim) = w1^* dnn_rank_table_exact(n, 1, s1) + w2^* dnn_rank_table_exact(n, 1, s2)`
- same idea for the asymptotic table.

This is the right level of abstraction because the TDNN weights are just signed linear combinations of two DNN tables.

## User-facing API options

There are two reasonable directions.

### Option A. Minimal additive API

Add TDNN-specific wrappers alongside the existing DNN wrappers:

- `weights_tdnn_*`
- `np_tdnn()`
- `cate_tdnn()`
- `plqr_tdnn()`

This is low risk and easy to explain.

### Option B. Generalize the higher-level wrappers

The stronger design is to let the user-facing estimators accept an arbitrary `rank_weight`, with DNN and TDNN just being convenience constructors.

That would move the package toward:

- one localization/ranking layer,
- one weight layer,
- one estimator layer,

which is already how `weights_matrix()` is structured.

I suspect this is the better long-run design, but it is a larger API change. If the goal is to land TDNN quickly, I would do Option A first.

## Fast jackknife: what is easy now

### Phase 1: reuse the existing fast delete-1 rank-update path

Once TDNN exists as a rank-weight spec, we can already get a decent fast delete-1 jackknife by copying the current logic in `R/user-dnn-core.R`:

1. compute full-sample ranks once;
2. for delete-one observation `i`, update ranks by decrementing entries larger than the deleted rank;
3. evaluate the TDNN rank-weight spec on the updated `(n - 1) x m` rank matrix;
4. aggregate the signal on the kept sample.

This is exact, simple, and already avoids recomputing distances.

In other words, the first TDNN jackknife implementation does **not** need new theory. It just needs:

- a TDNN rank-weight spec;
- a TDNN analogue of `.dnn_jackknife_se_loo_matrix()`, or better, a generic rank-weighted version of that helper.

This should already be reasonably performant if the weight tables are cached by `n`.

That path should work for both:

- exact TDNN;
- asymptotic TDNN.

So this is still the right low-risk first version.

## Fast jackknife: the exact TDNN shortcut

For the **exact** TDNN rule, there is a stronger statement than I wrote in the first draft: we can likely build a direct analogue of the current finite-support rank jackknife machinery.

Let

- `L = n - s1 + 1`

and write the exact TDNN rank-order coefficients as

- `a_1, ..., a_L`

with `a_r = a_r^TDNN(n, s1, s2)`.

Then delete-one only affects ranks from the deleted position onward, and the leave-one-out rule only needs the rank-`L + 1` pull-in term. That is exactly the structural pattern used by `jackknife_fast_rank()`:

- support length `L`;
- only the first `L + 1` nearest indices are needed;
- the delete-one shift is local in rank space.

So for the exact TDNN case, an implementation path parallel to the existing finite-support jackknife utility is available:

1. build the exact TDNN rank coefficients `a_1, ..., a_L`;
2. pass the first `L + 1` neighbors per target;
3. reuse the same delete-one shift algebra.

This is the strongest case for a package-native `fast_tdnn_exact` implementation with the same style of fallback structure as the current DNN code.

## Fast jackknife: what I think we actually want

### Phase 2: add a large-support / full-rank-sequence jackknife engine

If we want a jackknife implementation that is genuinely strong for the asymptotic TDNN rule, and also robust when the exact-support cutoff `L = n - s1 + 1` is close to `n`, I would still add a new engine for large-support / full-rank rules.

The key observation is that for TDNN we do not need to rebuild the whole leave-one-out weight matrix from scratch.

For a fixed target point, suppose the signal is in rank order:

- `z_(1), ..., z_(n)`

and let:

- `w_n[j]` be the TDNN weight for rank `j` with sample size `n`;
- `w_nm1[j]` be the TDNN weight for rank `j` with sample size `n - 1`.

Then the baseline estimator is:

- `theta0 = sum_{j = 1}^n w_n[j] z_(j)`

If we delete the observation currently sitting at rank `r`, the leave-one-out estimator is:

- `theta_{-r} = sum_{j < r} w_nm1[j] z_(j) + sum_{j = r}^{n - 1} w_nm1[j] z_(j + 1)`

This is the whole algebraic core.

From here the efficient implementation is straightforward:

1. For each target, build the signal in rank order using `ord`.
2. Precompute
   - prefix sums `P[r] = sum_{j = 1}^r w_nm1[j] z_(j)`
   - shifted suffix sums `S[r] = sum_{j = r}^{n - 1} w_nm1[j] z_(j + 1)`
3. Then
   - `theta_{-r} = P[r - 1] + S[r]`
4. Map rank `r` back to deleted observation `i` using the rank or order matrix.

This gives all delete-one aggregates for one target in `O(n q)` time after sorting, where `q` is the signal dimension.

That is the right analogue of a fast rank-based jackknife for TDNN.

### Why this is better than the current rank-update loop

The current high-level fast path still loops over deletions and reevaluates weights on an updated rank matrix.

The full-support jackknife engine above instead works directly in rank order and computes **all** delete-one values for a target from one pair of cumulative-sum passes.

So the expected complexity is much closer to:

- `O(n m q)` after ranks/orders are known,

rather than a repeated delete-one reweighting cost.

This is also generic:

- DNN can use it;
- TDNN can use it;
- any other full-support rank-only weight sequence can use it.

## Why I would not only rely on `jackknife_fast_rank()`

The existing `jackknife_fast_rank()` engine is designed around a finite-support structure:

- finite support `a_1, ..., a_s`;
- only the first `s + 1` neighbor indices matter;
- deleting one observation only causes a local shift inside that compact support.

Exact TDNN does fit that template, with support length `s = n - s1 + 1`.

But two caveats remain:

- in the usual `s1 = o(n)` regime, that support length is still almost the whole sample;
- the asymptotic TDNN rule does not have that compact cutoff.

So my revised view is:

- yes, an exact-TDNN analogue of the current finite-support jackknife is plausible and probably worth doing;
- but I would still want the larger-support engine because it is the cleaner abstraction for asymptotic TDNN and for the practically common regime where `n - s1 + 1` is very large.

A separate larger-support engine is still cleaner and more reusable.

## Theory-to-code translation from the paper and the jackknife repo

The paper proves jackknife consistency through the U-statistic representation of TDNN.

For implementation, though, even Demirkaya et al. note that the WNN/L-stat representation is the computationally useful one. That is exactly the representation we already have infrastructure for in `RdistNN`.

The jackknife repo pushes the same structural message:

- the TDNN object inherits from the single-scale DNN inputs;
- the real new object is the signed effective weight/localizer.

So in package terms the guiding principle should be:

- build TDNN as a linear combination of single-scale rank tables;
- keep all speedups on the rank-order side, not on the U-statistic side.

## Practical caveats

### Signed weights

TDNN weights can be negative.

That is fine for estimation and jackknife variance estimation, but it affects interpretation of some existing diagnostics:

- `sum(W) = 1` should still hold;
- `1 / sum(W^2)` can still be reported as a concentration diagnostic;
- but it should no longer be described as a literal positive-weight effective sample size without a caveat.

### Dimension argument

The bias-correction exponent uses the localization dimension.

For ordinary Euclidean ranking this is naturally `ncol(X)`, but for custom distances or embeddings the correct dimension may not be the raw number of columns. So I would:

- require `dim` explicitly in low-level TDNN weight constructors;
- let high-level wrappers default to `ncol(X)` while exposing an override.

### Delete-one legality

As with the current DNN jackknife logic, we need:

- `s2 < n`

For the exact compact-support TDNN representation, the relevant neighbor depth is:

- `L + 1 = n - s1 + 2`

so the "finite support" implementation is only materially sparse when `s1` is itself close to `n`.

so that the same scale can still be used after deleting one observation.

### Tuning should be separate from the core implementation

The Demirkaya sign-change and weighted LOOCV ideas are useful, but I would not couple them to the first TDNN implementation.

First land:

- deterministic fixed `(s1, s2, dim)` weighting,
- fast jackknife parity,
- user-facing wrappers.

Then add tuning helpers.

## Recommended implementation order

1. Add rank-table helpers for DNN exact/asymptotic weights and cache them by `n`.
2. Add TDNN rank-table helpers and `rank_weight_tdnn_exact/asymp`.
3. Add `weights_tdnn_*` wrappers.
4. Add an exact-TDNN jackknife path that mirrors the current finite-support structure, using support cutoff `n - s1 + 1`.
5. Add parity tests against brute-force leave-one-out on small examples.
6. Keep the generic rank-update fallback for both exact and asymptotic TDNN.
7. Only after that, add a larger-support C++ jackknife engine based on prefix/suffix sums in rank order.
8. Benchmark the exact finite-support path, the generic fallback, and the larger-support engine against one another.

## Files I would expect to touch

Likely new files:

- `R/wrank-tdnn-wrapper.R`
- `R/wrank-tdnn_exact.R` or a more generic rank-table helper file
- `R/user-tdnn-core.R` or a generalized replacement for parts of `R/user-dnn-core.R`
- `src/jackknife_rank_full.cpp`

Likely modified files:

- `R/api-weights.R`
- `R/user-dnn-core.R`
- `NAMESPACE`
- `README.md`
- relevant man pages
- `tests/testthat/*`

## My preferred design

If I had to choose one concrete path after this correction:

- implement TDNN first as a paper-faithful `k = 1` rank-weight object;
- make the rank-weight layer table-based and cache by sample size;
- for exact TDNN, add a direct finite-support jackknife path with cutoff `n - s1 + 1`;
- keep the existing delete-one rank-update strategy as the universal fallback;
- then add a dedicated larger-support rank-sequence C++ jackknife kernel if the benchmarks justify it.

That gives a low-risk first implementation and a clear path to a genuinely fast variance routine.
