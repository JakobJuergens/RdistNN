# 04 — Exact localization and the paired deletion engine

[Index](README.md) · [Previous](03_augmented_signals_and_validation.md) ·
[Next: inference](05_inference_confidence_sets_and_diagnostics.md)

**Status:** implementation design. The new owner is
`R/user-dnn-affine-ratio.R`, calling existing shared primitives.
No new compiled routine is needed to deliver the first complete feature.

## 1. Canonical computational inputs

After boundary validation, retain x_mat (n by k_x), x0_mat (m by k_x) and
signal = cbind(numerator,denominator) (n by 2). The original row order is
the canonical storage order. W and rank_matrix are n by m, and each rank
column is a permutation of 1..n. Do not pass an already sorted signal matrix
to the existing full-rank engine while retaining original-order ranks.

Use `.dnn_weights_with_extras()` once for the reporting targets and request
ranks internally whenever covariance is needed. Form raw component averages
with `.dnn_estimate_from_signal()`: its output is 2 by m, so transpose once
to the documented m by 2 public matrix. Preserve both signal roles and target
order through this conversion.

Compute neff = 1/colSums(W^2). For exact k=1 the numerical rank-weight values
depend only on n,s, so neff is the same across targets, up to accumulation
error. The physical radii, signal arrays receiving the weights and their
covariances can differ substantially.

## 2. Fixed geometry and coherent deletions

The covariance definition freezes the reported target, full nuisance state,
distance specification and order s. Full and retained weights use that same
order, not s-1, not a fresh call to the tuning rule at n-1.

Distance ties use input-row position. Retained order preserves the relative
order of surviving input rows. Row IDs are never substituted as a new sorting
key. Tests with tied points must compare this definition across engines.
Permutation invariance is asserted only with distinct distances or with an
explicitly preserved tie ordering.

The default l2 metric is deletion-invariant. For other existing distance
specifications, the algorithm still freezes the full distance/rank
realization. A distance whose setup recomputes scaling or learned geometry
on each retained sample would define a different pipeline jackknife.
Resolve any context once, retain the full ranks and derive deleted ranks
by removing rank r_j and subtracting one from ranks above it.

Consequently, the new wrapper must not blindly route arbitrary distance
specifications through the current `.dnn_jackknife_loo_array(...,
jackknife = "loo")` branch: that branch rebuilds weights from retained X.
For the new feature, use `.dnn_loo_loop(signal,ranks,rw)` for the explicit
rank-shift reference, or add a small adapter preserving the same frozen
context. For fixed l2, also compare against independently reranked retained
X in tests. Custom distance setup is invoked once, not once per deletion.

The software's fixed-array calculation is meaningful with a frozen
sample-derived metric. The primitive statistical theorem separately uses
its deterministic geometry and continuous-distance assumptions; the software
must not infer that those assumptions hold.

## 3. Exact weights and linear-time deletions

For sorted rank r,

$$
w_{n,s,1}=s/n,\qquad
w_{n,s,r+1}=w_{n,s,r}(n-r-s+1)/(n-r).
$$

Support ends at n-s+1. In particular, s does not bound the number of
nonzero weights. Use existing exact rank-weight constructors/table providers;
do not add another binomial-coefficient implementation to the wrapper.

Let v_r = w_(n-1,s,r) for retained sample size n-1 and f_r the two-column
signal in distance order. Deleting rank j gives

$$
L_j=\sum_{r<j}v_r f_r+\sum_{r>j}v_{r-1}f_r.
$$

Prefix sums of v_r f_r and shifted suffix sums of v_(r-1) f_r produce every
L_j in O(n) arithmetic after sorting. Existing
`.dnn_jackknife_loo_array(...,jackknife="fast",precomputed_ranks=...)`
already returns n by 2 by m for this signal dimension and dispatches through
`.dnn_loo_prefix()` and `jackknife_rank_full()`.

Use that path with method="exact",k=1,normalize=TRUE. Retain the explicit
rank-shift loop as the user-selected "loo" reference. The R and compiled
paths must implement the same arrays; compiled availability is an
optimization. Respect the existing compiled-dispatch option and restore
global options on exit in tests.

## 4. Joint covariance and scalar sandwich

Let d_R be the full average, L_Rj the retained average and
P_Rj = n d_R-(n-1)L_Rj. Define

$$
V_{RS}=\frac{1}{n(n-1)}
       \sum_j(P_{Rj}-\bar{P}_R)(P_{Sj}-\bar{P}_S).
$$

The equivalent stable implementation centers L directly and computes
`((n-1)/n) * crossprod(centered_loo)`. That is the existing
`.dnn_jackknife_vcov_from_loo(loo,q=2,n,m)` contract: it returns a list of
m 2 by 2 covariance matrices and component SEs.

At the full fitted ratio theta_hat, score pseudo-values satisfy
Q_j=theta_hat P_Dj-P_Yj. Therefore

$$
\widehat V_\theta=
\{V_{YY}+\widehat\theta^2V_{DD}
                    -2\widehat\theta V_{YD}\}/\widehat\Delta_D^2.
$$

Store this scalar as a named 1 by 1 matrix in theta_vcov[[target]].
Store neff[target] times that matrix in theta_omega. There is no extra
division by neff when obtaining se from theta_vcov.

For numerical stability, compute ratio variance primarily as a squared norm
of the centered deletion-score combination divided by the full denominator.
Use the three-entry quadratic formula as an independent consistency check,
not as the sole calculation when large terms nearly cancel.
Do not solve deleted ratios and do not fit deleted nuisance models.

## 5. Hand-computed acceptance fixture

Set X=matrix(1:4,ncol=1), x0=0, s=2, numerator=c(2,1,4,3),
denominator=c(1,-1,2,1). Then

| Quantity | Exact expected value |
| --- | --- |
| Full rank weights | (1/2,1/3,1/6,0) |
| Retained weights at n-1=3 | (2/3,1/3,0) |
| Components | (2,1/2) |
| theta | 4 |
| Numerator deletion averages | (2,8/3,5/3,5/3) |
| Denominator deletion averages | (0,4/3,1/3,1/3) |
| Numerator pseudo-values | (2,0,3,3) |
| Denominator pseudo-values | (2,-2,1,1) |
| Component covariance | rbind(c(1/2,1/2),c(1/2,3/4)) |
| Frozen score pseudo-values | (6,-8,1,1) |
| theta_vcov[[1]][1,1] | 34 |
| neff | 18/7 |
| theta_omega[[1]][1,1] | 612/7 |

The first deleted denominator is exactly zero. The specified frozen-score
covariance is nevertheless finite. A refitted/deleted-ratio jackknife cannot
serve as the reference for this fixture.

Complete-subset counting supplies further identities:
mean_j L_Rj=d_R, mean_j P_Rj=d_R and mean_j Q_j=0 when the empirical
equation is solved. Center numerically anyway. Independent small-array
enumeration verifies these identities without calling package weight or
deletion helpers.

## 6. Numerical failure handling

Compute theta only for representable component ratios with nonzero
denominators. Avoid a denominator-squared underflow by combining centered
deletion scores with division in a controlled, scaled norm calculation.
Check the resulting component matrices, scalar variance and neff for
finiteness before marking covariance usable.

A covariance matrix formed from a finite cross-product is positive
semidefinite mathematically. Tiny roundoff asymmetry can be symmetrized.
Do not silently replace a materially indefinite matrix or a negative
variance from inconsistent inputs by zero. An exactly zero score norm is a
computed zero variance, with the separate reporting status in Note 05.
Arithmetic overflow at a target is retained as numerical_failure.

Scale-aware tests use relative and absolute tolerances, for example
abs(a-b) <= atol+rtol*max(abs(a),abs(b)), with simple-array atol/rtol near
1e-12 and currency-scale tests chosen to preserve relative meaning.
Do not use one huge absolute tolerance that masks errors in small signals.

## 7. Time, storage and batching

For m targets, full sorting costs O(m n log n); the paired deletion map is
O(m n) for two fixed columns. Retaining W/ranks costs O(nm), and optional
deletions or pseudo-values cost O(2nm). Storing both raw and centered
versions unnecessarily multiplies memory.

Initial implementation can use the existing batch engine. A subsequent
internal target-chunk loop can cap peak working memory while preserving
the same result and optional-array shapes. Tests must compare chunked and
unchunked targets if that optimization is introduced. Benchmark runtime and
allocated memory separately; do not impose noisy wall-clock thresholds in CI.
There is no need to expand C++ before profiling shows the existing paired
engine is insufficient.
