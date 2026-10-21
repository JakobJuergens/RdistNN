# Weighted accumulation

## Ownership and contract

`R/utils-accumulation.R` owns compensated reductions, weighted signal products
and retained-state rank deletions. `src/accumulation.cpp` provides thin native
reductions; `inst/include/compensated_arithmetic.hpp` owns standalone C++
arithmetic and rank-table deletion passes. The header is mirrored in the Stata
companion. Its `_rdistnn_accumulation.ado` provides native Mata quad reductions
and an analogous retained-state fallback. Estimators continue to own scores,
validation, fitted geometry, statistical normalization and result assembly.
Covariance centering/range policy stays in its existing shared module.

Rank-weight mass and weighted signal means use Neumaier compensated addition
in input-row order. Each weighted product is rounded to double before addition;
native code explicitly prevents contraction into a fused multiply-add.
Weighted signal products return the same component-by-target shape and
dimnames as the previous cross-product calculation. The package diagnostic
switch `options(RdistNN.use_compiled = FALSE)` forces the R implementations.
Mata reduces explicitly rounded product vectors with native `quadsum(..., 1)`.
Shared numeric fixtures check the arithmetic; input and product-overflow missing
values propagate per component/target, without omitting valid rows elsewhere.

Mahalanobis whitening is another sensitive signed reduction. Its matrix product
uses the same ordered double products and compensated dot reduction for every
row, independently of scalar/batch BLAS kernel selection. The R fallback
vectorizes independent row states; the native path reuses the cross-product
kernel with a transposed input. Mata reuses its weighted-signal reduction.
This preserves projection identity when rows are subsetted or targets are
batched, while SVD factorization and its numerical rank cutoff remain separate.

Covariance keeps vectorized cross-products for ordinary entries. When an
off-diagonal's absolute value is at most `sqrt(eps)` times its two column norms,
recompute it with compensated, separately rounded products. Retain the
two-difference residual of mean subtraction and include its cross terms in
that reduction: centering a large value can otherwise erase a represented small
covariance. Full covariance requests use one additional signal-sized residual
matrix; marginal-only requests do not allocate it. Diagonal range policy and
statistical normalization stay unchanged. This correction addresses summation
and centering error, not error already incurred when multiplying two doubles.

## Deletions and workspace

The arbitrary-table deletion engine supports nonnegative DNN and signed TDNN
tables. Retain the existing union-of-deletion-support centering. Store the
forward sum's leading value and correction in the two prefix working arrays,
then stream the reverse accumulator. Merge leading prefix, leading suffix,
prefix correction, suffix correction and restored offset in that order, before
rounding the result. R/Mata process the reverse sweep in blocks of at most
2,048 rows, bounding vectorized merge temporaries. No additional full deletion
array is needed to retain the correction.

The finite-support target-by-neighbor jackknife adapter reuses this primitive;
target aggregation retains baseline corrections and adds replacement terms
separately. Its slow reference explicitly rebuilds affected neighbor sets. Repeated indices
retain the finite-support fast API's additive single-position convention; use
distinct indices per row to represent deletion of an observation.
The generic frozen-rank reference also uses compensated weighted means.
Candidate fits and saved prediction call the same estimator adapters as fits.

Keep vectorized covariance and general matrix algebra. Compensation is applied
to sensitive weighted reductions/deletions, not indiscriminately to every
addition, count, regression solve or BLAS call.

## Precision and verification

This is ordinary finite double-precision computation, without an arbitrary-
precision runtime dependency. A single correction does not guarantee correctly
rounded answers for every severely ill-conditioned sum. Compensation cannot
recover multiplication/subtraction error already present in the inputs and
does not replace range normalization or the host's representability limits.
R and Stata retain their documented inference and numerical-failure rules.

Check both accuracy and engine parity: literal binary cancellation identities,
explicit deletions, signed tables, constants, zero support, unused extremes,
small samples, reverse-buffer boundaries, reordered/duplicate targets, covariance,
candidates and restored prediction. Force each engine and require native
availability in native checks. Discrete decisions remain exact; continuous
library-dependent results use quantity-specific tolerances. Benchmark the
production fallback's runtime and peak memory, not only a scalar prototype.

## Weighted distributions and diagnostics

`R/utils-weighted-distribution.R` owns the diagnostic weighted CDF and population
SD. Quantiles stably sort by value then input position. Compute total and each
prefix using serial Neumaier double addition in that sorted order; select the
first index whose rounded compensated prefix divided by the rounded compensated
total is at least `prob`. No boundary epsilon is used. At zero probability the
first sorted observation is retained even if its weight is zero. Nonnegative
finite weights and finite positive mass are required. The native header owns
the same two-scan index rule; Mata spells out the same updates. CDF selection
is linear after sorting, with scalar scan workspace. Literal represented DNN
tables are checked against rounded exact rational totals/prefixes in the
fixture generator; this does not claim correct rounding for every input.

Weighted outcome SD uses positive-weight observations, an exact binary scale,
an observed origin and a scaled norm of square-root-weighted centered outcomes.
Restore the SD directly, so its squared variance need not fit in a double.
Zero-weight extreme outcomes cannot set the scale. Signed TDNN weights do not
define a probability-distribution SD or quantiles; those diagnostics are missing.
Physical CATE/CLATE radii use the stable builtin Euclidean distance in original
coordinates and the same CDF rule, independently of the localization metric.
Stata preserves the localization engine marker around this diagnostic call.
The pure CATE subset score also uses the shared compensated reduction in scalar,
native batch and R batch routes.

## Exact rank-table generation and numerical origins

`R/utils-rank-tables.R` owns exact rank tables; `src/rank_weights.cpp` is a thin
adapter to mirrored `inst/include/rank_weight_tables.hpp`. Mata uses
`_rdistnn_rank_weights.ado`. For each subset order statistic, anchor its PMF
at a mode, recurse outwards using combinatorial adjacent-rank ratios, normalize
with serial Neumaier addition in rank order, and average the first k PMFs with
retained correction state. Normalize the mixture once. Keep mantissa and binary
exponent separate during recurrence so subnormal multiplication cannot stick
at the smallest double. Workspace is O(n), work O(n + k*(n-s+1)); k=s is exactly uniform.
This avoids host-specific log-gamma rounding in exact tables. Asymptotic tables
retain host Beta routines: their continuous values can differ by rounding and
can still select different observations at an exceptionally close CDF boundary.
Identical represented tables use the exact common CDF rule; no epsilon merges
distinct selected observations.

`R/utils-centering.R` and Mata `rdistnn_signal_level()` share the observed
origin and exact-subtraction policy for SD, covariance and deletions.
Deletion origins use the first eligible localization rank if subtractions
are within the common numeric range and exact (checked by a two-difference
residual). Fall back to the established same-sign endpoint/zero range policy;
use zero if subtraction would lose a separately represented component. This
protects both large local offsets with nonzero remote tails and mixed-sign
cancellation fixtures. Never restore the origin before covariance centering.
Geometry-refitting reference deletions choose the observation with largest
total absolute deletion weight. Ordinary SD/covariance paths apply the same
exact-subtraction safeguard before observed-origin mean centering, after any
required binary normalization. Raw restored doubles can
still round away sub-ULP variation; inference uses the shifted values.

## Maintenance checks

The source modules own their helpers directly; jackknife callers use
`.numeric_signal_level()` without a second forwarding name. Verify mirrored
headers with the companion's `tools/ci/package_hygiene.py --r-root PATH`.
The test-suite maps group these contracts by arithmetic/geometry/consumer
behavior rather than by the development session that exposed a bug.
