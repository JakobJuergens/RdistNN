# Weighted accumulation

## Ownership and contract

`R/utils-accumulation.R` owns compensated reductions, weighted signal products
and retained-state rank deletions. `src/accumulation.cpp` provides thin native
reductions and `src/block_rank_loo.cpp` exposes original-signal block deletions;
`inst/include/compensated_arithmetic.hpp` owns standalone C++
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
tables. `R/distance-blocks.R` owns per-target group identification, weight
averaging and selected-policy validation. High-level construction averages assigned rank weights within exact
represented selected-metric distance groups, separately for each target.
Use compensated block mass reductions and preserve signed and unnormalized
masses, including zero entries beyond nominal support. Each deletion shrinks the deleted observation's block
before averaging the `n-1` table; it is not deletion from a permanently averaged
full-sample table. Retain union-of-deletion-support centering and expand that
union to entire intersecting target-specific distance groups before choosing an origin.

For contiguous rank blocks, accumulate original signals using three block
coefficients: the mean of the unshifted table entries for a deletion after the
block, the mean of the shifted entries for a deletion before it, and the mean
of the surviving entries for a deletion within it. Merge prefix and suffix
states strictly outside the deleted block with its actual survivor prefix and
suffix states. Add all four leading values, then their four corrections, before
rounding. Do not form rounded group-mean signals and subtract the removed row:
that can lose a finite survivor or turn signed-table variation into zero.
The R, C++ and Mata canonical sweeps perform O(n*m*q) work with O(n*q) scratch.
Caller-supplied ranks that interleave groups or reorder original members
within a group use direct deletion/regrouping within the selected implementation,
taking O(n^2*m*q) work. Forced native
execution uses the native block kernel, including that uncommon direct path.
Direct regrouping reduces each block's table entries in original observation
order, then reduces signal products in localization-rank order. Retain original
observation identities through native transport so rank sorting cannot change
the block-mass reduction order.

The rank-only engine stores the forward sum's leading value and correction in
two prefix working arrays, then streams the reverse accumulator. Merge leading prefix, leading suffix,
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
Form square-root weights before normalizing them, so an underflowed normalized
weight cannot erase a representable root contribution for unnormalized inputs.
Restore the SD directly, so its squared variance need not fit in a double.
Zero-weight extreme outcomes cannot set the scale. Signed TDNN weights do not
define a probability-distribution SD or quantiles; those diagnostics are missing.
Physical CATE/CLATE radii use the stable builtin Euclidean distance in original
coordinates and the same CDF rule, independently of the localization metric.
Stata preserves the localization engine marker around this diagnostic call.
The pure CATE subset score also uses the shared compensated reduction in scalar,
native batch and R batch routes.

Within-arm probability-weight ESS is `(sum(w))^2 / sum(w^2)`, independent of
a positive common scaling of that arm's weights. Keep the established ordinary
calculation when its squared mass is normal. If squared mass is subnormal or
zero despite positive weight, scale the arm weights before squaring and apply
the same concentration formula. An arm with exactly one positive weight has
ESS one; a genuinely zero-mass arm has unavailable ESS. Share this calculation
between arms and fresh or derived diagnostics wherever arm ESS is exposed;
this does not add an ESS field to a result schema that lacks one.
It changes no localization weights or estimator covariance. R's
[machine constants](https://stat.ethz.ch/R-manual/R-patched/library/base/html/zMachine.html)
and Mata's [numeric limits](https://www.stata.com/manuals/m-5mindouble.pdf)
distinguish the normal threshold from smaller representable values. The
exceptional scaling rule is the packages' application of that distinction.

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
also use the uniform identity when k=s, including the s=1 single-scale
component of TDNN; other cases retain host Beta routines.
Rounded interval endpoints can lose a tail mass or
produce different zero/nonzero support in R and Mata at subnormal scales.
Order-statistic averaging, signed TDNN blending and normalization can lose
additional tiny masses. This remains unresolved, including downstream effects
on estimates and deletion covariance for extreme signals. A component-only
probability fix would not establish correctness of the final mixed table.
Engine parity for a supplied table does not verify its construction, and
normalization cannot restore probability already lost. Exact and asymptotic
tables define different methods. Continuous values can also select different
observations at an exceptionally close CDF boundary.
Identical represented tables use the exact common CDF rule; no epsilon merges
distinct selected observations.

`R/utils-centering.R` and Mata `rdistnn_signal_level()` share the observed
origin and exact-subtraction policy for SD, covariance and deletions.
Deletion origins use the first eligible localization rank if subtractions
are within the common numeric range and exact (checked by a two-difference
residual). Fall back to the established same-sign endpoint/zero range policy;
use zero if subtraction would lose a separately represented component. Columns
whose entire supported amplitude is subnormal also use zero: translating them
can erase an unrepresentable weighted residual before restoring the origin.
When actual deletion weights are available, also use zero if the maximum
absolute coefficient times supported signal amplitude is subnormal. A nonzero
subnormal table mass, or a nonzero mass whose product with that amplitude is
subnormal, uses zero as well. This retains original separately rounded products
before any requested division by mass; it does not normalize weights first.
Exact zero mass alone does not disable centering. Maximum coefficient and mass
summaries are reused across components/targets; geometry-refitting references
collect actual maxima and the smallest nonzero mass in their existing support
pass. Ordinary tiny tails with normal coefficient amplitude and mass retain
their observed origin. These are underflow safeguards, not a guarantee of exact
product distributivity after arbitrary normal-table block averaging/cancellation.
Custom signed subnormal masses can still overflow during normalization even
when raw sums are finite; finite arithmetic does not guarantee partial recovery
of every aggregate. Keep those range limits separate from probability-table
construction and from covariance representability.
This exceptional policy is shared by fast and explicit deletion paths. These
choices protect both large local offsets with nonzero remote tails and mixed-sign
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
