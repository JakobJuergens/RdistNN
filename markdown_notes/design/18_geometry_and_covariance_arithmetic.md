# Geometry and covariance arithmetic

## Ownership

`R/dist_mahal.R` fits Mahalanobis geometry once and delegates projected norms
to the existing L2 distance specification. `R/utils-covariance.R` owns numeric
centering, range normalization, covariance assembly and safe symmetrization.
Estimator adapters supply statistical factors, choose marginal or joint
covariance, and assemble their existing classes and computational statuses.
Keep these responsibilities shared across NP, TDNN, CATE, paired ratios and
experimental moment estimators. Deletion engines continue to return aligned,
support-shifted deletion arrays; they do not implement separate covariance
policies. Prediction and candidates use the same adapters as initial fits.

The Stata counterparts are `_rdistnn_geometry.ado` and
`_rdistnn_covariance.ado`. Geometry fitting stays in R/Mata. The common native
or fallback L2 norm kernel evaluates transformed distances. Covariance assembly
stays in the host layer regardless of which deletion engine ran.

## Fitted Mahalanobis metric

Center each training column safely on its own binary scale around an observed
value, then choose a common scale from centered variation. Huge constant
levels must not erase small varying columns. For normalized centered design
`Z = U D V'`, retain values
`d > max(n,p) * 2^-52 * max(d)` and cache columns
`V[,keep] / (d[keep]/sqrt(n-1))`. Euclidean distances between projected rows
implement the sample-covariance Moore-Penrose metric without forming a Gram
matrix or inverse. Exact null directions contribute zero, including distances
to targets outside the training span. A rank-zero design gives zero distances.

This explicit design-scale cutoff replaces the old solve-then-truncate policy
on covariance singular values. New fits can intentionally choose different
neighbors when the old calculation lost a resolvable direction. Tests check
the resulting ranks and estimates against independent latent-coordinate
formulas. The cutoff is shared across packages, but SVD fitting and projection
rounding can change represented distance equality or ordering even for
well-conditioned full-rank designs far from the cutoff. Ranks, exact-distance
groups, weights and estimates therefore need not agree across packages or
platforms at this boundary. Applying the same norm after projection does not
make independently fitted geometry bitwise equivalent. Group averaging uses
each fit's computed distances, with no tolerance or additional rounding; it
cannot reconcile different computed distances.

For example, the centered rows `(2,1)`, `(-2,-1)`, `(1,2)`, `(-1,-2)` are
mathematically equidistant from zero, with squared Mahalanobis distance `3/2`.
A one-ULP difference in computed distances can split this group even though
the design has full rank. With outcomes `(0,0,10,10)`, `k=1` and `s=3`,
weights `(.5,.5,0,0)` give estimate zero, while uniform weights give five.

Contexts retain scale, a subtraction-safe raw origin, factor, version and rank.
The origin may be zero when an observed-value shift would erase represented
differences. Evaluation
must respect this state for subsets, duplicate/reordered targets and restored
fits. Supplied contexts must use the current whitening format; older inverse
contexts require refitting. Stata packs scale metadata and origin above a
zero-padded square factor, distinguished from an inverse by its row count;
this needs no change to the outer saved-model version. Projection skips exact
null coordinates before normalization and recovers overflowing raw differences
on the fitted scale. Exceptional rows project normalized finite differences
first and then restore the binary ratio; use operand normalization only when
raw subtraction overflows. This avoids a large common level erasing tiny
query differences. Stata's existing
single-covariate estimator support remains; its one-column header omits rank.

## Covariance arithmetic and failures

`.covariance_matrix(x, factor, center, diagonal, allow_failure)` computes
`factor * crossprod(centered_x)`. The factor remains an estimator decision;
jackknife callers use `(n-1)/n`. Weighted moment callers pass square-root
weighted signals with `center = FALSE`, preserving their uncentered meat.

Ordinary columns retain vectorized arithmetic. Exceptional columns use exact
binary normalization before centering or squaring. High-amplitude detection
uses the common conservative gate `sqrt(2^1020/n/max(1,factor))`; positive
amplitudes or centered spreads below `sqrt(DBL_MIN)` also normalize. A shifted
observed origin protects local variation at large offsets. Restore each
cross-product using two binary exponent factors, so their product is never
formed first. Full and marginal requests use the same diagonal reduction;
marginal prediction allocates no quadratic target matrix. Symmetrization
preserves the diagonal and uses a guarded midpoint without doubling entries.

Range protection and small cross-covariance accuracy need separate checks.
Column normalization can erase an original small entry or product before
restoration. For sensitive off-diagonals involving a normalized column, retain
original-scale data and reuse the observed-origin/two-difference centering
policy. Cache each recovered column once, then compensate its separately
rounded products. Use that result only when the direct products, their sum
and final scaled value are finite and every nonzero-operand product remains
normal. Separately rounded nonzero subnormal products can cancel to a false
zero even when the normalized cross-entry remains representable; zero-only
underflow checks are insufficient. The normal-range boundary uses
[R's `double.xmin`](https://stat.ethz.ch/R-manual/R-devel/library/base/html/zMachine.html).
Otherwise retain the normalized calculation, including aggregate subnormal
recovery. Ordinary and marginal-only calls retain no original-data matrix or
second product pass; exceptional full covariance uses O(n*q) extra scratch.

This is an adaptive range repair, not an exact dot-product algorithm or a
correct-rounding guarantee. [LAPACK DLASSQ](https://netlib.org/lapack/explore-html/d8/d76/group__lassq_gae8f40b0a34771b4f2d9c863de3af7be5.html)
documents safe scaling for sums of squares; that norm contract does not certify
every small signed covariance entry. [Ogita, Rump and Oishi (2005)](https://ogilab.w.waseda.jp/ogita/math/doc/2005_OgRuOi.pdf),
Sections 2 and 5, separate overflow assumptions and multiplication underflow
from accurate accumulation. The package applies that distinction through
independent binary fixtures for erased inputs/products, centering residuals,
exact cancellation and protected overflow, including public paired covariance.

A genuinely constant column has zero variance. Nonzero variation whose final
variance is zero, nonfinite or negative is unavailable. Ordinary covariance
requests raise an actionable rescaling error; point-only fits remain possible.
Paired-ratio adapters request failure masks and retain their existing
`numerical_failure` status, preserving unaffected targets. The final value must
fit the host's representation: Stata's nonmissing finite range is smaller
than R's IEEE range. Scaling cannot recover a mathematically unrepresentable
answer. These guards change numerical computation, not inference assumptions.

## Score comparisons and contribution ranges

`R/utils-binary.R` owns exact operand decomposition into a mantissa and a
power of two. PLQR comparisons and contribution-range recovery reuse it;
geometry and covariance retain their separately documented scaling policies.

`R/utils-plqr-comparison.R`, the native `linear_threshold_comparison.hpp` and
Stata's `_rdistnn_plqr_mata.ado` implement the same adaptive comparison for
`Y <= D*theta + g`. A finite threshold alone does not ensure a correct sign:
ordinary cancellation can depend on the BLAS reduction, and finite products
can overflow a partial sum. Let `L = ncol(D)+2`, `eps = 2^-52`,
`eta = 2^-1074`, and let `A` add `abs(g)`, `abs(Y)` and the absolute rounded
products in column order. Retain the ordinary comparison only when the
threshold and margin are finite, `A < 2^1022`, and the absolute rounded margin
strictly exceeds `8*(L+1)*eps*A + 32*L*eta`. The conservative error envelope
includes product, sum, margin and gradual-underflow rounding. The range gate
supplies the no-overflow condition needed for ordinary dot-product error
analysis; see [Ogita, Rump and Oishi (2005)](https://ogilab.w.waseda.jp/ogita/math/doc/2005_OgRuOi.pdf).

Uncertain comparisons use an exact fixed-point sum of products, following the
[exact-accumulation approach of Kulisch and Snyder](https://publikationen.bibliothek.kit.edu/1000018822/1375602)
and the [adaptive-predicate principle of Shewchuk](https://people.eecs.berkeley.edu/~jrs/papers/robustr.pdf).
Separate positive and negative accumulators each use 170 base-`2^26` limbs
with bit origin `-2252`. Split significands and their products before shifting;
every integer operation remains exactly representable in binary64. Immediate
carry normalization bounds intermediates independently of the number of terms.
The storage extends below every finite binary64 product and above the sum of
`INT_MAX+2` maximal products. Comparing the highest differing limbs gives the
exact sign; equality satisfies the inclusive comparison. This needs no runtime
multiprecision dependency. The contract assumes round-to-nearest binary64,
gradual underflow and builds without unsafe floating-point reassociation.

`R/utils-score-concentration.R` shares finite-operand contribution recovery
between CATE concentration and CATE/CLATE influence plots. Zero-weight rows
are inactive before risky arithmetic. Ordinary finite weighted residuals keep
their represented values. Recover a range failure only when both residual
operands, including the already rounded paired product, are finite; retain
mantissa/exponent descriptors until the weighted contribution or dimensionless
share is formed. A finite concentration can remain available when dimensional
contributions exceed the host range. This does not infer the magnitude of an
overflowing paired product or change paired jackknife concentration formulas.

## Verification

Independent binary-power tests exercise overflowing raw sums, recoverable
subnormal aggregates, unrepresentable nonzero variance, exact constants and
large symmetric matrices. Latent orthogonal coordinates certify the
near-collinear Mahalanobis counterexample and exact neighbor choices.
Cover native/fallback dispatch, vector/scalar/batch results, singular and wide
designs, fitted covariance, candidate views and restored predictions. Shared
R-exported Stata fixtures complement analytic references; equal random seeds
are not cross-language fixtures. The numerical parity policies in notes 09
and 10 apply.
