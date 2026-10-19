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
formulas. The cutoff is shared across packages; SVD libraries can still differ
at its rounding boundary. Do not claim universal bitwise SVD equivalence.

Contexts retain scale, raw observed origin, factor, version and rank. Evaluation
must respect this state for subsets, duplicate/reordered targets and restored
fits. Legacy R `Si` contexts and Stata square inverse snapshots continue to
use their original metric. Stata packs scale metadata and origin above a
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

A genuinely constant column has zero variance. Nonzero variation whose final
variance is zero, nonfinite or negative is unavailable. Ordinary covariance
requests raise an actionable rescaling error; point-only fits remain possible.
Paired-ratio adapters request failure masks and retain their existing
`numerical_failure` status, preserving unaffected targets. The final value must
fit the host's representation: Stata's nonmissing finite range is smaller
than R's IEEE range. Scaling cannot recover a mathematically unrepresentable
answer. These guards change numerical computation, not inference assumptions.

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
