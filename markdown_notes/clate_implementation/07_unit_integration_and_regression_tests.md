# 07 — Unit, integration and regression-test plan

[Index](README.md) · [Previous](06_learners_bounded_model_and_studies.md) ·
[Next: documentation](08_documentation_examples_and_vignettes.md)

**Status:** test specification for future code. The named files below are
proposed additions, not tests run or delivered by this planning collection.
Use testthat edition 3 and the existing [testing policy](../design/09_testing_and_regression_policy.md).
No test should need the paper checkout, the Literature vault or downloaded
SIPP data.

## 1. Independent references before optimized comparisons

Create `tests/testthat/helper-clate-reference.R` with small, deterministic
functions for augmentation, subset enumeration, literal retained-sample
averages and a direct paired covariance. Helpers define functions only;
they do not fit models or run expensive work when loaded.

The reference weight oracle enumerates combn(seq_len(n),s), selects each
subset's closest observation using the stated tie rule and averages winners.
It does not call a package rank-weight table. The retained reference removes
one observation, preserves original relative order and repeats that
enumeration at the same s. Limit exhaustive cases to 2<=n<=8 and 1<=s<n.
Use fixed signal arrays, so results do not depend on a learner.

A separate augmentation oracle expands each term of the formula directly.
The covariance oracle constructs pseudo-values explicitly and takes centered
cross-products. Production prefix/suffix helpers, shared covariance helpers
and future constructor internals must not be called from these expectations.
Comparing two wrappers over the same faulty helper is insufficient.

## 2. Proposed suite map

| New file under tests/testthat/ | Contract and named tests |
| --- | --- |
| `test-clate-signals.R` | "augmented responses match explicit formulas"; "product remainder has both positive regression-error terms"; "one-sided and unrestricted signals agree" |
| `test-user-affine-ratio-dnn-core.R` | "hand fixture returns paired ratio covariance"; "subset enumeration matches full and retained arrays"; "zero deleted denominator does not invalidate frozen variance" |
| `test-user-affine-ratio-dnn-validation.R` | "unsupported weights fail before distance evaluation"; "covariance controls and n/s boundaries are explicit"; "optional returns are additive" |
| `test-user-clate-dnn.R` | "wrapper equals signal constructor plus ratio core"; "perfect compliance reduces to CATE"; "instrument relabeling preserves ratio uncertainty" |
| `test-user-clate-dnn-data-interface.R` | "matrix data-frame and tibble inputs agree"; "column order and single-target dimensions are preserved"; "mixed input routes are rejected" |
| `test-user-clate-dnn-validation.R` | Prediction/binary/range/identity validation and precise error messages |
| `test-user-affine-ratio-methods.R` | coef/vcov/confint/summary/predict shapes, NA targets, serialized prediction state |
| `test-user-clate-diagnostics.R` | Components, denominator SE, N, clipping, instrument support and pseudo-value concentration |
| `test-clate-engine-parity.R` | Fast R, explicit loop and available compiled results; frozen geometry and ties |
| `test-clate-demo-data.R` | Later generator formulas, bounds, truth and reproducibility |
| `test-clate-ridge-learners.R` | Later ridge objectives, retained normalizers and deletion identities |
| `test-clate-confidence-sets.R` | Optional inversion geometry, PSD validation, membership and perturbations |

Existing `test-dnn-jackknife-prefix.R`,
`test-compiled-fallback-parity.R`, `test-dnn-scalar-reference.R` and
`test-result-x0-target-matrix.R` remain regression anchors. Add focused
CLATE cases instead of renaming or duplicating those suites wholesale.
Update [tests/README.md](../../tests/README.md) when real suites are added.

## 3. Exact scalar acceptance test

The first end-to-end numeric test should have these assertions:

~~~r
fit <- affine_ratio_dnn(
  X = matrix(1:4, ncol = 1L), x0 = 0,
  numerator = c(2, 1, 4, 3),
  denominator = c(1, -1, 2, 1),
  s = 2L, return = c("deletions", "pseudo_values")
)
expect_s3_class(fit, "affine_ratio_result")
expect_equal(unname(fit$components), matrix(c(2, .5), nrow = 1L))
expect_equal(fit$theta, 4)
expect_equal(unname(fit$component_vcov[[1L]]),
             matrix(c(.5, .5, .5, .75), 2L))
expect_equal(unname(fit$theta_vcov[[1L]]), matrix(34, 1L))
expect_equal(unname(fit$theta_omega[[1L]]), matrix(612 / 7, 1L))
expect_equal(fit$se, sqrt(34))
expect_equal(fit$neff, 18 / 7)
expect_equal(as.numeric(fit$deletions[, 2L, 1L]),
             c(0, 4 / 3, 1 / 3, 1 / 3))
~~~

Use explicit small-array tolerances of about 1e-12. Also assert score
pseudo-values (6,-8,1,1), mean pseudo-values equal the two full components
and no nuisance learner is called. This fixture fails an implementation
that divides each deleted numerator by its deleted denominator.

## 4. Algebraic and transformation tests

Test the generic affine transformation numerator_new=a*numerator+b*denominator:
theta_new=a*theta+b and variance_new=a^2*variance. Include a<0 and a=0;
the latter produces a degenerate score variance with an explicit status.
Test common scaling of both signal columns: the ratio and ratio variance
are unchanged for representable nonzero scales.

At the IV layer, instrument relabeling transforms A to 1-A, e to 1-e and
swaps all four arm means. Both signals negate and the ratio/variance stay
the same. This is an observed-functional identity; the reversed orientation
need not retain the original positive-monotonicity interpretation.

Adding a constant to Y and both outcome means leaves phi_y unchanged.
Scaling them scales theta and its SE appropriately. Perfect compliance
D=A,m_d0=0,m_d1=1 makes phi_d=1 and reduces the ratio, covariance and
normal-reference interval to `cate_dnn()` with the same frozen AIPW signal,
geometry and s. Expect singular component covariance in that case, not an
error requiring a positive first-stage variance.

Enumerate finite candidate/true propensities and regression means to
verify the exact nuisance-bias identity. Include known e, correct response
regressions, and a case where both are wrong and the remainder is nonzero.
Do not infer causal identification from those algebra tests.

## 5. Input and control failure matrix

Exercise every relevant branch through exported calls:

| Family | Concrete cases |
| --- | --- |
| Interface selection | data with X; data with direct Y; numerator_col without denominator_col; both X0 and x0; neither target |
| Columns | Absent named column; duplicate data column names rejected; reordered or repeated reporting selectors preserve the existing shared extractor's semantics |
| Shapes | Transposed X0; one-column vector supplied as a matrix; wrong signal length; zero rows; m=1 and k_x=1 |
| Values | Separate NA, NaN, Inf and -Inf in each input role; factor A; A=-1/1; D=2 |
| Predictions | e outside [0,1]; finite e at endpoints; m_d outside [0,1]; missing m_d0 unrestricted; conflicting nonzero m_d0 one-sided |
| One-sided model | D>A error; eligible nonparticipant with negative phi_d retained |
| Scalars | clip=0,.5,NA or vector; s=0/noninteger/>n; k=2; method="asymp"; invalid denominator_tol |
| Identity | Duplicate/NA/empty row IDs; numeric non-finite IDs; valid character IDs differing from input positions |
| Returns | Unknown option warns once; omitted optional arrays remain absent; covariance fields NULL without active covariance |

Use meaningful error substrings naming the argument and expected relation,
not full punctuation snapshots. Use a counted distance object or mock only
where it proves early validation precedes heavy work. Avoid tests pinning
incidental private call order.

Boundary cases include n=1,s=1,se=FALSE; n=2,s=1 with covariance;
s=n with covariance disabled; s=n rejected with active covariance; and
`se=TRUE,jackknife="none"` returning NULL covariance. Test all controls
through both high-level wrappers where their routes differ.

## 6. Denominator and numerical statuses

Construct separate targets with ordinary positive, negative, exactly zero
and nearly canceling denominators. Assert no target is dropped, original
target IDs survive and finite negative ratios remain computed.
At exact zero, components and component_vcov remain while theta/SE/ratio
covariance are NA. The near-zero tolerance changes only the flag.

Use a constant-proportional numerator to produce zero score variance and
assert se=0, zero_score_variance and NA Wald bounds. Use a constant
denominator with a varying numerator to verify ordinary nondegenerate
inference with singular paired covariance.

Test currency-scale signals and correlated columns with near cancellation.
The direct score-norm variance and the paired quadratic formula must agree
where both are numerically resolvable. Extreme finite inputs that overflow
derived components or ratios produce the documented numerical statuses;
non-finite supplied inputs are rejected before computation. Numerical
failure at one profile leaves other profiles intact.

## 7. Engine parity and geometry

For n in a small set including 2,4,17 and m in 1,3, test s=1,s=n-1 and an
interior s. Compare:

1. independent subset/retained reference where enumeration is feasible;
2. explicit rank-shift loop;
3. fast prefix/suffix with compiled algorithms disabled;
4. fast prefix/suffix with compiled symbols available.

Compare full components, all deletion arrays, paired covariance, ratio,
SE, statuses and shapes. Include identical distances, repeated X,
targets equal to observed X and input-row permutations. Under distinct
distances a common row permutation preserves results; under ties compare
the documented order, not a false permutation-invariance requirement.

Use an X-dependent distance setup fixture to confirm geometry is resolved
once and retained ranks are shifted from the frozen full ranks. Do not
accidentally compare a frozen rank deletion to a reestimated metric.

Restore `options(RdistNN.use_compiled)` with on.exit(). Skip only the
compiled-specific comparison if the symbols are unavailable; R/reference
tests must still run. New C++ is unnecessary initially. Any later new
compiled path needs the same tests on the supported Linux/Windows/macOS
architectures.

## 8. Reporting, prediction and diagnostics

Assert coef names and target order; vcov's single-target matrix and
multi-target list; valid target selection; all confidence levels;
NA interval rows; summary component labels and status counts.
Print methods return the object invisibly and do not mutate it.

Predicting at the original x0 must reproduce the original point estimate
and covariance. New target prediction equals an explicit fresh generic
engine call on the same frozen arrays, distance and s. It must not refit,
retune at another n or alter propensity clipping. Test matrix newdata and
data-frame newdata using retained x_cols, conflicting target arguments,
missing model state and saveRDS/readRDS.

For diagnostics, reconstruct N, first-stage SE and largest pseudo-value
shares from explicit arrays. Generic fixed-signal results have no invented
instrument/propensity diagnostics. CLATE objects expose those quantities
only from the retained observed/prediction state.

## 9. Later confidence-set and learner tests

The internal optional `.ratio_quadratic_set(a,b,c,tol)` needs all
quadratic geometries. Examples include (1,0,-1), (1,-2,1), (1,0,1),
(-1,0,1), (-1,0,0), (-1,0,-1), (0,2,-4), (0,-2,4),
(0,0,-1) and (0,0,1). These are algebra fixtures; some cannot arise from
coherent nonzero-denominator paired inputs.

At the public level test PSD and symmetry validation, zero-denominator
degeneracies, fitted-ratio membership, scale transformations, stable roots
and near-boundary status. Verify interval membership by substitution in
the original inequality at points away from roots. Infinite endpoints and
empty sets must survive serialization and printing.

For the bounded generator, test structural moments and true contrasts.
For the ridge adapter, explicitly refit every one-/two-deleted small
sample and compare normal equations, clipping and the fixed-normalizer
rank-one identity. Include zero retained arm observations and unoccupied
cells. These are finite algorithm tests, not proof of asymptotic rates.

## 10. Running and accepting the suite

Focused examples, after implementation:

~~~bash
Rscript -e 'devtools::test(filter = "^(clate-signals|user-affine-ratio|user-clate|clate-engine-parity)")'
Rscript -e 'devtools::test(filter = "^(clate-demo-data|clate-ridge-learners)$")'
Rscript -e 'devtools::test(filter = "^clate-confidence-sets$")'
~~~

Run the full package suite when shared validation/weight/deletion code
changes and before feature release. A green focused suite does not certify
unrelated existing estimators. Reuse successful verification while the
content tree is unchanged.

No CI test should require a coverage estimate near .95 from a handful of
Monte Carlo draws, an expensive forest fit, internet access or exact timing.
Tests prove formulas, implementation and documented behavior. Production
coverage, learner validity and causal assumptions have separate evidence
owners.
