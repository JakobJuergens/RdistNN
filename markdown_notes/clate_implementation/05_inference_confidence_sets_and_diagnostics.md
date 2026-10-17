# 05 — Inference, confidence sets and diagnostics

[Index](README.md) · [Previous](04_localization_and_deletion_engine.md) ·
[Next: learners and studies](06_learners_bounded_model_and_studies.md)

**Status:** implemented reporting, diagnostics and separate ratio confidence-set
inversion. The original staged design below is retained; implementation does
not establish uniform weak-instrument validity.

## 1. Three covariance objects with distinct scales

Let Xi_x be the conditional covariance of (phi_Y,phi_D) at x and
Gamma_x=Delta_D(x). At the identifying target,

$$
\Sigma_x=\Xi_{YY}+\theta_C^2\Xi_{DD}-2\theta_C\Xi_{YD},\qquad
\Omega_x=\Sigma_x/\Gamma_x^2.
$$

The conditional score variance Sigma and local-limit covariance Omega are
O(1) objects in the fixed-law theorem. The package's component_vcov and
theta_vcov are instead on the estimator scale. Their normalizations are
approximately Xi_x/N_x and Omega_x/N_x. theta_omega is N_x times
theta_vcov; se is sqrt(theta_vcov[[j]][1,1]).

Neither the numerator SE divided by the first stage nor a calculation that
sets V_YD to zero reproduces the specified covariance. The generic scalar
theorem requires positive score variance, not positive definiteness of the
entire two-signal covariance. Perfect compliance gives a deterministic
first-stage signal and singular paired covariance, while CATE inference
can remain nondegenerate.

## 2. Pointwise Wald output

For finite theta_hat and strictly positive finite V_theta, report
theta_hat +/- qnorm((1-level)/2,lower.tail=FALSE)*sqrt(V_theta).
The upper-tail form avoids rounding to probability one at extreme valid
confidence levels. Validate 0<level<1.
Keep every requested target in the result and confint matrix.

At exact zero denominator, retain components/component_vcov but set the
ratio and its variance/SE/interval to NA. At a finite negative denominator,
retain the computed ratio and covariance. Print a negative-denominator flag;
do not take abs(denominator) before solving. Standard errors use its square
through the sandwich.

A zero score norm yields a computational SE of zero. Record
zero_score_variance and set wald_available=FALSE; return NA interval bounds
for that target. A point-mass numeric result does not satisfy a
strictly-positive-variance normal-reference theorem. If all covariance was
disabled, confint errors with an instruction to refit using se=TRUE and an
active jackknife.

The near-zero flag in Note 02 leaves representable outputs intact.
Its scale is numerical. First-stage standardized magnitude is a separate
diagnostic, and no universal F-statistic cutoff is supplied. Runtime code
does not silently filter targets or select a broader neighborhood because
their intervals are long.

## 3. A diagnostic contract that exposes the mechanism

Add `diagnostics.affine_ratio_result(object, ...)` and
`diagnostics.clate_result(object, ...)` in
`R/user-dnn-clate-diagnostics.R`. Reuse the existing diagnostics class and
common weight calculations after inspecting their required payload.
If a common class assumes CATE-only fields, add an estimator-specific
payload rather than manufacture a tau field or fake treatment array.

Return a target-indexed data frame with:

- numerator, denominator, both component SEs, their covariance and correlation;
- theta, ratio SE, neff, max weight and available physical-radius summaries;
- denominator sign, numerical cancellation flag, denominator standardized
  magnitude when its SE is positive, and all computational statuses;
- largest and largest-five centered score-pseudo-value-square shares when
  covariance arrays are available or can be reconstructed from stored state;
- instrument-arm weight mass and propensity clipping mass for CLATE
  objects with observed A and propensity payload;
- component roles, training declaration, method, s and frozen-geometry rule.

For zero component variance, correlation or standardized magnitude is NA
with a status; do not report an infinite first-stage statistic as ordinary
evidence. A nonzero deterministic denominator can be explicitly described
as having zero estimated first-stage variation.

Pseudo-value concentration attributes variation to deletion indices,
including changes in retained ranks. It is not the same as a household's
raw weighted score contribution. Name and document those diagnostics
separately. Effective information is a measure of rank-weight concentration,
not an observed or estimated count of compliers.

Predictive state should retain frozen A and e only when available from the
CLATE path. Generic fixed-array fits can still provide geometric and paired
covariance diagnostics; IV-specific fields are unavailable rather than
invented from the denominator signal.

## 4. Optional fixed-candidate confidence-set inversion

Propose a separate lower-level export for a later stage:

~~~r
ratio_confidence_set(components, component_vcov,
                     level = 0.95, tol = 64 * .Machine$double.eps)
~~~

Here components is the same m by 2 matrix and component_vcov the same list
of m 2 by 2 matrices defined in Note 02. The function returns one set record
per target with type, intervals, quadratic coefficients, discriminant,
numerical_status and level. intervals is an r by 2 numeric matrix with
inclusive endpoints and permitted infinite endpoints; r=0 represents empty.
Document closed finite endpoints and open infinity notation in printing.
Keep this function separate from the initial Wald-only confint method.

For candidate theta, invert

$$
(\widehat\Delta_Y-\theta\widehat\Delta_D)^2
\le z^2(V_{YY}-2\theta V_{YD}+\theta^2V_{DD}).
$$

This is a theta^2+b theta+c<=0 with

$$
a=\widehat\Delta_D^2-z^2V_{DD},\quad
b=-2\widehat\Delta_Y\widehat\Delta_D+2z^2V_{YD},\quad
c=\widehat\Delta_Y^2-z^2V_{YY}.
$$

| Geometry | Set |
| --- | --- |
| a>0, discriminant>0 | Closed bounded interval |
| a>0, discriminant=0 | Singleton |
| a>0, discriminant<0 | Empty |
| a<0, discriminant>0 | Two closed rays |
| a<0, discriminant<=0 | All real numbers |
| a=0,b!=0 | One half-line |
| a=b=0,c<=0 | All real numbers |
| a=b=0,c>0 | Empty |

With coherent positive semidefinite covariance and a nonzero fitted
denominator, the fitted ratio satisfies the inequality. Empty output then
signals inconsistent inputs or a numerical problem. At exact zero
denominator, empty output can be mathematically legitimate.

## 5. Numerical contract for inversion

Validate dimensions, symmetry, finiteness and positive semidefiniteness of
each supplied component covariance using a scale-aware tolerance. Do not
repair materially indefinite inputs. Preserve exact zeros in deterministic
fixtures. Normalize coefficients by a common positive scale to reduce
overflow, retaining their original values if representable.

For separated real roots use the stable quadratic construction
q=-.5*(b+sign_choice(b)*sqrt(discriminant)), with sign_choice(0)=1, then
q/a and c/q, sorting endpoints. Handle q=0 separately. Use the direct
linear branch when a is exactly zero. A relative cancellation bound on
b^2 and 4ac identifies near-zero discriminants.

Tolerance can make a nearly degenerate finite quadratic numerically
ambiguous. Record near_boundary or unresolved status rather than claiming
an exact geometry obtained by silently rounding meaningful coefficients
away. Tests of exact geometric cases use exactly representable fixtures;
perturbation tests additionally check membership at well-separated
candidate points. Rescaling Y should rescale the confidence-set endpoints
with the corresponding sign/order convention.

Zero variance at a candidate is an algebraic branch of the inequality,
not an automatically justified standardized test. Document that
distinction in the returned numerical record and the conceptual article.

## 6. Statistical scope of confidence-set coverage

At a fixed candidate theta, define H_theta=phi_Y-theta phi_D. The auxiliary
constant-coefficient score a-H_theta fits the affine theorem for its
conditional mean if that response's design, learner, bias and
positive-variance assumptions hold. At a fixed true theta_0, inclusion
probability of theta_0 is the acceptance probability of that one test;
simultaneous validity over all candidate theta is unnecessary for that
pointwise statement.

Uniform coverage over increasingly weak first stages or changing
data-generating processes requires a new sequence-specific analysis.
The implemented quadratic geometry does not prove it. A fixed-law
ratio theorem with nonzero Gamma is also not such a result.
Title the software output a ratio confidence set, and qualify its
statistical justification precisely.

## 7. Future plots and selection

An initial useful display is a three-panel components/ratio plot with
first-stage uncertainty and per-target statuses. Reuse plot-data conventions
and preserve finite negative estimates and NA targets. An unbounded
confidence set needs arrows/rays rather than a truncated finite interval
drawn as if bounded. Every plot should have an accessible tabular equivalent.

Implementing plotting or grid selection is not required for the initial
calculation. A future sensitivity grid must reuse both signal columns and
ranks, report all prespecified candidates and keep the default explicit.
Favorable-denominator or favorable-sign selection is not an inferential
policy supplied by these notes.
