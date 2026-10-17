# 01 — Estimator, scientific scope and source map

[Index](README.md) · [Next: interfaces](02_public_interfaces_and_results.md)

**Status:** future feature specification, grounded in the inspected snapshots
listed in the index. Existing functions below are present at that baseline;
new names are proposals. Re-read the exact sources before implementation.

## 1. The question the feature answers

Let A be a binary instrument, D binary treatment, Y an outcome and V = (X,Z)
the full set of pre-instrument adjustment covariates. X contains the few
coordinates at which the user reports heterogeneity. Z remains an adjustment
variable; it is not reused as the instrument's name.

Under consistency, exclusion, conditional instrument independence given V,
monotonicity, instrument overlap and positive compliance at x, the target is

$$
\theta_C(x)=E[Y(1)-Y(0)\mid D(1)>D(0),X=x]
           =\Delta_Y(x)/\Delta_D(x).
$$

For R = Y,D, define m_Ra(v) = E[R | A=a,V=v] and e(v) = P(A=1 | V=v).
The standardized contrasts are

$$
\Delta_R(x)=E[m_{R1}(V)-m_{R0}(V)\mid X=x].
$$

Both expectations integrate over the same law of Z given X=x. A raw
instrument-arm difference conditional only on X can average over different
adjustment distributions. An estimator of a full-V ratio followed by an
unweighted average can also answer a different question.

The computational object remains a standardized observed-data Wald functional
without a causal IV interpretation. Software validates its arrays; it cannot
verify exclusion, monotonicity or conditional randomization from those arrays.
Help and output must distinguish the functional from its causal interpretation.

## 2. Exact estimating equation

For a candidate nuisance state define

$$
\phi_R=m_{R1}-m_{R0}
  +\frac{A}{e}(R-m_{R1})
  -\frac{1-A}{1-e}(R-m_{R0}).
$$

The paper's score is psi(W;theta,eta) = theta phi_D - phi_Y. It has
B = phi_D, b = -phi_Y and point Jacobian Gamma_x = Delta_D(x). Thus

$$
\widehat\Delta_R=\sum_iw_i(x)\widehat\phi_{R,i},\qquad
\widehat\theta=\widehat\Delta_Y/\widehat\Delta_D.
$$

This is already an affine coefficient-average inverse. No new general ratio
theorem is needed for its algebra. The first stage gives the varying
coefficient an economic interpretation, and the existing varying-coefficient
covariance argument retains movement of the fitted target.

Do not compute an observation-level phi_Y/phi_D array. Individual augmented
first-stage signals can be negative or zero even with monotonicity and
positive population compliance. Do not clip those signals to [0,1], floor the
aggregate denominator or discard negative estimated first stages.

## 3. Inferential contract

Use one deterministic exact-DNN subsample order s and k=1. At distance rank r,

$$
w_{n,s,r}=\binom{n-r}{s-1}/\binom{n}{s},\quad
N_x=(\sum_iw_i^2)^{-1}\sim2n/s.
$$

The first complete implementation restricts the user interface to this exact
rule. Other weight rules remain available elsewhere in RdistNN. Exposing
asymptotic DNN, k>1, signed TDNN or custom rank weights here would need its
own scope decision and tests.

The affine theorem requires local design and oracle regularity, fixed
nonzero Gamma_x, point-target undersmoothing, bounded score coefficients and
the retained-fit prediction/deletion contract. In particular, orthogonal
mean-bias control and prediction convergence do not alone justify full-sample
nuisance reuse. The theorem's one- and two-deleted nuisance fits are proof
devices. The implemented standard error holds all nuisance predictions and
the target fixed and only deletes rows from the local averaging operation.

The limit variance is

$$
\Omega_x=\operatorname{Var}(\theta_C(x)\phi_D-\phi_Y\mid X=x)/
                 \Delta_D(x)^2.
$$

Estimator variance is Omega_x/N_x. The paired finite-array implementation in
Note 04 is the specified feasible calculation. A covariance of deleted ratios,
a retrained-pipeline jackknife and a weighted score-cross-product meat are
different calculations.

## 4. Existing engineering owners

| Existing source | Reusable behavior | Required care |
| --- | --- | --- |
| [user-dnn-core.R](../../R/user-dnn-core.R) | Column extraction, target resolution, weights, tuning resolution, matrix-signal deletions and covariance | Point estimates from `.dnn_estimate_from_signal()` are q by m; convert once to the public m by 2 convention |
| [api-weights.R](../../R/api-weights.R) | Validated weight wrappers and optional ranks/distances | W and ranks are n by m; columns correspond to targets |
| [wrank-dnn_exact.R](../../R/wrank-dnn_exact.R) | Exact without-replacement weights | s is subsample competition size, not the number of nonzero weights |
| [jackknife-rank-full.R](../../R/jackknife-rank-full.R) | R/C++ full-rank matrix-signal deletions | Signal rows stay in original input order; ranks map each row to a distance rank |
| [utils-firststage.R](../../R/utils-firststage.R) | Existing propensity validation/clipping and AIPW formulas | Propensity validation rejects values outside [0,1] before clipping; avoid accidental double clipping |
| [utils-validate.R](../../R/utils-validate.R) | Specific shape, length, binary and finiteness errors | Resolve invalid inputs before weights or C++ |
| [user-dnn-cate-methods.R](../../R/user-dnn-cate-methods.R) | Reporting and prediction conventions | Reuse conventions without inheriting CATE-specific fields or assumptions |
| [user-dnn-diagnostics.R](../../R/user-dnn-diagnostics.R) | Weight concentration, radius and metadata patterns | CLATE needs instrument/first-stage labels rather than treatment overlap labels copied blindly |
| [utils-compiled.R](../../R/utils-compiled.R) | Global compiled/fallback switch | Respect `options(RdistNN.use_compiled = FALSE)` |

The paper's `sipp_401k_frozen_rank()` is a useful independent sorted-array
reference. It lives outside the package and must not become a runtime import.
The installed replication copy is hash-protected; new reusable functionality
belongs in R/, with tests under tests/testthat/.

## 5. Why existing local IV-GMM is a separate workflow

[gmm_model_dr_iv()](../../R/api-gmm-models.R) evaluates residualized instrument
moments of the form
(A-e(V)){Y-m_Y(V)-theta(D-m_D(V))}. With heterogeneous stratum effects it
introduces e(V){1-e(V)} weighting in addition to compliance weighting.
Its companion [GMM wrapper](../../R/user-dnn-gmm.R) also has its own
optimization and covariance route.

For two equally likely strata with (c,tau,e) equal to (1/4,0,1/2) and
(3/4,2,1/4), CLATE is 3/2, the average stratum LATE is 1, and the residualized
IV ratio is 18/13. This fixture belongs in tests and the conceptual vignette.
Do not rename the current IV-GMM API or describe its existing outputs as the
new estimator. A custom augmented GMM moment does not automatically select
the frozen rank-jackknife variance.

## 6. Delivery boundaries

| Stage | Delivered feature | What its verification establishes |
| --- | --- | --- |
| Fixed arrays | `affine_ratio_dnn()`, result methods and paired covariance | Defined local ratio and covariance computed correctly |
| Prediction arrays | `iv_augmented_signals()`, `clate_dnn()`, dual interfaces | Augmented signals and ratio agree with explicit formulas |
| User documentation | Help, examples, two vignettes and site entries | A user can reconstruct and interpret every reported component |
| Bounded reference | Generator and designated ridge adapter | A reproducible example corresponding to the notes' sufficient construction |
| Optional extensions | Confidence-set inversion, learner adapters, grid views | Their separately stated algorithms and contracts |

The first three rows form a complete initial user feature. Production Monte
Carlo, practical forest validity and a new SIPP participation analysis are
separate deliverables. Preserve the accepted eligibility analysis and the
frozen affine replication bundle, and keep numerical checks separate from
proof of learner rates.

## 7. Acceptance questions

An implementer should be able to point to one source owner for each signal,
one common rank/weight operation, one paired covariance calculation and one
result-shape contract. Every explanatory statement must identify whether it
concerns deterministic arrays, statistical sufficient conditions or causal
assumptions. The index's scientific references and these formulas supply a
portable specification; no sibling checkout becomes an installation
dependency.
