# 06 — Learners, a bounded benchmark and numerical studies

[Index](README.md) · [Previous](05_inference_confidence_sets_and_diagnostics.md) ·
[Next: tests](07_unit_integration_and_regression_tests.md)

**Status:** the bounded generator and internal ridge reference adapter are
implemented and tested. The practical-learner, production simulation and
empirical programs below remain future work. Small execution tests are
not accepted manuscript evidence.

## 1. External nuisance fitting first

The [live roadmap](../rdistnn_roadmap.md) deliberately keeps nuisance fitting
outside current high-level wrappers. Preserve that architecture for
`clate_dnn()`: five externally fitted functions form one nuisance state.
They can be fitted separately or jointly, but both augmented responses
share the same instrument propensity.

An adapter must eventually supply Y,D,A,e,m_y0,m_y1,m_d0,m_d1 with a common
row map, clipping specification and training record. It must distinguish
literal full-sample predictions, OOB predictions and genuinely held-out
cross-fitted predictions. Outcome and treatment regressions are indexed by
instrument arm. In a one-sided model m_d0 is known zero and only four
functions need fitting.

Do not hide fitting inside predict(): it evaluates new reporting profiles
by relocalizing the original frozen signals. A future new-household
prediction task would require a different API and target definition.

## 2. Bounded latent-compliance benchmark

A proposed `clate_demo_data(n, seed, heterogeneous = FALSE)` can later
generate both observed data and explicitly named oracle nuisance columns.
It is a synthetic teaching/reference function, not a replacement for
`sim_iv_data()` or a new alias for its existing designs.

Let X be uniform on [0,1]^2, Z uniform on [0,1]^4, H uniform on [0,1] and
epsilon take values +/-1/8 equiprobably, mutually independent. With
S_j=sin(2*pi*Z_j), set

$$
e=1/2+S_4/8,\quad a=1/4+S_1/32,\quad
c=1/3+S_2/32,\quad b=a+c.
$$

Draw A conditional on Z as Bernoulli(e), independently of H,epsilon and X.
Define D(0)=1{H<=a}, D(1)=1{H<=b}, D=D(A), and

$$
g=(S_1+S_4+X_1S_2+X_2S_3)/8,\quad
\vartheta(X)=(X_1+X_2)/8,
$$
$$
Y(d)=g+\vartheta(X)d+(H-1/2)/8+\epsilon.
$$

This construction has bounded Y, monotonicity, instrument overlap and a
first stage Delta_D=1/3. The latent H induces treatment endogeneity:
Cov(D,H|V)=-{(1-e)a(1-a)+eb(1-b)}/2<0. The shared S_4 term also makes
adjustment beyond X matter for instrument exogeneity.

True regressions are m_d0=a,m_d1=b,m_y0=g+vartheta*a and
m_y1=g+vartheta*b. Conditional CLATE equals vartheta in this core model.
Its oracle constants are

$$
\Sigma_x=13/(48\sqrt{15}),\qquad
\Omega_x=39/(16\sqrt{15})\approx 0.6293597938.
$$

The unadjusted conditional Wald ratio is vartheta+3/32, supplying a
meaningful negative control. Do not put the generated latent H in the
learner adjustment set: it is unobserved in the intended experiment.

For heterogeneous=TRUE, replace vartheta by
tau(X,Z)=vartheta(X)+sin(2*pi*Z_2)/16. Then CATE is vartheta but CLATE is
vartheta+3/1024. Oracle variance must use the heterogeneous construction's
formula, not the core constant.

## 3. Generator return and verification

Return data containing observed X,Z,A,D,Y and oracle columns, truth metadata
with the target function, first-stage function and core covariance when
applicable, plus configuration and seed. Keep latent variables in an
explicit optional audit payload rather than in the default regression
column set. Do not label a single realized draw as causal evidence.

Verify formulas by deterministic integration over latent compliance
intervals and the two epsilon values. Numerical quadrature over Z can
verify covariance constants independently of a learner fit. Unit tests
check finite/binary/bounded output and deterministic seed behavior, and
preserve the caller's RNG state if the package's existing demo helpers use
that convention. Monte Carlo sample means are secondary loose checks,
not brittle replacements for analytic identities.

## 4. A designated ridge-series adapter

Start with an internal study adapter
`.clate_ridge_predictions(data, M, lambda, normalizer_n, delete = integer())`
in a future `R/sim-clate-learners.R`. It returns the five aligned prediction
vectors and a training record. It need not become a general public
automatic-learner interface.

Use q(X)=(1,sqrt(12)(X_1-.5),sqrt(12)(X_2-.5)) and J=M^4 equal adjustment
cells. The dictionary p=q tensor (sqrt(J) cell indicators) has K=3J and
supremum squared norm <=7J. Fit propensity with all observations and fit
the two outcome/treatment pairs within instrument arms:

$$
\widehat\beta_j^{[-I]}=\arg\min_\beta
\{n^{-1}\sum_{i\notin I}L_{ji}(R_{ji}-p_i^\top\beta)^2+
                                       \lambda_n\|\beta\|^2\}.
$$

Keep the ambient normalizer n, dictionary, penalty and clipping fixed for
every deletion. For a solver averaging over r retained arm observations,
the corresponding penalty is n*lambda_n/r. Handle r=0 explicitly through
the defined penalized normal equations rather than dividing by zero.
Clip e to [1/4,3/4], outcome predictions to [-1,1] and treatment predictions
to [0,1] in this designated model.

Five response fits use only three Gram matrices: full, A=0 and A=1.
Test the actual weighted normal equations and one-/two-deleted predictions,
including unoccupied cells. A Sherman–Morrison comparison is a validation
oracle only when the chosen full/deleted inverse and residual match the
exact identity; compare against explicit refits.

## 5. Rates and theorem boundary

The notes reuse the paper's weighted projection and common-Gram proof:

$$
E\rho_n^2\le C\{K^2/n+K^{-1/2}+K\lambda_n^2+P(E_n^c)\},
\quad \delta_n=O(K^2/n),
$$
$$
P(E_n^c)\le CK\exp(-cn/K).
$$

Choose M=max(1,floor(n^(1/12))), lambda=n^(-1/2) and s=floor(n^(3/4)).
Then N is order n^(1/4), expected squared prediction radius is
O(n^(-1/6)), and the product, deletion and undersmoothing payments are
n^(-1/24), n^(-5/24) and n^(-1/4).

These conservative sufficient rates are note-level verification for this
construction, not a theorem for an arbitrary adapter. The four learner
conditions still include expected retained prediction convergence,
sqrt(N) times the full-fit mean remainder tending to zero, sqrt(N)delta
tending to zero and N P(E_n^c) tending to zero, with measurable continuous
conditional versions and bounded envelopes.

Report realized M,K,s,N in every study. The floor in M makes finite
refinement very slow; the asymptotic powers are not an adequate description
of every small sample. Successful deletion-refit comparisons test the
algorithm, not uniform asymptotic learner rates.

## 6. Future practical learner adapters

A later forest adapter can use the already Suggested grf dependency, checked
only when invoked. It must document instrument-arm training, all-tree/OOB/
held-out prediction, seed offsets, common row IDs and clipping. A five-fold
benchmark fits up to 25 regressions in the unrestricted model, or 20 in a
one-sided model, unless joint learning changes the architecture.

For pooled cross-fitting, pool numerator and denominator averages before
division. Do not average fold ratios arithmetically. Freeze fitted values,
geometry and the chosen s for covariance. Cross-fitting training details
do not alter the fixed-array variance definition.

No forest adapter inherits the bounded-series theorem merely because it
returns the same prediction names. A full training convenience layer should
be separately reviewed against the roadmap after the supplied-array feature
is complete.

## 7. Study sequence and saved evidence

Start with exact-array checks, then a small oracle benchmark, then bounded
ridge full-sample/known-propensity/five-fold comparisons on paired draws.
Use the core and heterogeneous models before adding smaller fixed compliance,
shrinking compliance, reporting curvature or nuisance misspecification.
Shrinking compliance must be proportional, such as
c_n(Z)=c0_n*(1+eta*S_2), 0<=eta<1; retaining a fixed sine amplitude as c0_n
vanishes can create invalid negative probabilities.

Save components, their full 2 by 2 covariance, ratio, both variance scales,
N, max weight, radius, training identity, clipping rates, statuses and
draw/seed IDs. Separate nuisance/oracle gaps, oracle sampling error and
localization bias. A realized structural target is
sum(w_i*c(V_i)*tau(V_i))/sum(w_i*c(V_i)); it is not generally
sum(w_i*theta_C(X_i)).

Report bias, Monte Carlo SD, average estimated variance, studentized
quantiles, coverage and interval length together. Retain failed and
zero-denominator draws. Give unconditional success-and-coverage as well as
explicitly conditional summaries. Set inversion can return rays/all-real/
empty sets; those are categories to report, not failed finite endpoints
to discard. CI smoke runs assert schema and reproducibility, not coverage.
Production counts and claims need a separately approved study contract.

## 8. SIPP participation is a new analysis

The accepted SIPP code uses D for eligibility and P for participation.
In this feature's notation A=eligibility and D=participation. A future
adapter must use explicit roles and cannot mutate the accepted variable
meaning or frozen replication source.

The saved eligibility outcome signal can be the CLATE numerator only after
checking identical rows, outcome, adjustment, predictions, clipping and
geometry. Treatment-by-instrument-arm predictions and paired covariance are
new work. Verify the actual saved schema and row IDs; vector lengths alone
are inadequate.

The inspected accepted sample has 9,915 households and no observed
participation among ineligible households, according to the scientific
packet. Recheck the actual selected data before using one-sided compliance.
Raw assets, practical forests, sample-derived geometry and integer age do
not become the bounded primitive example. Eligibility exogeneity, exclusion,
timing and wealth accumulation need substantive justification.

Keep participation outputs, estimand labels and acceptance manifests in a
new study directory. Do not add them to the existing affine replication
manifest or replace eligibility tables through this package feature.
