# Rcpp and compiled code policy

Compiled code accelerates the documented R contract. R wrappers own public
validation, interface/argument resolution, diagnostics, fallback selection and
object assembly; native routines own explicit numeric inner computations.
Keep wrappers thin and pass validated numeric matrices/vectors and integer
controls. Native boundaries still protect their own invariants and reject
malformed inputs consistently rather than crashing.

Document indexing, matrix axes, shapes, parameter ranges, sizes and memory
layout. Avoid silent narrowing. Edit source/attributes and regenerate
`RcppExports.cpp`/`RcppExports.R`; never hand-edit them, and commit required
generated changes together. See [CONTRIBUTING.md](../../CONTRIBUTING.md) for tools.

## Semantics and arithmetic parity

When an R fallback and a compiled path both exist, they must agree on:

- valid-input semantics;
- output shapes;
- return conventions;
- edge-case behavior where documented;
- ranks, tie resolution, support, computational statuses and missing-value masks;
- numerical values within the quantity-specific tolerances defined by tests.

For methods shared with Stata, apply this same contract across R, C++ and Mata
after mapping documented indexing, matrix layout and unavailable-value
representations. Backend selection and automatic fallback are performance and
availability choices; they must preserve the selected estimator's semantics.

### Arithmetic and reproducibility

Specify operation order and precision where rounding can change a discrete
decision. Distance kernels subtract before applying fitted inverse scales,
accumulate coordinates in order using doubles, and round each L2 squared
coordinate before addition. L2 uses a shared exceptional-row recurrence when
the accumulated squared sum reaches `2^1023`, a nonzero coordinate's square is
subnormal or zero, or a scaled raw difference reaches `2^1023`. This common
boundary accounts for Stata's missing-value range; waiting for IEEE overflow
in R/C++ can change ranks relative to Mata. For a huge raw difference and an
inverse scale below one, both L1 and L2 compute the difference of binary-halved operands and
multiply by twice the inverse scale. Recover the norm with scaled sum of
squares and constant accumulator state; rescan only exceptional rows. Squared
R APIs still return squared results, including genuine squared overflow.
Keep native and fallback implementations aligned;
break exact distance ties lexicographically by original raw covariates, in
their supplied column order, then by row position for identical covariate rows.
Distance-only sorting primitives without covariates retain stable index ties;
explicit supplied ranks/orders remain authoritative. Outcomes, treatment,
nuisances and row identifiers never define covariate-aware tie order. No tie
averaging or tolerance grouping changes rank weights or ESS. Scale fitting and saved prediction
must reuse the same fitted geometry contract.
Current saved fits reuse their retained geometry and cached ranks. Older
development fit schemas do not require compatibility readers.

Fitted scale helpers reject nonfinite raw scales before epsilon clamping.
SD calculations use common conservative range gates before intermediate
overflow/underflow. Exceptional columns use exact binary normalization, an
observed-value shift and compensated mean/squared-residual sums. This preserves
local variation around large offsets. Ordinary SD also uses the shared
observed-origin/exact-subtraction policy before vectorized centering.
The Stata companion uses the same recovery convention and R-compatible guarded
quantile interpolation. Sensitive weighted/prefix sums use shared compensated
arithmetic; [note 19](19_weighted_accumulation.md) defines ownership, correction
state, workspace and cross-engine verification.

Align engines toward stable mathematical reference calculations. Agreement
between implementations is insufficient if they share an accuracy defect.
Prefer algebraically equivalent formulas that avoid cancellation, including
`log1p`/`expm1` for small differences and direct tail probabilities instead of
subtracting CDF values near one. Use scaled sums of squares where intermediate
squaring would overflow or underflow. Consider compensated sums in existing
accumulation loops when cancellation is consequential. Assess constant runtime
and allocation costs, especially before replacing vectorized R/Mata operations
with scalar loops. Apply a stability change consistently to affected engines
and preserve the documented statistical method.

Apply the same discipline to rank-weight normalization, support-aware centering
and deletion/covariance accumulation. Reuse shared helpers within each package
so estimator-specific paths cannot acquire their own arithmetic conventions.
Compiler contraction, wider accumulators, reassociation and parallel reductions
must not silently change an ordering or computational status. Protect sensitive
operations explicitly and test compiler/architecture-sensitive counterexamples.

Require exact discrete behavior and pursue exact numerical agreement where
operations can be controlled. General bit-for-bit identity across all systems
is not a package-wide promise: platform math libraries and linear algebra may
produce small differences. Tests must use justified numeric tolerances while
still checking discrete behavior exactly. A discrepancy that affects a decision
requires investigation rather than a broader tolerance. See
[the testing policy](09_testing_and_regression_policy.md) for acceptance coverage.

The package-wide diagnostic switch is
`options(RdistNN.use_compiled = FALSE)`. Internal compiled dispatch must
respect that switch and use the corresponding R implementation. This is an
algorithm-selection control, not a no-compiler installation mode: source
installation still compiles the package's registered C++ code.

## Verification and performance

Every OpenMP parallel region uses `rustat::thread_limit()` in its `num_threads`
clause. The default cap is two; `RdistNN.num_threads` accepts a positive integer
and remains bounded by runtime limits. Resolve the option on the calling R
thread. Do not change the global OpenMP runtime defaults. Serial builds retain
the same calculations. Tests establish actual native teams and unchanged results.

Every new native path needs reference comparisons covering ordinary/edge
inputs, shapes, exact discrete decisions, justified continuous tolerances and
unavailable-fast-path behavior. Native/fallback suites run on Linux x86-64,
Windows x86-64, macOS arm64 and macOS x86-64. Force and establish actual engines
as required by the [testing policy](09_testing_and_regression_policy.md).

Benchmark representative runtime/allocation costs. A faster implementation must
preserve the statistical method and readable ownership; performance does not
justify hidden API policy or weaker validation.
