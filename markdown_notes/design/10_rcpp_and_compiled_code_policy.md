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

## Diagnostic signaling

Follow [the reporting design](05_error_messages_and_warnings.md). Numerical
kernels return facts/statuses; wrappers own user context, warning aggregation
and actionable errors. Do not format messages or signal host conditions from
OpenMP worker loops. Signal R errors/warnings on the main thread after workers
join. Keep native invariant and direct-call guards, and contain every exception
at the Stata plugin ABI with allocation-safe failure handling.

Distinguish unavailable acceleration from actual execution failure. Automatic
fallback needs a known unavailable capability or a specifically documented
equivalent recovery route. Unexpected invariants, allocation failures and
interrupts must not silently become "native unavailable". Preserve applicable
host resource codes. Do not choose recovery by matching English error text or
add a formatting wrapper to every successful inner call for cosmetic parity.
The [testing policy](09_testing_and_regression_policy.md) requires separate
capability, forced-engine, genuine-failure and cleanup controls.

The companion ABI maps allocation failures to Stata's resource code 909 and
unexpected native failures to 498, retaining specific input/transport codes.
Host bridges clean temporary matrices and files before propagating the original
failure. Only declared missing capability permits automatic portable dispatch.
R's existing capability selector likewise leaves real native exceptions visible;
do not add a generic retry around it.

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
decision. Ordinary distance kernels subtract before applying fitted inverse
scales, accumulate coordinates in order using doubles, and round each L2 squared
coordinate before addition. Unscaled L2 uses a shared exceptional-row recurrence
when the accumulated squared sum reaches `2^1023` or a nonzero coordinate's
square is subnormal or zero. Scaled L1/L2 also recover rows whose raw difference
reaches `2^1023` in magnitude or whose nonzero scaled coordinate product becomes
subnormal, including zero. The common boundary accounts for Stata's
missing-value range for both signs; waiting for IEEE overflow in R/C++ can
change ranks relative to Mata.

Exceptional scaled rows retain each coordinate product as a binary mantissa
and exponent before forming the L1 sum or L2 scaled sum of squares. This keeps
individually unrepresentable scaled coordinates available when their complete
norm is representable. For huge raw differences, this recovery uses separately
rounded binary quarters before subtraction. Restore the final binary power in
two balanced steps to avoid premature subnormal rounding. Recovery uses constant
accumulator state and rescans only exceptional rows; ordinary rows keep their existing
arithmetic. A coordinate is active when its raw difference and inverse scale
are both nonzero. With exactly one active coordinate, retain its directly
rounded scaled magnitude; squared R L2 APIs then round its square separately.
With multiple active coordinates, an unsquared recovered norm is at least the
largest directly rounded scaled magnitude. No additional lower bound is applied
to multi-coordinate squared results. Recovery does not guarantee correct rounding
for every norm. Squared R APIs retain genuine squared overflow. Keep native and
fallback implementations aligned.

High-level ranks use exact represented distance then original row index.
Per-target groups contain every observation at exactly the same selected-metric
distance, including distinct raw covariate profiles and zero-weight ranks beyond
nominal support. Average assigned rank weights over observations in each group;
never evaluate a weight rule at a mean rank. No tolerance, additional rounding, outcome, treatment or
nuisance defines the groups. Signed and unnormalized mass is preserved.
Bare supplied ranks/orders remain strict orders. Supplied distances retained
alongside a rank/order object define its target-specific groups without changing
the supplied permutation. Deletions remove the actual observation and average
the selected `n-1` table over each surviving group, dividing by the survivor
count. Scale fitting and saved prediction reuse the fitted geometry and the
selected `exact_distance_block_averaging_v1` policy.

R, C++ and Mata must implement the same block arithmetic through shared weight
and deletion helpers. Accumulate original signals with each block's actual
unshifted, shifted or survivor-averaged weights. Retain the leading and correction
states until the outside-block and within-block survivor sums are merged.
Rounded group-mean signals can erase a finite survivor after actual deletion.
Canonical contiguous blocks retain linear work in sample size; explicitly
interleaved supplied groups use direct regrouping within the selected engine.
Tests establish actual engine execution and compare independent explicit
deletions. Public calls warn that inference with persistent distance ties is not
established. The paper assumes a continuous distance distribution, under which
the previous and current policies agree almost surely; averaging does not
establish a uniform MSE improvement.
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

The shared observed-origin policy accepts actual weight summaries for deletion
calculations. Zero is used for subnormal weighted amplitude, nonzero subnormal
mass, or nonzero mass times supported amplitude below the normal range. Native
and R paths preserve original rounded products and divide by mass afterward;
zero signed mass alone retains centering. Reuse table summaries rather than
adding pure-R table scans to compiled paths. Finite-support native flattening
uses `size_t` before multiplication so valid long-vector buffers do not overflow
an intermediate `int` index.
