# Rcpp and Compiled Code Policy

This note standardizes how compiled code should be designed, integrated, and
tested in `RdistNN`.

## 1. General rule

Compiled code exists to improve performance, not to create a second package
semantics.

The authoritative behavior of the package is the documented R-level contract.
Compiled paths must implement that contract, not redefine it.

## 2. Thin-wrapper principle

Rcpp entry points should be as thin as reasonably possible.

Prefer this division of labor:

### R layer
- user-facing validation;
- interface resolution;
- argument normalization;
- friendly errors and warnings;
- fallback dispatch;
- object assembly.

### C++ layer
- performance-critical inner computation;
- simple, explicit numeric contracts;
- no hidden API policy.

## 3. Canonical input contract

Compiled code should receive already validated canonical inputs whenever possible.

Preferred inputs to compiled code:
- numeric matrices;
- numeric vectors;
- validated integer scalars;
- explicit shape assumptions.

Do not push user-facing ambiguity into C++ if it can be resolved cleanly in R.

## 4. Generated files

Generated files are not hand-edited.

This includes:
- `RcppExports.cpp`
- `RcppExports.R`

Rules:
- modify source code or attributes;
- regenerate exports through tooling;
- commit generated updates in the same change that requires them.

## 5. Semantics parity with R fallbacks

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
retain stable original-row tie breaking. Scale fitting and saved prediction
must reuse the same fitted geometry contract.

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

## 6. Error-surface policy

User-facing errors should usually be raised from R, not only from C++.

Prefer:
- validating early in R;
- entering C++ only after the contract is already clear.

If C++ must detect an internal impossible state, surface it clearly and avoid
cryptic crashes or undefined behavior.

## 7. Indexing and shape discipline

Any nontrivial compiled code should document:

- indexing convention;
- row/column interpretation;
- shape assumptions;
- scalar parameter constraints.

Comments near the wrapper boundary should make the data contract easy to audit.

## 8. Type discipline

Be explicit about type assumptions.

Rules:
- do not rely on silent narrowing where it can matter;
- use validated integer scalars at the R boundary;
- document any assumptions about sizes, ordering, or memory layout that affect
  correctness.

## 9. Missingness and finiteness

Compiled code should not be the first place where basic user-input pathologies
are discovered.

Rules:
- NA / non-finite checks should usually happen in R;
- if compiled code still checks them, it should do so consistently with the R
  contract;
- fallback and compiled paths should not diverge on malformed-input handling.

## 10. Testing requirement

Any new compiled path must come with tests that compare it to the reference path.

The supported CI matrix must execute the compiled/fallback parity suite on
Linux x86-64, Windows x86-64, macOS arm64, and macOS x86-64.

At minimum, test:
- one ordinary valid case;
- one edge case;
- shape agreement;
- exact discrete decisions and numerical agreement at justified tolerances;
- behavior when the fast path is unavailable, if applicable.

## 11. Performance changes

A performance-motivated change is not justified if it makes semantics murkier.

Prefer:
- clarity first;
- then speed;
- then micro-optimization.

Benchmarking is useful, but not at the expense of maintainability.

## 12. Review checklist

Before merging compiled-code changes, ask:

- Does the compiled path implement the same contract as the R path?
- Is the R wrapper still the main validation boundary?
- Were generated files updated correctly?
- Are indexing and shape assumptions documented?
- Is there a fallback-equivalence test?
- Does this improve speed without making maintenance harder?
