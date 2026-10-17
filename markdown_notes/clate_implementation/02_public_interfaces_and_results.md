# 02 — Public interfaces, shapes and result methods

[Index](README.md) · [Previous](01_estimator_scope_and_source_map.md) ·
[Next: signals](03_augmented_signals_and_validation.md)

**Status:** recommended signatures for future implementation. None is an
existing export. Names follow the package's estimator suffix convention:
`affine_ratio_dnn()` refines the paper notes' provisional
`dnn_affine_ratio()` spelling. Settle this contract before writing code.

## 1. Three layers with distinct responsibilities

~~~r
affine_ratio_dnn(
  data = NULL, x_cols = NULL,
  numerator_col = NULL, denominator_col = NULL,
  X = NULL, X0 = NULL, x0 = NULL,
  numerator = NULL, denominator = NULL,
  distance = "l2", method = "exact", k = 1L,
  s = NULL, tuning = NULL,
  se = TRUE, jackknife = c("fast", "loo", "none"),
  denominator_tol = sqrt(.Machine$double.eps),
  row_ids = NULL, ctx = NULL, return = c("W")
)

iv_augmented_signals(
  Y, D, A, e, m_y0, m_y1, m_d0 = NULL, m_d1,
  clip = 1e-8, one_sided = FALSE, row_ids = NULL,
  training = "supplied"
)

clate_dnn(
  data = NULL, x_cols = NULL,
  y_col = NULL, d_col = NULL, a_col = NULL, e_col = NULL,
  m_y0_col = NULL, m_y1_col = NULL,
  m_d0_col = NULL, m_d1_col = NULL,
  X = NULL, X0 = NULL, x0 = NULL,
  Y = NULL, D = NULL, A = NULL, e = NULL,
  m_y0 = NULL, m_y1 = NULL, m_d0 = NULL, m_d1 = NULL,
  clip = 1e-8, one_sided = FALSE,
  distance = "l2", method = "exact", k = 1L,
  s = NULL, tuning = NULL,
  se = TRUE, jackknife = c("fast", "loo", "none"),
  denominator_tol = sqrt(.Machine$double.eps),
  row_ids = NULL, training = "supplied", ctx = NULL,
  return = c("W")
)
~~~

The low-level constructor takes arrays only. The high-level wrappers support
matrix/vector and data-frame paths, including tibbles. In data-frame mode,
sample variables come from the named columns; do not allow simultaneous
direct sample vectors or X. X0/x0 still specify reporting targets in either
mode. Reject incomplete or mixed routes before extracting columns.

The generic fixed-array estimator has two ordinary numeric vectors; it does
not interpret their causal roles. `clate_dnn()` constructs signals once,
calls that engine and adds IV-specific metadata. It never reuses an
observation-level ratio as a CATE pseudo-outcome. An `iv_signals` object can
be passed through its documented `phi_y` and `phi_d` fields to the generic
estimator; an additional polymorphic `signals=` argument is unnecessary.

## 2. Restricted controls and default resolution

Keep method and k in the grammar but accept only `method = "exact"` and
`k = 1L` in this feature's first version. Reject unsupported values early.
Exact weights are normalized; no public normalize switch is needed. Common
core calls pass `normalize = TRUE`.

Use `.dnn_resolve_s()` for either explicit s or a tuning specification.
If neither is supplied, it currently uses
`tuning_dnn_paper()`, whose default multiplier is 0.4. Record the resolved
integer and tuning payload. A default is a reproducible rule, not an
automatic undersmoothing certificate for arbitrary smoothness or learners.
Reject s>n explicitly even if the shared resolver changes. s=n is allowed
only for point estimation with covariance disabled. Covariance requires
n>=2 and s<n; n=1 point estimation can use s=1.

Use effective covariance availability
`se && jackknife != "none"`. If it is false, all covariance fields below
are NULL and metadata records `jackknife = "none"`. Do not treat
`se = TRUE, jackknife = "none"` as having standard errors.

Grid-mode evaluation is a later stage. A tuning object can supply the default
candidate now, but non-default candidates are not evaluated implicitly.
When grid mode is added, follow [design note 14](../design/14_integrated_tuning_grid_reuse.md)
and `candidate_fit()`; do not create a separate grid-estimator export.

## 3. Stable public result contract

For m reporting targets and n input observations, use the following core
payload. No dimensions disappear when m=1.

| Field | Type and meaning |
| --- | --- |
| `theta` | Numeric length m; NA for exact-zero or unrepresentable ratios |
| `se` | NULL without covariance; otherwise numeric length m |
| `components` | Numeric m by 2 matrix, columns `numerator`, `denominator` |
| `component_vcov` | NULL or list of m named 2 by 2 estimator-scale covariance matrices |
| `theta_vcov` | NULL or list of m named 1 by 1 estimator-scale covariance matrices |
| `theta_omega` | NULL or list of m 1 by 1 local-limit covariance estimates, N times theta_vcov |
| `x0` | Numeric m by k_x matrix, preserving target names with `.dnn_x0_with_names()` |
| `neff` | Numeric length m, 1/colSums(W^2), always returned |
| `status` | Data frame with one row per target and fields defined in section 4 |
| `row_ids` | Unique aligned input identities, length n |
| `meta` | n, m, reporting dimension p, k, s, method, distance, covariance mechanism, component roles, tie rule and tuning source |
| `.model` | Validated reporting design, two frozen signals, resolved controls and x_cols needed to evaluate new targets |

Use lists of matrices for covariance to match existing PLQR/GMM storage and
the shared covariance helper. This deliberately refines the scientific
notes' provisional array schema. Matrix names remain
`numerator`, `denominator`; CLATE metadata maps those roles to
`reduced_form`, `first_stage`. The ratio coefficient is named `theta`.

`neff` is always present because the feature explains both variance scales.
For compatibility with shared extras, a requested `return = "ess"` adds an
identical `ess` field. It never rescales covariance a second time.
Default W retention follows the existing wrappers; `return = character()`
omits optional W while preserving prediction state and all core fields.

Generic optional returns are `W,ranks,ord,dx,ess,signals,deletions,pseudo_values,tuning`.
CLATE additionally accepts `nuisance`. Unknown entries warn once and are
ignored through the shared policy. Signal matrices are n by 2; deletion and
pseudo-value arrays are n by 2 by m. Both signal columns remain in original
input-row order. A nuisance return uses the constructor's raw/used prediction
payload. Optional arrays can be large and must not be stored by default.

## 4. Per-target numerical statuses

Retain independent flags rather than encode all problems in one ordered
string:

| Field | Rule |
| --- | --- |
| `target_id` | 1..m, never derived from rounded target coordinates |
| `ratio_status` | `ok`, `zero_denominator`, `component_overflow`, `ratio_overflow` |
| `negative_denominator` | Finite denominator <0, independent of ratio_status |
| `near_zero_denominator` | abs(denominator) <= denominator_tol times max(abs(input denominator signal)); excludes exact zero |
| `covariance_status` | `not_computed`, `ok`, `zero_score_variance`, `numerical_failure` |
| `wald_available` | Finite ratio and finite strictly positive ratio variance |

Require denominator_tol to be a finite nonnegative scalar; zero disables only
the near-zero flag. The scale uses the supplied denominator signal and
preserves positive common scaling. An identically zero signal has an exact
zero denominator, not an undefined relative tolerance. A near-zero flag does
not suppress a representable ratio or change its standard error. It is a
numerical cancellation diagnostic, not a universal weak-IV threshold.

Components and their covariance remain available at an exactly zero
denominator. Set theta, se and the corresponding ratio covariance entries
to NA; never replace them with zero. A computed zero score variance gives
se=0 but no normal-reference interval. Overflow receives an explicit status
and NA derived output; invalid supplied non-finite arrays fail validation
globally. A numerical failure at one target must not erase other targets.

## 5. Classes and methods

The generic object has class `affine_ratio_result`. The IV wrapper has
classes `c("clate_result", "affine_ratio_result")`. It does not inherit
`cate_result`: CATE methods assume tau and a scalar signal.

Proposed S3 methods live in `R/user-dnn-affine-ratio-methods.R`:

- `coef.affine_ratio_result(object, ...)` returns a named length-m vector.
  Names are target_1 through target_m in original target order.
- `vcov.affine_ratio_result(object, target = NULL, ...)` follows existing
  PLQR/GMM: a single matrix for one target or a selected target, otherwise a
  list. These are marginal target variances, not cross-target covariance.
- `confint.affine_ratio_result(object, parm = NULL, level = .95, ...)`
  returns an m_selected by 2 matrix with conf.low/conf.high, retaining NA
  rows where wald_available is false. Missing covariance is an actionable
  error. Validate level and target indices as in existing CATE methods.
- `summary.affine_ratio_result(object, level = .95, ...)` returns a classed
  table containing both components, denominator SE, ratio, ratio SE, interval,
  N and statuses. Print never drops failed targets.
- `predict.affine_ratio_result(object, newdata = NULL, X0 = NULL, x0 = NULL,
  se.fit = FALSE, interval = c("none","confidence"), level = .95, ...)`
  relocalizes the stored two-signal state with the original order; it fits no
  nuisance functions and never evaluates new households' outcomes.

The summary table has classes
`c("summary.affine_ratio_result","data.frame")` and columns target_id,
numerator, denominator, numerator_se, denominator_se, theta, se,
conf.low, conf.high, neff, ratio_status, negative_denominator,
near_zero_denominator, covariance_status and wald_available. Without
covariance, summary still returns point estimates with NA uncertainty
columns; confint remains an error in that case.
`print.affine_ratio_result(x,...)` and
`print.summary.affine_ratio_result(x,...)` return x invisibly after
displaying component roles and status counts. Document and register both.

Prediction follows the existing CATE return grammar: numeric vector without
SE/interval, list(fit,se.fit) with SE, or fit/lwr/upr matrix with intervals.
Invalid targets retain NA entries. For full components/statuses at new targets,
make a documented fresh `affine_ratio_dnn()` call using stored input arrays;
do not hide metadata as incidental numeric-vector attributes.

CLATE-specific print and summary labels may delegate to the base methods with
component-role metadata. Neither print method claims causal identification.
`diagnostics.clate_result()` and the generic diagnostic method are specified
in Note 05. Plot methods and candidate views are separately gated extensions.

## 6. Contract review before coding

Check matrix and data-frame calls on the same aligned data; m=1, m>1 and
k_x=1; two targets with different numerical statuses; omitted optional
returns; prediction after saveRDS/readRDS; and all covariance scales. Update
these notes if a signature changes before implementation. A changed spelling
must propagate to tests, examples, help and the task sequence together.
