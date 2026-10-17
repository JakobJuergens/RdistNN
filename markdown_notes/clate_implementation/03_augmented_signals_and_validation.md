# 03 — Augmented signals, clipping and input validation

[Index](README.md) · [Previous](02_public_interfaces_and_results.md) ·
[Next: computation](04_localization_and_deletion_engine.md)

**Status:** proposed constructor and validation contract. The implementation
owner is `R/score-clate.R`; wrapper interface resolution belongs in
`R/user-dnn-clate.R`. No learner is fitted in either file.

## 1. Constructor output and naming

`iv_augmented_signals()` takes the signature in Note 02. Return a classed
`iv_signals` list with length-n phi_y and phi_d vectors, row_ids,
propensity = list(raw,used,clipped), nuisance = list(m_y0,m_y1,m_d0,m_d1),
and meta = list(n,clip,one_sided,training). The clipped field is a length-n
logical vector. Retain raw and used propensity even when no value changes.
The constructor must not keep a data frame whose columns can be mutated
independently of the canonical arrays.

The a suffix indexes instrument assignment. For example m_d0 predicts D
given A=0,V, not the outcome among untreated participants. A is the
instrument, D the treatment, e the instrument propensity. This mapping must
appear in argument help and examples.

Training is a descriptive scalar from
`supplied,oracle,full_sample,oob,cross_fit`. It records the caller's
declaration, not a fact inferred from fitted values. An oracle label does not
make estimated predictions oracle; a cross_fit label does not prove held-out
construction. Future adapters must supply their own reproducible training
records and row/fold identities.

## 2. Validation sequence

Follow the package's numbered input-validation policy, in this order:

1. Resolve mutually exclusive sample interfaces and required arguments.
2. Validate containers and column selectors through the shared extractors.
3. Determine n and require nonempty reporting/sample inputs.
4. Check numeric scalar-vector shape and equal lengths. A matrix with one
   column is not silently accepted as an arbitrary sample vector.
5. Distinguish NA from NaN/Inf using the existing helpers.
6. Check binary A,D; propensity range; prediction ranges; one-sided coding;
   clip and other scalar controls.
7. Canonicalize to numeric matrices/vectors, preserving row order and IDs.
8. Construct signals and verify finite results before localization.

Use `.validate_binary01()`, `.validate_vector_nrow()`,
`.abort_if_na_or_nonfinite()`, `.dnn_col_vec()` and
`.dnn_col_mat()` where their exact signatures fit. New shared validation
is limited to a small IV-specific resolver and row-identity helper.
Do not add stopifnot() or duplicate the entire existing extractor family.

For lower-level constructor vectors, errors name Y,D,A,e or the exact
prediction argument. Data-frame wrapper errors name the corresponding
*_col selector and missing column. Require numeric/integer 0/1 indicators;
reject factors, character labels, logical arrays and -1/1 coding. Apply
numeric/type validation before the existing binary helper, which itself
coerces with as.numeric(). Otherwise a factor or logical vector could
acquire unintended semantics before the binary check.

The current vector-length helper also coerces numeric matrices to vectors.
Add an explicit is.null(dim(x)) guard for constructor/sample-vector inputs
before using it. This enforces the proposed vector contract without
changing that shared helper's behavior for existing estimators.

## 3. Propensity and response-prediction conventions

Require finite e in [0,1] before clipping. Endpoints are admissible raw
predictions and are mapped to [clip,1-clip]. Values outside [0,1] are errors.
Require a finite scalar 0<clip<1/2, initially 1e-8 to match existing wrappers.
Reusing `.cate_validate_and_clip_propensity()` is acceptable through a
small internal adapter; a name-neutral refactor is optional and must not
change CATE behavior or errors incidentally.

Both phi_y and phi_d use exactly the same e_used array. Do not clip a
propensity again under a different default inside one branch. Test
idempotence with the same clip and exact recording of endpoint changes.

Treatment-regression predictions m_d0,m_d1 must be in [0,1] in the initial
constructor contract. Caller-side projection is explicit. A later adapter
may project raw regression outputs and record that transformation; do not
silently alter only one supplied nuisance pair in the core constructor.
Outcome predictions and Y must be finite, with no default winsorization.
This permits computing practical financial-outcome signals, but does not
certify the bounded-learner theorem. Prediction clipping cannot make an
unbounded Y bounded.

If a bounded learner is used, its adapter owns projection of outcome
predictions and the model's bounded response. Transforming Y, capping Y or
changing its population changes the estimand and belongs in the caller's
analysis specification. No rows are dropped by the constructor.

## 4. Formula and independent audit expression

Compute for each response R:

$$
\phi_R=m_{R1}-m_{R0}
     +A(R-m_{R1})/e_{\mathrm{used}}
     -(1-A)(R-m_{R0})/(1-e_{\mathrm{used}}).
$$

The code may share an internal augmentation helper between the two
responses. Independent test expectations must expand the expression
directly, rather than call that same helper twice and compare it with itself.
Reject non-finite resulting signals with an argument-specific explanation
that finite extreme inputs caused arithmetic overflow. Do not silently
replace them by capped signals.

The exact fresh-mean error for fixed candidate functions is

$$
E[\phi_R(\eta)-\phi_R(\eta_0)\mid V]
=(e-e_0)\{(m_{R1}-m_{R1,0})/e+
          (m_{R0}-m_{R0,0})/(1-e)\}.
$$

Both regression-error terms inside the braces have positive signs.
The complete score has remainder theta R_D-R_Y. Tests enumerate a small
finite response/assignment distribution and evaluate both sides
independently. Correct e or both correct response regressions sets this mean
remainder to zero. These tests establish an identity; they do not prove
inference with arbitrary inconsistent nuisance learners.

## 5. Explicit one-sided compliance

If one_sided=TRUE, require D<=A at every row. Omitted m_d0 is replaced by
an explicit zero vector. If m_d0 is supplied, require it to be identically
zero; do not ignore conflicting predictions with a warning.
The simplified signal is

$$
\phi_D=m_{D1}+A(D-m_{D1})/e_{\mathrm{used}}.
$$

The outcome still needs both instrument-arm regressions. Outcome Y is not
zeroed among ineligible observations. The simplified and general formulas
must agree when D=0 for A=0 and m_d0=0.

Do not infer the structural restriction from a sample with no observed
violations. One-sided compliance must be selected by the user. With
A=1,D=0,e=.5,m_d1=.6, the simplified first-stage signal is -.6.
Retain it. If one_sided=FALSE, a non-NULL length-n m_d0 is required.

## 6. Alignment and reproducibility

Default row_ids are seq_len(n). Accept unique nonmissing character or
integer/numeric identifiers with finite numeric values; reject duplicate,
missing, empty-character or malformed IDs. IDs are identities, not distance
tie-breakers. The current rank rule uses input-row position.

One row_ids vector proves the identity of the resulting aligned bundle; it
cannot prove that every independently supplied vector was previously
ordered correctly. Documentation must say this explicitly. A later
prediction-object adapter can require each input's IDs and reject mismatches
before construction. Do not sort independent arrays by their own ID vectors
without an explicit shared join contract.

Column extraction preserves data-frame order, and no missing-value handling
deletes a subset of rows. If a user filters a sample, all observed arrays,
predictions and reporting design must be filtered together before the call.
Every optional output uses the same order and IDs.

## 7. Concrete validation and invariance tests

| Input or transformation | Required behavior |
| --- | --- |
| Missing A or one outcome-arm mean | Error naming the missing argument before weights are computed |
| e=c(0,1), clip=.1 | e_used=c(.1,.9), both clipped flags TRUE |
| e outside [0,1] | Error before clipping |
| m_d1 outside [0,1] | Error; an adapter may project it explicitly |
| one_sided=TRUE with D>A | Error retaining the original row identity in the diagnostic |
| Missing m_d0 in unrestricted mode | Error |
| Full instrument relabeling | Swap all arm means, A to 1-A and e to 1-e; both signals negate |
| Perfect compliance D=A, m_d0=0,m_d1=1 | phi_d is identically 1 for every admissible e |
| Add c to Y,m_y0,m_y1 | phi_y unchanged |
| Multiply Y,m_y0,m_y1 by b | phi_y multiplied by b; include b<0 and b=0 |
| Eligible nonparticipant in the example above | Negative phi_d retained exactly |
| Independently permuted predictions without IDs | Not detectable from lengths alone; explain caller responsibility |

Validation errors are global misuse errors. Statistical/numerical statuses
from a valid local ratio belong to Note 05. Keeping those two failure classes
separate prevents accidental exclusion of unfavorable targets or simulations.
