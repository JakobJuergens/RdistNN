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
- numerical intent up to reasonable tolerance.

Compiled code should not quietly drift from the reference behavior.

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

At minimum, test:
- one ordinary valid case;
- one edge case;
- shape agreement;
- numerical agreement within tolerance;
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
