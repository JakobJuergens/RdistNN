# Meaningful Test-Coverage Audit Design

Date: 2026-07-22

## Objective

Improve confidence in `RdistNN` by testing important public and numerical
contracts that are currently untested or silently skipped. Coverage percentage
is a diagnostic for finding gaps, not the acceptance target.

The clean baseline from `covr::package_coverage()` is 86.38%. The audit also
ran the test suite against an installed copy of the package, because
source-tree testing exposes unexported symbols that users and `R CMD check` do
not see.

## Audit findings

### Installed-layout gaps

`test-user-cate-methods.R` is skipped in its entirety against an installed
package. Its guard checks whether `summary.cate_result`,
`confint.cate_result`, and `predict.cate_result` are visible as ordinary
functions. They are correctly registered S3 methods, but they are not exported
symbols. The public generic calls work; the guard prevents the tests from
proving that.

The custom-distance tests for CATE and PLQR also skip their main expectations
in an installed layout. Their public extension-point logic only needs exported
functions, but the test implementation also asks for the internal
`ord_to_ranks()` helper.

### Public scientific contracts with weak coverage

- `plqr_score_matrix()` has little direct testing beyond two theta-validation
  branches. Its defining observation-level orthogonal-score formula is a
  scientific contract.
- `np_demo_data()` and `np_demo_subsets()` have no dedicated test file. The
  public CATE and PLQR subset helpers are also uncovered.
- CATE and GMM result printing has weak or zero coverage. These methods are
  user-facing object contracts, including the distinction between available
  and unavailable standard errors and the display of overidentification
  diagnostics.

### Numerical-path gap

`score_plqr_component()` normally uses compiled score construction and subset
aggregation when compiled symbols are available. Existing tests establish the
formula in the active path, but do not compare the compiled path with the
pure-R fallback selected by `options(RdistNN.use_compiled = FALSE)`.

### Deliberate non-targets

The audit will not add tests solely to cover benchmark timing utilities,
package-startup hooks, optional third-party benchmark integrations, or every
defensive line. In particular, `.validate_same_n()` currently has no call
sites; directly testing an unused internal helper would raise coverage without
protecting active behavior.

## Proposed changes

### 1. Make installed-package tests honest

- Remove the symbol-visibility guard from the CATE method test file. The file
  will exercise `summary()`, `confint()`, and `predict()` through their public
  generics, so missing registration becomes a failure rather than a skip.
- Replace the custom-distance tests' dependency on `ord_to_ranks()` with a
  small test-local conversion from order vectors to rank vectors. Keep the
  contract expressed entirely through exported extension-point functions.

The before/after evidence is the installed-layout reporter: currently the CATE
method file produces zero expectations, and the CATE/PLQR custom tests report
"required APIs not available".

### 2. Add direct public-contract tests

For `plqr_score_matrix()`:

- compare every score column with the documented formula;
- check equivalent name-based and index-based column selection;
- check dimensions and stable numeric output;
- exercise meaningful schema, dimension, probability, and non-finite-input
  failures.

For the demo-data and subset helpers:

- check reproducibility under a fixed seed;
- check the documented data schema and target/truth dimensions;
- use zero noise to verify the oracle regression exactly for each supported
  NPR design;
- check NA/Inf injection and invalid input diagnostics;
- check that NPR, CATE, and PLQR subset matrices have valid indices,
  deterministic seeds, and the requested dimensions.

For result methods:

- verify that CATE and GMM print methods return their input invisibly and expose
  stable identifying lines;
- verify the GMM summary coefficient and overidentification schemas;
- verify output branches with and without standard errors and with and without
  finite overidentification statistics.

### 3. Add PLQR compiled/fallback parity

Construct one deterministic PLQR score and evaluate it with compiled dispatch
enabled and disabled. Compare:

- observation-level score contributions stored by `score_setup()`;
- one subset evaluation; and
- a batch of subset evaluations.

The test will skip only when the relevant compiled symbols are genuinely
unavailable. It will restore the user's option on exit.

## Test organization

- Repair the existing CATE method and CATE/PLQR custom-distance files rather
  than duplicating their scenarios.
- Add focused tests alongside the existing PLQR wrapper and compiled-parity
  suites.
- Add a dedicated NPR demo-data test file, using the existing IV demo-data
  tests as the structural pattern.
- Extend the existing GMM method tests with result-object assertions.

## Verification and acceptance

The changes are accepted when:

1. targeted tests first demonstrate the current skip or missing-contract gap
   and then pass after the test repair/addition;
2. the repaired files execute in an installed-package layout;
3. the full source-tree test suite passes;
4. the test dependency and lint checks pass;
5. a local installed-package check has no new errors or warnings; and
6. coverage is rerun and reported as an outcome, without treating any target
   percentage as a release gate.

Coverage that remains in optional benchmark integrations, inactive helpers, or
low-value formatting branches will be documented rather than filled with
artificial tests.
