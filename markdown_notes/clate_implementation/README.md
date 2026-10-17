# Conditional LATE implementation contracts and acceptance

Prepared October 17, 2026; implementation authorized subsequently in the same
session. **Status: implemented; local acceptance complete.** The nine design
notes below specify the augmented conditional local average treatment effect
(CLATE) feature. Tasks 1–8 are implemented, with 5,395 passing assertions,
built installed-package guides and a checked source archive. Independent
review, hosted platform acceptance and publication remain the final gates. No CLATE
production campaign or SIPP participation analysis is included. This navigation
index is not a tenth substantive design note.

The implementation computes two augmented signals, localizes both
with the same exact distributional nearest-neighbor weights, solves their
ratio equation, and estimates uncertainty from paired frozen-score deletions.
The plan includes a complete user workflow with externally supplied nuisance
predictions. Learner-training adapters are a subsequent, separately testable
stage. The bounded reference adapter is implemented internally; practical
automatic fitting remains future work.

## Read and execute in this order

| Note | Responsibility and concrete deliverable |
| --- | --- |
| [01 — Estimator, scope and source map](01_estimator_scope_and_source_map.md) | Scientific specification, estimand, reuse boundaries, current implementation owners |
| [02 — Public interfaces and result contracts](02_public_interfaces_and_results.md) | Implemented signatures, dual interfaces, shapes, metadata, S3 methods and statuses |
| [03 — Augmented signal construction](03_augmented_signals_and_validation.md) | Five prediction arrays, clipping, one-sided compliance, alignment and validation |
| [04 — Exact localization and deletion engine](04_localization_and_deletion_engine.md) | Shared-core adapters, rank conventions, paired covariance, complexity and numerical stability |
| [05 — Inference and diagnostics](05_inference_confidence_sets_and_diagnostics.md) | Covariance scales, Wald reporting, denominator behavior, optional confidence-set inversion |
| [06 — Learners, model and study program](06_learners_bounded_model_and_studies.md) | Explicit bounded benchmark, ridge verification contract, future adapters, simulations and SIPP separation |
| [07 — Unit and integration tests](07_unit_integration_and_regression_tests.md) | Named suites, independent oracles, exact fixture, failures, fallback parity and acceptance rules |
| [08 — Help, vignettes and site](08_documentation_examples_and_vignettes.md) | Roxygen topics, executable examples, two articles, pkgdown, accessibility and documentation checks |
| [09 — Staged implementation and release](09_implementation_sequence_and_release.md) | File ownership, dependency order, test-first tasks, commands and completion gates |

Read 01–03 before settling an API. Notes 04–05 define the numeric contract;
07 is its acceptance specification. Notes 06 and 08 distinguish a usable
package feature from a theorem claim or empirical finding. Note 09 is the
implementation handoff and acceptance checklist. Its original notes-only
authorization boundary was superseded by the explicit implementation request;
original prospective language below remains design provenance, not current
export or test status.

## Authority and inspected snapshots

The inspected RdistNN baseline is
`906789382e813c7a87c33be906d44755d2474436`, version `0.1.0.9000`, on
`master`. The scientific packet is in the sibling Cond_Z_Estimation
repository, at `5242cc1`: `markdown_notes/affine_scores/CLATE_Affine_DNN_Research_Note.md`
and the six files indexed by `markdown_notes/affine_scores/clate_exploration/README.md`.
Its mathematical TeX is unchanged from `1501105` through that inspected head.
The [affine setup](https://github.com/JakobJuergens/Cond_Z_Estimation/blob/1501105577ac48d9052f6b3193c6a148160c7fae/Paper/Parts/affine_scores/02_Setup.tex),
[conditions](https://github.com/JakobJuergens/Cond_Z_Estimation/blob/1501105577ac48d9052f6b3193c6a148160c7fae/Paper/Parts/affine_scores/04_Conditions.tex)
and [theorem](https://github.com/JakobJuergens/Cond_Z_Estimation/blob/1501105577ac48d9052f6b3193c6a148160c7fae/Paper/Parts/affine_scores/05_Results.tex)
are the mathematical dependencies. Their availability may require repository
access. Essential formulas and contracts are reproduced here, so implementation
does not require the sibling checkout or a network connection.

[AGENTS](../../AGENTS.md) and the numbered [design corpus](../design/01_interface_standardization_decisions.md)
govern engineering. This collection owns the feature-specific contract;
the implementation adds the new exports without superseding that policy. The
[live roadmap](../rdistnn_roadmap.md) owns open-work routing.

## Fixed decisions and later extensions

The implemented feature has a scalar target, binary instrument and
treatment, externally supplied predictions, exact without-replacement DNN,
`k = 1L`, one common order, paired frozen-score covariance and pointwise
reporting. It supports both matrix and data-frame inputs. A ratio with a finite
negative denominator remains a computed ratio; exact-zero and numerical
failures remain visible with their components.

Practical learners, adaptive tuning, weak-instrument-uniform guarantees,
multiple instruments, continuous treatment, simultaneous bands and generalized
affine systems require additional contracts. They are not silently included.
Optional confidence-set inversion has its own stage and geometry tests.

Existing CATE/regression onboarding, experimental PLQR/IV-GMM behavior and the
frozen `inst/paper-replication/affine/` bundle remain separate owners. The
CLATE feature has an experimental lifecycle designation.

## Requirements and evidence

| Requirement | Authoritative implementation and acceptance evidence |
| --- | --- |
| Augmented signals, alignment and one-sided checks | `R/score-clate.R`, `test-clate-signals.R`; exact product remainders, clipping and all public vector checks |
| Common exact weights and paired covariance | `R/user-dnn-affine-ratio.R`; subset enumeration, four-row variance 34, ESS 18/7, limit covariance 612/7, explicit-loop/R/compiled parity |
| Public supplied-prediction wrapper | `R/user-dnn-clate.R`; matrix/data-frame/tibble, composition, instrument relabeling, perfect compliance and signed/zero denominator suites |
| Methods, retained prediction and diagnostics | `R/user-dnn-affine-ratio-methods.R`, `R/user-dnn-clate-diagnostics.R`; target shapes, serialization, unavailable uncertainty and independent concentration references |
| Bounded benchmark and designated ridge objective | `R/demo-data-clate.R`, `R/sim-clate-learners.R`; latent integration, core constants, heterogeneous target, dense normal equations and all one-/two-deleted refits |
| Set-valued moment inversion | `R/clate-confidence-sets.R`; every quadratic geometry, PSD/singular inputs, membership, rescaling, extreme confidence levels and subnormal covariance |
| Executable documentation | Both installed-package articles, generated help/README and full pkgdown build pass; 131 public topics indexed; 574 local links across 13 changed pages checked; browser-reviewed tables/formulae/set output |
| Existing workflows and frozen bundle | 882 test blocks, 5,395 passing assertions, no failures/errors/warnings, one opt-in visual skip; source archive has no errors/warnings and one 11.4 MB size NOTE; alignment, lint and pinned bundle verification pass |
| Supported platforms and publication | Hosted R-CMD-check matrix and direct `master` push pending final review |

Completed source slices: `a6a2504` (signals), `a84097b` (paired engine/wrapper),
`01c2756` (methods/diagnostics), `c95f4a9` (bounded generator/ridge), `8a6a139`
(confidence sets and boundary regressions), `53da721` (zero-score, namespace
and expanded edge-grid regressions). Focused test commands and the
remaining acceptance checks are in Note 09. These computational results do
not establish a practical learner's asymptotic rates or new empirical evidence.
