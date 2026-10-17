# Conditional LATE implementation planning

Prepared October 17, 2026. **Status: planning only.** This collection contains
nine substantive notes for a future RdistNN implementation of the augmented
conditional local average treatment effect (CLATE) estimator. No proposed
export, learner, test suite, vignette or numerical campaign is implemented by
these notes. The navigation index is not a tenth substantive note.

The intended implementation computes two augmented signals, localizes both
with the same exact distributional nearest-neighbor weights, solves their
ratio equation, and estimates uncertainty from paired frozen-score deletions.
The plan includes a complete user workflow with externally supplied nuisance
predictions. Learner-training adapters are a subsequent, separately testable
stage; the plan does not stop at an undocumented array helper.

## Read and execute in this order

| Note | Responsibility and concrete deliverable |
| --- | --- |
| [01 — Estimator, scope and source map](01_estimator_scope_and_source_map.md) | Scientific specification, estimand, reuse boundaries, current implementation owners |
| [02 — Public interfaces and result contracts](02_public_interfaces_and_results.md) | Proposed signatures, dual interfaces, shapes, metadata, S3 methods and statuses |
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
implementation handoff and does not authorize starting implementation.

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
govern engineering. This collection proposes a feature-specific contract;
it does not change existing exports or supersede that policy. The
[live roadmap](../rdistnn_roadmap.md) owns open-work routing.

## Fixed decisions and later extensions

The proposed first complete feature has a scalar target, binary instrument and
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
proposed CLATE feature begins with an experimental lifecycle designation.
The present task publishes these planning documents only.
