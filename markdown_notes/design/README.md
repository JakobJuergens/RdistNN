# Maintainer design reference

These documents define the maintained interface, numerical and engineering
contracts for RdistNN. Use [AGENTS.md](../../AGENTS.md) for editing rules,
[CONTRIBUTING.md](../../CONTRIBUTING.md) for commands and checks, and
[the package scope](../rdistnn_roadmap.md) for current capabilities and limits.

The [paper](https://jakobjuergens.com/research/#orthogonal-inference-conditional-z)
owns statistical derivations. These notes retain formulas, numerical
choices and scope limits needed to implement the package correctly. Read the
contract for the surface being changed; the numbered files are references,
not a release task list.

## Interface and validation

1. [Interface decisions](01_interface_standardization_decisions.md)
2. [Code design](02_code_design_principles.md)
3. [Input validation](03_input_validation_contract.md)
4. [Shapes and input formats](04_api_shapes_and_input_formats.md)
5. [Errors and warnings](05_error_messages_and_warnings.md)

## Engineering policy

6. [Namespace and dependencies](06_namespace_and_dependency_policy.md)
7. [Lifecycle and deprecation](07_lifecycle_and_deprecation_policy.md)
8. [Comments and function documentation](08_comments_and_function_docs.md)
9. [Testing and regression](09_testing_and_regression_policy.md)
10. [Rcpp and compiled code](10_rcpp_and_compiled_code_policy.md)
11. [Examples, vignettes and pkgdown](11_examples_vignettes_and_pkgdown_policy.md)
12. [Release and CI](12_release_and_ci_policy.md)

For reporting changes, read notes 05, 09 and 10 together: they record the
intended message identities/severity, selective recovery, grouping reuse,
native boundaries and behavioral/performance tests. These requirements guide
implementation; they do not claim that every existing call site has already
been migrated. Source, help and the affected regression checks establish
implemented behavior.

## Architecture and numerical contracts

13. [Estimator and result architecture](13_estimator_and_result_architecture.md)
14. [Integrated tuning-grid reuse](14_integrated_tuning_grid_reuse.md)
15. [PLQR variance and nuisance notation](15_plqr_variance_estimation_rationale.md)
16. [Affine replication integrity](16_affine_replication_contract.md)
17. [CLATE and paired-ratio inference](17_clate_ratio_inference.md)
18. [Geometry and covariance arithmetic](18_geometry_and_covariance_arithmetic.md)
19. [Weighted accumulation](19_weighted_accumulation.md)

User-facing estimators accept the documented matrix/vector or data-frame
interface and then compute with canonical numeric inputs. Shared generics and
builders preserve distinct estimands and uncertainty scales. New contributions
must maintain these contracts and their independent numerical tests.

Keep this reference focused on current decisions and rationale. Session
instructions, task checklists, local paths and execution receipts belong in
private working records rather than public maintainer documentation.
