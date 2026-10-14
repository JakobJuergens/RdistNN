# RdistNN Design Notes

This directory contains the canonical ordered design-note corpus for the
`RdistNN` package.

The numbered files are the authoritative development policy set for this repo.
Older unnumbered notes have been retired so that humans and coding agents can
rely on a single explicit read order.

## Read order

1. `01_interface_standardization_decisions.md`
2. `02_code_design_principles.md`
3. `03_input_validation_contract.md`
4. `04_api_shapes_and_input_formats.md`
5. `05_error_messages_and_warnings.md`
6. `06_namespace_and_dependency_policy.md`
7. `07_lifecycle_and_deprecation_policy.md`
8. `08_comments_and_function_docs.md`
9. `09_testing_and_regression_policy.md`
10. `10_rcpp_and_compiled_code_policy.md`
11. `11_examples_vignettes_and_pkgdown_policy.md`
12. `12_release_and_ci_policy.md`
13. `13_existing_codebase_standardization_plan.md`
14. `14_integrated_tuning_grid_reuse.md`
15. `15_plqr_variance_estimation_rationale.md`
16. `16_affine_replication_contract.md`

## Structure

### Foundation
- `01_interface_standardization_decisions.md`
- `02_code_design_principles.md`
- `03_input_validation_contract.md`
- `04_api_shapes_and_input_formats.md`
- `05_error_messages_and_warnings.md`

### Package engineering
- `06_namespace_and_dependency_policy.md`
- `07_lifecycle_and_deprecation_policy.md`
- `08_comments_and_function_docs.md`
- `09_testing_and_regression_policy.md`
- `10_rcpp_and_compiled_code_policy.md`
- `11_examples_vignettes_and_pkgdown_policy.md`
- `12_release_and_ci_policy.md`

### Execution
- `13_existing_codebase_standardization_plan.md`
- `14_integrated_tuning_grid_reuse.md`
- `15_plqr_variance_estimation_rationale.md`
- `16_affine_replication_contract.md`

## Core package-wide decision

User-facing estimators may accept either:

- explicit matrices and vectors; or
- `data.frame`-style inputs, including tibbles.

After validation and interface resolution, all internal computation should proceed
on canonical numeric matrices and vectors.
