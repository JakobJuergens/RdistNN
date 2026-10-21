# Interface standardization decisions

This note owns package-wide interface decisions. See the
[shape contract](04_api_shapes_and_input_formats.md) for dimensions and input
routes, and [validation](03_input_validation_contract.md) for boundary checks.

## Public and internal boundaries

High-level estimators support documented matrix/vector and `data.frame`
interfaces. Tibbles use the ordinary data-frame route; there is no separate
tibble API. Document both routes for each high-level export.

Validate and resolve inputs once at the public boundary, then compute with
canonical numeric matrices/vectors. Lower-level primitives prefer explicit
numeric inputs and should not repeatedly branch on container type. Use shared
column extractors and target resolvers rather than wrapper-specific versions.
When `X0` and `x0` are supported, require exactly one and resolve it to a target
matrix. Preserve any explicitly documented default target behavior.

## Internal names

| Role | Names |
| --- | --- |
| Training and targets | `x_mat`, `x0_mat` |
| Weights and ranks | `W`, `R` or `rank_matrix` |
| Outcomes, treatments and signals | `y_vec`, `t_vec`, `psi`, `signal` |
| Nuisance vectors | `e`, `m`, `m0`, `m1` |
| Controls | `k`, `s`, `method`, `jackknife`, `distance` |

Follow these conventions in new wrappers/helpers rather than inventing synonyms.

## Results and optional returns

Return classed estimator objects with a stable main payload and documented
methods. Requested `return` fields add extras without changing the estimate.
Warn about unknown entries, list valid choices and ignore the unknown values
consistently across wrappers. Result classes and retained-state ownership are
specified in [the architecture](13_estimator_and_result_architecture.md).
