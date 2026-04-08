# Interface Standardization Decisions

This note records the package-wide interface decisions that should be treated as
the default going forward.

## 1. Official input policy

### Public user-facing estimators
User-facing estimators may accept either:

- explicit matrices and vectors; or
- `data.frame`-style inputs, including tibbles.

This applies to functions such as:

- `np_dnn()`
- `cate_dnn()`
- `plqr_dnn()`
- future high-level estimation wrappers

### Internal computation
After validation and interface resolution, all internal computation should
proceed on canonical numeric matrices and vectors.

This is the main standardization decision for the code base.

## 2. Tibbles

Tibbles are allowed only through the `data.frame` pathway.
There is no separate tibble-specific API.

Rule:
if an object satisfies the `data.frame` interface and column extraction works
correctly, it is admissible as a user-facing tabular input.

Internally, tabular inputs should be converted to numeric matrices and vectors.

## 3. Public boundary vs internal boundary

### Public boundary
Flexible and user-friendly.

### Internal boundary
Strict and canonical.

This means:

- the outer API may accept multiple formats;
- deeper helpers should not keep branching on container type.

Once inputs are validated, convert them to canonical matrix/vector forms and
stay there.

## 4. Canonical internal naming

Standardize the following names across wrappers and helpers.

### Matrices
- `x_mat` for the training/design matrix
- `x0_mat` for resolved target matrix
- `W` for weight matrix when returned or constructed explicitly
- `R` or `rank_matrix` for batch rank matrices when relevant

### Vectors
- `y_vec`
- `t_vec`
- `psi`
- `e`
- `m`
- `m0`
- `m1`
- `signal`

### Scalars / control arguments
- `k`
- `s`
- `method`
- `jackknife`
- `distance`

Do not invent one-off naming conventions inside new wrappers.

## 5. Target resolution

Whenever a function supports both `X0` and `x0`:

- exactly one should be allowed;
- target resolution should happen through a shared helper;
- the result should be a canonical target matrix form.

Do not reimplement target resolution ad hoc in multiple wrappers.

## 6. Data-frame column extraction

All column-name-based extraction should pass through shared internal helpers.

Recommended helper family:

- single-column numeric vector extractor and validator
- multi-column matrix extractor and validator

This avoids divergent missing-column behavior across estimators.

## 7. Lower-level function preference

Lower-level functions should prefer matrix/vector inputs, not `data`.
That keeps the lower-level API explicit and reduces branching.

As a default rule:

- high-level estimators may accept `data=...`;
- lower-level computation helpers should not.

## 8. Return-value standardization

High-level wrappers should follow the same `return=` policy.

Recommended default while the package remains experimental:

- unknown `return` values produce a warning and are ignored;
- valid `return` values add optional fields without changing the main result.

This behavior should be consistent across all wrappers that expose `return=`.

## 9. Result object policy

User-facing estimators should return classed objects with:

- a stable main estimate payload;
- optional extra fields added predictably;
- stable method support where appropriate.

Avoid returning raw unclassed lists from exported estimation workflows.

## 10. Documentation rule

Every exported high-level estimator should document both:

- the matrix/vector interface; and
- the `data.frame` interface.

The docs should explicitly state that tibbles are accepted via the
`data.frame` path if they behave like ordinary data frames for the required
columns.
