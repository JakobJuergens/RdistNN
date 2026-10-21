# PLQR Variance-Estimation Rationale and Nuisance Notation

PLQR coefficient-level inference through `confint()`, `predict(se.fit = TRUE)`
and `summary()` requires a density-backed coefficient sandwich. This document
specifies its variance scales and supplied nuisance notation. For statistical
assumptions and derivations, see the
[conditional Z-estimation paper](https://jakobjuergens.com/research/#orthogonal-inference-conditional-z).

## 1. Moment and coefficient uncertainty

Scalar regression/CATE uncertainty and PLQR coefficient uncertainty have
different roles. CATE exposes exact/asymptotic single-scale DNN; two-scale
TDNN remains NP-regression functionality. A supported numerical calculation
must not be interpreted as a theorem for another score or weighting scheme.

For the PLQR quantile score, a componentwise jackknife SE for the localized
moment is not by itself a coefficient-level SE. Coefficient inference requires
a bread/meat/bread sandwich. The bread remains density-backed in all
coefficient-inference routes. The meat can be estimated either by the
analytic PLQR formula
`S = tau(1-tau) sum_i w_i tilde_D_i tilde_D_i'` or by a **full vector
jackknife covariance** of the localized PLQR moment (`moment_vcov`),
which is a valid candidate meat once sandwiched with the density-backed
bread. The package must — and does — distinguish moment diagnostics from
jackknife-meat coefficient inference.

For the analytic route, the local coupled-score sandwich uses the
following meat and density-backed Jacobian:

```
S(bx) = tau(1-tau) * sum_i w_i * tilde_D_i tilde_D_i'
J(bx) = sum_i w_i * q_i * tilde_D_i * D_i'        (q_i estimates the local conditional density)
```

`.plqr_direct_sandwich()` in `R/user-dnn-plqr-core.R` implements this
bread/meat/bread calculation. Preserve the distinction between its
limit-covariance scale and the estimator covariance after division by `neff`.

**Consequence for the package's public API:** `confint.plqr_result()` and
`predict.plqr_result(..., se.fit = TRUE)` only produce coefficient SEs/CIs
when `density`/`density_col` was supplied at fit time, because every valid
coefficient route needs the density-backed bread. Two meat estimators are
supported through `plqr_dnn(theta_meat = )`:

| Route | Bread | Meat | Coefficient inference? | Package fields |
|---|---|---|---|---|
| Analytic direct sandwich (`theta_meat = "analytic"`, default) | density-backed `jacobian` | analytic `score_covariance` | yes | `theta_vcov`, `theta_se`, `theta_omega` |
| Jackknife-meat sandwich (`theta_meat = "jackknife"`) | density-backed `jacobian` | full `moment_vcov` (delete-1 covariance of the localized vector moment) | yes | `moment_vcov`, `theta_vcov`, `theta_se`, `theta_omega` |
| Moment diagnostic only (no density input) | none | componentwise `moment_se` | **no** | `moment_se` |

The componentwise moment-level jackknife SE (`fit$moment_se` when
`density` was not supplied) remains available as a moment diagnostic, but
must **not** be treated as a valid basis for a coefficient-level CI on
`theta` by itself — the covariance structure across moment components is
lost, and no theorem covers that shortcut. The scaling conventions are
aligned so that under `theta_meat = "jackknife"`,
`moment_se = sqrt(diag(moment_vcov[[j]]))` and
`theta_vcov[[j]] = J^{-1} moment_vcov[[j]] J^{-T}` with
`theta_omega[[j]] = theta_vcov[[j]] * neff[[j]]`.

## 2. Nuisance notation mapping

The following table maps the current PLQR arguments to the coupled-score
notation and specifies the adjustment state for each supplied nuisance:

| RdistNN name         | Paper symbol                  | Meaning                                                                 | Domain (per paper)          |
|-----------------------|-------------------------------|--------------------------------------------------------------------------|------------------------------|
| `g_col` / `g`          | `g_0(U)`                       | Nonparametric partially-linear shape nuisance; enters the threshold as `D'theta_eval + g_0(U)`. It is **not itself a quantile function** — the sum `D'theta_0 + g_0` is the conditional quantile. | `U = (X, Z)`, the full adjustment state — not `X` alone. |
| `mD_cols` / `mD`       | `m_0(U) = E[D \| U]`            | Conditional mean used to residualize the design: `tilde_D = D - m_0(U)`. | `U = (X, Z)` |
| `density` / `density_col` | `q_{0,bx}(D, Z)`            | Conditional density of `Y` at the fitted quantile value; enters **only** the Jacobian `J(bx)`, never the moment or `S(bx)`. | function of `(D, Z)` given `X = bx` |

Both `g_col`/`g` and `mD_cols`/`mD` are supplied externally by the user
(no nuisance is estimated internally by `plqr_dnn()`). They can depend on
the full adjustment state, rather than on reporting coordinates `X` alone.
Keep help, vignettes and implementations aligned with this table.

## 3. Scope note

This numerical contract does not prove Hájek dominance or jackknife
consistency for a PLQR score kernel. Statistical inference requires the
corresponding design, learner, density, bias and covariance assumptions.
Preserve the `confint`/`predict` gate on `density`/`density_col` and avoid
**componentwise moment-SE shortcuts** to coefficient CIs. The full
jackknife-meat route requires the density-backed bread, its appropriate
statistical justification, independent tests and consistent scale documentation.
