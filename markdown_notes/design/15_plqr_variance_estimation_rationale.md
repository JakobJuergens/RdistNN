# PLQR Variance-Estimation Rationale and Nuisance Notation

This note records why `plqr_dnn()`'s coefficient-level inference (`confint()`,
`predict(..., se.fit = TRUE)`, `summary()`) is gated on the direct coefficient
sandwich (`density`/`density_col`) rather than accepting the ordinary
delete-1 jackknife SE the way `np_dnn()` and `cate_dnn()` do, and maps
RdistNN's PLQR argument names onto the theory papers' notation. This is a
permanent, theory-driven design decision, not a temporary limitation.

## 1. Why PLQR requires the direct sandwich

`np_dnn()` and `cate_dnn()`'s delete-1 jackknife standard errors are backed
by `JK_GenUStat`'s jackknife ratio-consistency theorems
(`thm:JK_Consistency_compl`, `thm:JK_Consistency_incompl`,
`Main_GenUStat.tex`), specialized to the DNN/TDNN kernel in
`thm:JK_Consistency_TDNN` (`Main_TDNN.tex`). That specialization works
because the underlying score (`Y` for `np_dnn()`, the AIPW/Robinson
pseudo-outcome for `cate_dnn()`) is a function of a single observation —
the same "selector times a function of a single observation" shape the
DNN/TDNN Hájek-dominance proof verifies
(`Supp_DNN_HajekDominance.tex:104-108`).

`JK_GenUStat`'s own conclusion is explicit that this does **not** extend to
conditional Z-estimation for free:

> "A related direction is conditional Z-estimation, where DNN localizer
> weights enter a score-based rather than regression-based kernel; the
> Hájek-dominance machinery applies once the score-kernel first-projection
> variance is controlled ... The present paper does not resolve those
> broader questions." (`Main_End.tex:9-13`)

The PLQR moment (a quantile check-function score) is exactly this
unresolved case. There is no proof anywhere in this repository's source
papers that a delete-1 jackknife SE on the PLQR *moment* has a valid
translation into a coefficient-level (`theta`) confidence interval.

`Cond_Z_Estimation`'s active "coupled-score" theory
(`Paper/Parts/coupled_score/05_Applications/05b_Coupled_Nonsmooth_Applications.tex`)
instead derives a complete, self-contained asymptotic theory for PLQR that
does *not* route through jackknife at all: Corollary
`cor:coupled_plqr_worked_point_inference` gives
`sqrt(neff(x)) (hat theta(bx) - theta_0(bx)) -> N(0, Omega_PLQR(bx))`,
`Omega_PLQR(bx) = J(bx)^{-1} S(bx) J(bx)^{-T}`, with condition **PLQR-E**
specifying the direct plug-in estimators

```
S(bx) = tau(1-tau) * sum_i w_i * tilde_D_i tilde_D_i'
J(bx) = sum_i w_i * q_i * tilde_D_i * D_i'        (q_i estimates the local conditional density)
```

`.plqr_direct_sandwich()` (`R/user-dnn-plqr.R:56-114`) implements exactly
this: `S` at lines 71-72, `J` at line 74, `Omega = J^{-1} S J^{-T}` at lines
90-91, `vcov = Omega / neff` at line 92. This was verified line-by-line
against the paper's formulas during the PLQR-maturity planning pass
(2026-07-08) and matches exactly.

**Consequence for the package's public API:** `confint.plqr_result()` and
`predict.plqr_result(..., se.fit = TRUE)` only produce coefficient SEs/CIs
when `density`/`density_col` was supplied at fit time (the direct-sandwich
route, reported in `theta_se`/`theta_vcov`), because that is the only route
this repository's theory has proven valid. The moment-level jackknife SE
(`fit$moment_se` when `density` was not supplied) remains available as a
moment diagnostic, but must **not** be treated as a valid basis for a
coefficient-level CI on `theta` — no theorem in this repository's paper
corpus covers that translation.

## 2. Nuisance notation mapping

RdistNN's PLQR argument names were aligned with the paper's current
"coupled-score" notation in the 2026-07-08 pre-public cleanup
(`markdown_notes/rdistnn_public_release_cleanup_note.md`). This table is the
authoritative mapping for anyone connecting code to theory:

| RdistNN name         | Paper symbol                  | Meaning                                                                 | Domain (per paper)          |
|-----------------------|-------------------------------|--------------------------------------------------------------------------|------------------------------|
| `g_col` / `g`          | `g_0(U)`                       | Nonparametric partially-linear shape nuisance; enters the threshold as `D'theta_eval + g_0(U)`. It is **not itself a quantile function** — the sum `D'theta_0 + g_0` is the conditional quantile. (Before the 2026-07-08 pre-public cleanup this input was misleadingly named `q_col`/`q`.) | `U = (X, Z)`, the full adjustment state — not `X` alone. |
| `mD_cols` / `mD`       | `m_0(U) = E[D \| U]`            | Conditional mean used to residualize the design: `tilde_D = D - m_0(U)`. | `U = (X, Z)` |
| `density` / `density_col` | `q_{0,bx}(D, Z)`            | Conditional density of `Y` at the fitted quantile value; enters **only** the Jacobian `J(bx)`, never the moment or `S(bx)`. | function of `(D, Z)` given `X = bx` |

Both `g_col`/`g` and `mD_cols`/`mD` are supplied externally by the user
(no nuisance is estimated internally by `plqr_dnn()`), so nothing in the
implementation actually restricts them to depend on `X` alone — the
restriction was only in the documentation's phrasing. The roxygen and
vignette wording were corrected to match this table in the same pass that
introduced this note; this note is the source of truth going forward.

## 3. Scope note

This note does not attempt to prove Hájek dominance or jackknife
consistency for the PLQR score kernel — per `Main_End.tex`, that remains
open in the source papers. It only records why the package's existing
`confint`/`predict` gating on `density`/`density_col` is the theoretically
correct behavior, not a placeholder limitation, so a future contributor
does not "fix" the asymmetry with `cate_dnn()` by adding an unproven
jackknife-based coefficient CI.
