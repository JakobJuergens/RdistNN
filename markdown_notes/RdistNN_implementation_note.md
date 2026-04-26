# RdistNN package audit note: application completeness, remaining gaps, and implementation priorities

*Created: 2026-04-09. Updated: 2026-04-26 after the diagnostics layer landed and the manuscript-facing simulation orchestration moved into `Cond_Z_Estimation/Codes/Simulations/`.*

## Purpose

This note summarizes the current state of the `RdistNN` package relative to the functionality needed for the paper applications, and identifies the remaining gaps in an order that is useful for implementation planning with codex.

The goal here is not to restate everything the package already does in abstract terms. The goal is to answer the practical question:

> Does the package now basically support the estimator workflows actually used in the paper applications, and if not, what still needs to be built next?

---

## Bottom line

The package is now **functionally complete for the current baseline simulation workflows used by the paper**, but it is **not yet fully complete as a polished end-to-end paper-facing reporting layer**, and it is **not yet a complete package implementation of the full general conditional Z-estimation framework**.

That distinction matters.

### For the applications currently shown in the paper

The answer is basically:

- **yes on the estimator side**,
- **mostly yes on the diagnostics/reporting side**, because `diagnostics()` is now a live public API for `np_dnn()`, `cate_dnn()`, and `plqr_dnn()`,
- and **partially no on the higher-level simulation-reporting side**, because there is still no package-level tuning-sensitivity helper that reruns fits across the candidate `\kappa` grid and summarizes estimate / SE / CI / diagnostics in one object.

### For the full general framework

The answer is still:

- **not yet**.

The package now has substantial modular infrastructure and several strong user-facing wrappers, but it does **not yet expose a single polished general conditional Z-estimation front end** that takes a score / moment description, localizes, solves, and reports inference in a unified way.

---

## Key evidence from the recent repo state

A recent major commit added a unified kernel-order tuning API and integrated it into the main DNN wrappers and simulation workflow. In particular, the commit titled:

- `Add kernel-order tuning API and registry`
- commit SHA: `7890822365548546f212c281751e1ba98999980f`

adds:

- a tuning registry and tuning object API,
- `tuning_dnn_paper()` as the built-in paper-default undersmoothing rule,
- `tuning_eval()` plus grid support,
- integration of `tuning` into `np_dnn()`, `cate_dnn()`, `plqr_dnn()`, and `sim_cate_run()`,
- attachment of tuning metadata to returned fit objects,
- tests and vignettes for this interface.

This is a real completeness improvement, not just a cosmetic cleanup.

The current package surface also includes:

- `np_dnn()` for localized nonparametric regression,
- `cate_dnn()` for localized CATE estimation,
- `plqr_dnn()` for localized PLQR orthogonal moments,
- `diagnostics()` plus package vignettes for diagnostics/reporting,
- a score API (`new_score()`, `as_score()`, `score_eval()`, `score_diff_path()`, etc.),
- simulation and replication code for the CATE paper workflow.

---

# 1. What is already in place

## 1.1 Main paper-facing estimator wrappers

The package now has credible user-facing wrappers for the main application-style workflows currently visible in the paper and vignettes:

- `np_dnn()` for localized conditional mean estimation,
- `cate_dnn()` for localized treatment-effect estimation,
- `plqr_dnn()` for localized orthogonal PLQR moments.

These are no longer just low-level building blocks. They provide a reasonably coherent user-facing layer.

### What these wrappers already support

Across the current wrapper layer, the package supports:

- matrix and data-frame interfaces,
- exact and asymptotic DNN localization,
- delete-1 jackknife standard errors,
- integration of the new tuning API,
- optional return of weights and related intermediate objects,
- optional return of effective sample size and neighborhood-level information.

That means the package already covers most of the estimator mechanics that are needed for the paper’s existing applications.

---

## 1.2 Tuning is now much more complete

The new tuning layer is a meaningful step forward.

Before this addition, the package was closer to “the estimator works if you manually supply `s` and know what you are doing.”

Now it has:

- a formal tuning object abstraction,
- a built-in paper-default rule,
- support for candidate grids,
- registry-based reuse,
- wrapper integration,
- tuning metadata attached to outputs.

This improves both:

- user experience, and
- reproducibility of the paper-facing workflows.

It also gives a clean place to build diagnostics and sensitivity summaries next.

---

## 1.3 The CATE simulation workflow is substantially packaged

The package includes a real simulation and replication layer for the CATE setting, including:

- simulation DGPs,
- benchmark methods,
- `sim_cate_run()`,
- per-method performance summaries,
- and an installed replication script under `inst/paper-replication/cate/run.R`.

This is a strong sign that the package is no longer just a low-level methods sandbox. At least for the CATE simulation application, the repo now contains a real paper-facing workflow.

---

## 1.4 Diagnostics and reporting are now materially in place

The earlier version of this note treated diagnostics as the main missing layer.
That is no longer accurate.

The current repo now exposes:

- a public `diagnostics()` generic with methods for `np_result`, `cate_result`, and `plqr_result`,
- weight-concentration and distance-scale diagnostics,
- CATE overlap diagnostics,
- PLQR residual-Gram conditioning diagnostics,
- compact tuning metadata in fitted-object diagnostics, and
- dedicated diagnostics/reporting documentation in the README and vignettes.

So the package no longer stops at "estimate plus standard error."
It now has a real fitted-object diagnostics surface.

---

## 1.5 The modular architecture is now strong

The package has a modular decomposition across:

- distance / ranking,
- rank weights,
- tuning,
- scores,
- wrappers.

This is good architecture and should be preserved.

In particular, the score API is already flexible enough that it could serve as the basis of a more general conditional Z-estimation front end later.

That said, the current modularity is stronger than the current high-level abstraction. In other words:

- the building blocks are ahead of the final user-facing framework.

---

# 2. What is still missing

The remaining gaps fall into two different categories and should not be mixed up.

## 2.1 Gaps relative to the *paper applications*

These are the things that still make the package feel incomplete even for the concrete workflows you already show.

## 2.2 Gaps relative to the *full conditional Z-estimation framework*

These are broader framework gaps. They matter if the package is supposed to embody the full estimator class described in the paper, not just the currently implemented wrappers.

---

# 3. Main remaining gap for the paper-facing version: diagnostics (historical heading; see status update below)

Status update, 2026-04-26:
the public diagnostics/reporting layer has now landed. The live simulation-facing gap is no longer "add diagnostics at all." It is the higher-level tuning-sensitivity / stability workflow that would rerun a fit over the candidate `\kappa` grid and summarize estimate / SE / CI / diagnostics in one object.

Read the detailed discussion in this section as the design history for why that diagnostics layer was important, not as the current top-priority task list.

At the time this note was first drafted, this was the most important practical gap.

The package now estimates the main objects, computes jackknife standard errors, and exposes some helpful ingredients such as:

- weights,
- ranks,
- distances,
- effective sample size,
- neighbors,
- tuning payloads.

At that point there was **not yet a real diagnostics/reporting layer** that turned these ingredients into paper-relevant diagnostics.

That was the single most important thing to build next.

## Why this mattered at the time

At that point the package could often produce the estimate and standard error, but it still did not fully help the user answer questions like:

- Is the fit meaningfully local?
- Is the local neighborhood too concentrated on a handful of observations?
- Is local overlap poor?
- Is the local moment system nearly singular?
- How sensitive is the result to reasonable tuning variation?

Those are exactly the questions a serious empirical user, a referee, or even your future self will care about.

So this is not “nice to have.” It is the missing paper-facing layer that sits between raw computation and serious applied use.

---

# 4. Secondary paper-facing gap: sensitivity-grid workflow

This section is now the live package-side follow-up for simulation work. The manuscript repo already owns the baseline and fixed-`\kappa` run orchestration; what remains package-side is a reusable sensitivity/reporting helper rather than new baseline-simulation infrastructure.

The tuning infrastructure now supports candidate grids, but there is not yet a high-level workflow that uses those grids to evaluate stability of the estimate and its inferential objects.

In other words, the package now has the machinery for tuning metadata, but not yet the paper-facing workflow:

- run estimator on the default paper tuning grid,
- summarize estimate / SE / CI / diagnostics across grid points,
- report whether the result is stable.

This is the natural next step after diagnostics.

It should probably not be built as a separate ad hoc script. It should be built into the diagnostics/reporting layer.

---

# 5. Main framework gap: no general conditional Z-estimation front end yet

This is the biggest “framework completeness” gap.

The package now has a substantial score API, including parameter-aware evaluation and gradient support. That means the modular pieces needed for general moment-based estimation are starting to exist.

However, there is still no single polished front end that does something like:

1. define a localized moment / score,
2. supply or fit nuisances,
3. build localized estimating equations,
4. solve the local Z-equation,
5. return estimate, standard error, Jacobian/bread information, and diagnostics.

Instead, the current state is still closer to:

- there are strong wrappers for specific use cases,
- and there are modular pieces from which a more general framework could be built.

That is good progress, but it is not the same as full framework completeness.

So if the package is supposed to *fully represent the paper’s general estimator class*, that work still remains.

---

# 6. Usability gap: nuisance fitting is still largely external

This is lower priority than diagnostics, but still worth noting.

For example:

- `plqr_dnn()` expects nuisance objects like `q_tau(X)` and `m_D(X)` to be supplied.
- `cate_dnn()` is more flexible, but the user still largely provides nuisance inputs rather than calling one fully automated empirical workflow.

For the paper itself this may be totally fine. In fact, keeping nuisance fitting external may be the right design choice for now.

Still, from a package-completeness perspective, this means the package is currently strongest when used by someone who already understands the estimator and can prepare the right inputs.

That is acceptable for now. It is just not the first thing to build next.

---

# 7. Possible empirical-application gap

The package clearly includes a real installed replication workflow for the CATE simulation study.

What I did **not** see in the inspected material is an equally mature installed package-level workflow for a real empirical application, assuming such an application remains in the paper.

So if you want the package to fully support the actual paper applications in a replication sense, there may still be a gap here.

This is not as urgent as the tuning-sensitivity/reporting gap unless the application section is a major part of the paper/package story.

---

# 8. Recommended priority order

This is the implementation order I would recommend.

## Priority 1: tuning-grid sensitivity summary

This is now the most important missing piece.

Reason:

Updated rationale, 2026-04-26:
the diagnostics layer already exists, the tuning grid already exists, and the missing package-side piece is the structured stability/sensitivity workflow across that grid.

- it is the main remaining gap even for workflows already implemented,
- it improves the package’s credibility immediately,
- it gives empirical users something they can actually inspect and report,
- it naturally complements the new tuning API.

## Priority 2: wrapper summary integration

Once diagnostics exist, integrate them into user-facing summaries / print methods rather than forcing users to pull raw objects apart.

Reason:

Updated rationale, 2026-04-26:
the issue is no longer whether diagnostics can be computed, but whether they surface naturally in routine user-facing summaries.

- the tuning grid now exists,
- but there is no structured way to report stability across it,
- and that is one of the most natural paper-facing deliverables.

## Priority 3: empirical-application packaging, if needed

Only prioritize this highly if the empirical application is part of the package story.

This matters only if the empirical application is part of the package story rather than remaining manuscript-only replication code.

## Priority 4: generic conditional Z-estimation front end

This is the big framework expansion.

It matters, but it is not the first bottleneck for the current paper-facing state.

## Priority 5: nuisance-fitting convenience

This is useful, but should come later.

---

# 9. What diagnostics should be included

Historical note, 2026-04-26:
the detailed diagnostics design material below remains useful as the route record for the API that ultimately landed. It should no longer be read as the live top-priority task list. The live package-side simulation follow-up is the higher-level tuning-sensitivity/reporting layer described above.

The diagnostics should be organized in a way that matches both the theory and the applied use case.

I would distinguish between:

- **core diagnostics** that every DNN-localized fit should report,
- **application-specific diagnostics** that depend on the object being estimated.

---

## 9.1 Core diagnostics

These should be available for all main wrappers.

### (a) Effective local sample size

This is the most important diagnostic.

Use the standard quantity already present in the package:

\[
n_{\mathrm{eff}}(x) = \left(\sum_i w_i(x)^2\right)^{-1}.
\]

This should be reported prominently, not buried.

Reason:

- it is interpretable,
- it is already available,
- it connects directly to the localization and self-normalization story.

---

### (b) Weight concentration

`n_eff` alone is not enough. You can have a seemingly moderate `n_eff` but still have a fit that is effectively driven by very few observations.

So report at least:

- `max_weight`
- `top1_weight_share`
- `top5_weight_share`
- `top10_weight_share`

Potentially also:

- number of weights above simple thresholds like `1%`, `5%`, etc.

This gives the user an immediate sense of whether the local fit is overly concentrated.

---

### (c) Realized neighborhood scale

You should report some empirical neighborhood-size quantity.

The cleanest candidate is:

- distance to the `s`-th nearest observation, i.e. a local radius proxy.

Potentially also:

- weighted mean distance,
- weighted distance quantiles.

Reason:

- `s` is not directly interpretable across densities,
- the realized radius gives a geometric sense of locality at the target point.

---

### (d) Tuning payload and tuning sensitivity

For fits using the built-in paper tuning or any grid-based tuning rule, diagnostics should expose:

- the default chosen `s`,
- the candidate grid,
- feasibility flags,
- and ideally a convenience summary of the resulting fit across candidate grid values.

This is already partially available as raw metadata. The missing part is presentation and standardization.

---

## 9.2 CATE-specific diagnostics

These matter specifically for treatment-effect estimation.

### (a) Local overlap

This is probably the most important application-specific diagnostic on the CATE side.

At a minimum report:

- weighted treated mass:
  \[
  \sum_i w_i(x) D_i
  \]
- weighted control mass:
  \[
  \sum_i w_i(x) (1 - D_i)
  \]
- weighted treated share

Reason:

- even if the formal nuisance quantities look fine, poor effective overlap in the localized neighborhood is a first-order issue,
- this is directly interpretable,
- and it is easy to report.

---

### (b) Optional propensity extremeness summary

If propensity values are available, you could also report:

- weighted mean propensity,
- weighted lower and upper quantiles of estimated propensity,
- perhaps min/max after clipping.

This is useful, but secondary to treated/control weight mass.

So this should not be the headline diagnostic.

---

## 9.3 PLQR / local moment-system diagnostics

The key application-specific issue here is local identification and numerical conditioning.

### (a) Local Jacobian / bread strength

Report something like:

- smallest eigenvalue of the estimated local Jacobian or bread,
- condition number of the estimated local Jacobian,
- or, if computationally more natural, the condition number of the relevant local weighted Gram matrix.

Reason:

- a user can otherwise get a point estimate and standard error without realizing the local problem is nearly singular,
- this is exactly the kind of thing that matters for localized moment systems.

This should be one of the central diagnostics for `plqr_dnn()` or any future general conditional Z front end.

---

## 9.4 Secondary diagnostics

These are optional for a first pass.

### (a) Local covariate anisotropy / geometry

Examples:

- weighted covariance condition number,
- local Mahalanobis vs Euclidean shape diagnostics.

These are interesting but not first priority.

### (b) More elaborate localization summaries

For example:

- entropy of the weight vector,
- Gini coefficient of weights,
- alternative concentration measures.

These are potentially useful, but probably not needed in the first implementation pass.

---

# 10. Recommended diagnostics API

The diagnostics should be implemented as a real API layer, not just ad hoc extracted values.

---

## 10.1 Low-level helper

Add a low-level internal helper that computes generic diagnostics from weights and optional distance information.

For example:

- `.dnn_weight_diagnostics()`

### Inputs

- weight matrix or vector,
- optional distance matrix/vector,
- optional tuning payload,
- optional `k`, `s`,
- optional target identifiers.

### Outputs

A stable tabular object or list containing:

- `n_eff`
- `max_weight`
- top-k cumulative shares
- optional radius / distance summaries
- basic tuning information

This helper should be purely generic.

It should know nothing about treatment or Jacobians.

---

## 10.2 Wrapper-specific diagnostics methods

Then add wrapper-aware diagnostics methods, for example:

- `diagnostics(np_result, ...)`
- `diagnostics(cate_result, ...)`
- `diagnostics(plqr_result, ...)`

or, if you do not want to commit to S3 generics yet, separate helpers such as:

- `diagnose_np_dnn()`
- `diagnose_cate_dnn()`
- `diagnose_plqr_dnn()`

### For `np_result`

Use the generic weight diagnostics and tuning diagnostics.

### For `cate_result`

Add:

- treated mass,
- control mass,
- treated share,
- optional propensity summaries.

### For `plqr_result`

Add:

- local bread/Jacobian strength diagnostics,
- condition number / eigenvalue summaries.

---

## 10.3 Summary integration

After the diagnostics objects exist, integrate them into:

- `summary()` methods,
- maybe custom `print()` methods,
- possibly `as.data.frame()` for convenience.

The key point is: diagnostics should not remain only as raw return extras that the user must manually assemble.

---

# 11. Recommended tuning-sensitivity workflow

The new tuning API is already a strong base. The missing part is a standardized sensitivity workflow.

I would recommend a helper along the lines of:

- `sensitivity_dnn()`
- or `diagnose_tuning_stability()`

This helper should:

1. take an estimator call or fit object plus a tuning grid,
2. re-estimate across candidate grid values,
3. summarize:
   - estimate,
   - standard error,
   - CI length,
   - `n_eff`,
   - weight concentration,
   - overlap / Jacobian diagnostics where relevant.

This gives you a compact report of whether the result is robust across plausible localizations.

That is much more useful than exposing only the raw candidate `s` values.

---

# 12. Proposed output fields for diagnostics

Below is a concrete first-pass diagnostics schema.

## 12.1 Universal fields

- `target_id`
- `k`
- `s`
- `method`
- `distance`
- `n_eff`
- `max_weight`
- `top1_weight_share`
- `top5_weight_share`
- `top10_weight_share`
- `radius_s`
- `weighted_mean_distance`
- `weighted_median_distance`

## 12.2 Tuning fields

- `tuning_source`
- `tuning_name`
- `tuning_default_s`
- `tuning_grid_available`
- maybe `tuning_kappa_default`

## 12.3 CATE-specific fields

- `treated_mass`
- `control_mass`
- `treated_share`
- `weighted_mean_propensity` (optional)
- `weighted_propensity_q10` / `q90` (optional)

## 12.4 PLQR-specific fields

- `bread_min_eig`
- `bread_max_eig`
- `bread_condition_number`

---

# 13. Recommended first implementation pass

If you want to keep the next step manageable, I would implement diagnostics in the following staged way.

## Phase 1: generic weight diagnostics

Implement for all wrappers:

- `n_eff`
- `max_weight`
- `top1_weight_share`
- `top5_weight_share`
- `top10_weight_share`
- `radius_s` if distance information is present

This should be enough to materially improve the package.

## Phase 2: treatment-effect overlap diagnostics

Add for `cate_dnn()`:

- treated mass
- control mass
- treated share

This is probably the highest-value application-specific addition.

## Phase 3: tuning sensitivity summary

Add a helper that runs the current fit across a candidate tuning grid and returns a table of:

- estimate
- SE
- CI length
- `n_eff`
- weight concentration diagnostics

## Phase 4: PLQR identification diagnostics

Add:

- local bread/Jacobian eigenvalue or condition-number summaries.

This may take a little more design work, but it is important.

---

# 14. What should *not* be done next

These are the things I would explicitly not prioritize before diagnostics.

## 14.1 Do not immediately build a huge general front end

The score API is promising, but the package will benefit more right now from becoming a better version of the workflows it already supports than from jumping too early to a maximally abstract general solver.

## 14.2 Do not prioritize automatic nuisance fitting over diagnostics

Nuisance-fitting convenience is nice, but the package is already usable for technically informed users who can supply nuisances.

By contrast, diagnostics are currently missing even for that informed user.

## 14.3 Do not bury diagnostics in optional raw return extras only

If diagnostics are implemented, they should become part of the visible user-facing reporting layer.

---

# 15. Specific codex-facing implementation plan

Below is a direct implementation checklist that should be suitable to hand to codex.

## Task A: add generic DNN diagnostics helper

Create an internal helper that takes weights and optional distances and returns a standardized diagnostics object with:

- `n_eff`
- `max_weight`
- top-k cumulative shares
- optional radius / weighted distance summaries
- basic tuning metadata if available

Make the helper vectorized over multiple targets.

## Task B: add diagnostics methods for fit objects

Add user-facing diagnostics helpers for:

- `np_result`
- `cate_result`
- `plqr_result`

The methods should call the generic helper and then add object-specific diagnostics.

## Task C: add CATE overlap diagnostics

For `cate_result`, compute and report:

- weighted treated mass
- weighted control mass
- weighted treated share

Use the actual weight matrix returned by the fit object, or reconstruct it if needed in a stable way.

## Task D: add PLQR local identification diagnostics

For `plqr_result`, compute and report a local Jacobian / bread conditioning summary.

Prefer a simple and interpretable first-pass diagnostic such as:

- smallest eigenvalue,
- condition number.

## Task E: add tuning-stability helper

Create a helper that:

- takes a fit specification or fit object,
- evaluates the estimator across the tuning grid,
- returns a tidy table with:
  - estimate
  - SE
  - CI length
  - `n_eff`
  - concentration diagnostics
  - application-specific diagnostics if available

## Task F: integrate into summaries

Update summary / print behavior so that diagnostics are visible by default in a compact way, without overwhelming the user.

At minimum:

- show `n_eff`,
- show a concentration summary,
- show overlap or conditioning diagnostics where relevant,
- show whether tuning metadata is available.

## Task G: tests

Add tests for:

- correct `n_eff` and concentration calculations,
- stable behavior with one or multiple targets,
- correct overlap diagnostics for CATE,
- stable PLQR conditioning diagnostics,
- tuning-grid sensitivity helper,
- summary output not breaking existing object structure.

---

# 16. Longer-run framework task

After diagnostics are in place, the major framework task is to build a true general conditional Z-estimation front end on top of the score API.

That future front end should likely do something like:

- take a score / moment description,
- take nuisance inputs and possibly derivative information,
- localize with DNN,
- solve the local estimating equation,
- report estimate, standard error, and diagnostics.

But this should be treated as a separate later project, not as the immediate next patch.

---

# 17. Final assessment

The package is now in a good place structurally.

It has moved beyond being just a collection of low-level DNN components. It now contains:

- real user-facing wrappers,
- paper-default tuning integration,
- a public diagnostics/reporting API,
- CATE simulation replication support,
- modular score infrastructure.

So the package is **not far away** from being a strong paper-facing implementation.

However, the most important remaining simulation-facing gap is no longer diagnostics itself. It is the **tuning-grid sensitivity / stability layer** that would sit on top of the current diagnostics and tuning APIs.

That is the next thing to fix.

If that layer is implemented well, the package will feel much closer to complete for the actual paper applications and for the second-stage fixed-`\kappa` simulation agenda.

If you later want the package to fully embody the general framework of the paper, the next major step after that is a **general localized conditional Z-estimation front end** built on the score API.

---

## Short priority summary

1. **Add tuning-grid sensitivity summaries**
2. **Integrate diagnostics into summary/print more fully**
3. **Package empirical application workflows if needed**
4. **Later: build a full general conditional Z-estimation front end**
5. **Only after that: think about nuisance-fitting convenience**
