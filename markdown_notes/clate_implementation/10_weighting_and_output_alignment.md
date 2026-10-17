# CLATE weighting and output alignment

Status: implemented, verified and published on `master` in `d655282`.
All six hosted package-check platforms, lint and alignment passed.
The author's follow-up authorizes alignment with the existing CATE and
nonparametric-regression interfaces, followed by direct commit and push.
This extends the original nine-note exact, single-neighbor baseline.
The author also explicitly requests plotting capabilities analogous to CATE,
with matching user interactions. Include estimate plots and generic dispatch,
weight/overlap/influence diagnostics, tuning plots and dashboards. CLATE
overlap refers to the instrument propensity, not treatment propensity.

The common summary contract is a classed list with `coefficients`, `level`,
`has_se` and `meta`. The coefficient table uses `estimate`, `std.error`,
`statistic`, `p.value`, `conf.low` and `conf.high`, with one row per target.
CLATE retains its component uncertainty and per-target computational status
in separate summary fields. Undefined or degenerate targets remain visible.
Raw fitted objects retain `theta` and paired component covariance.

Weighting parity includes exact and asymptotic single-scale DNN weights,
general `1 <= k <= s`, the existing normalization control and distance
specifications. The author's subsequent clarification excludes two-scale
TDNN from CATE and CLATE because its use there is not theoretically
justified. No TDNN ratio/CLATE entry point is added; the existing NP TDNN
interface remains separate. Exact DNN and k=1 remain the defaults. Prediction
and diagnostics must preserve the selected DNN scheme and all its parameters.
The existing `tuning_mode = "grid"` sensitivity workflow also uses the
common candidate/estimate/diagnostic tables, retaining the default fit as
the primary result and reusing frozen ranks across candidates. It does not
select a candidate or establish post-selection inference.
All full and deleted calculations use the same paired signals and frozen
full-sample geometry; covariance remains on the estimator scale.

Independent verification covers finite-subset enumeration for exact DNN,
literal/polynomial beta-increment examples for asymptotic DNN,
paired deletions and covariance, reduction
to existing CATE/NP estimates, R/compiled parity, tabular interfaces,
serialization, prediction, diagnostic reconstruction and summary schemas.
The frozen affine paper replication bundle is outside the edit scope.

This software extension does not establish that the original exact-k=1
theorem applies to every weighting choice. Documentation must distinguish
computable alternatives from their additional asymptotic obligations.

Final acceptance evidence will be recorded after implementation and review.

## Shared interfaces and extension points

The public result classes describe distinct estimands. CLATE has class
`c("clate_result", "affine_ratio_result")`; it does not impersonate CATE.
S3 dispatch selects the appropriate numerical adapter. Plotting builders
consume common tidy tables rather than depending on a renamed estimate field.

| Responsibility | Shared implementation | Method-specific responsibility |
| --- | --- | --- |
| Rank weights | `.dnn_rank_weight_spec()` in `R/user-dnn-core.R` | Selected exact/asymptotic rule, k, s, normalization |
| Orders and grids | `.dnn_resolve_s()`, `.dnn_grid_eval()` | Estimate/SE callbacks; ratio target statuses |
| Paired ratios | `.affine_ratio_evaluate()` | CLATE constructs IV signals and component-role metadata |
| Summaries | CATE-compatible coefficient schema | CLATE component/status tables and rowwise Wald availability |
| Tidy estimates | `plot_data_estimates()` S3 generic, common frame builders | Estimate field/scale and unavailable bounds |
| Candidate views | `candidate_fit()` / `candidate_fits()` S3 generics, common metadata/extras | Paired components and covariance reconstructed from retained signals |
| Estimate/weight plots | `R/plot-estimates.R`, `R/plot-weights.R` | Tidy extractor; zero line for effect estimates |
| Overlap | `R/plot-overlap.R` | CATE treatment versus CLATE instrument and propensity |
| Influence | `R/plot-influence.R` | CATE residual versus CLATE paired score residual |
| Tuning/dashboard | `plot_tuning()` / `plot_dashboard()` S3 generics, shared panel builders | Frozen-state grid adapter and diagnostics |

Numerical equality is not the only contract. Matching arguments, candidate
selection grammar, return shapes, optional-dependency guards, captions,
serialization and predictions are all tested. The common six-column summary
schema does not require unavailable ratio uncertainty to mimic CATE's
zero-variance behavior. Components and statuses remain accessible separately.
CLATE influence plots display `w * (phi_y - theta * phi_d)`; they do not
claim to show deletion pseudo-values or observation-level causal effects.

To extend an estimator:

1. Retain its scientifically distinct class and canonical numeric state.
   Document each quantity's scale before adding a common adapter.
2. Implement the existing S3 extractor/reporting methods and common schemas;
   use shared validators, rank weights and frame builders. Add a public generic
   only if no existing interface expresses the operation.
3. Supply a numerical callback for grid evaluation, retaining frozen state and
   separate candidate/target statuses. Rebuild every candidate-dependent field
   in candidate views; do not copy stale covariance or deletion arrays.
4. Reuse panel builders. Keep different scores, uncertainty calculations and
   labels in explicit adapters. Do not change an object's class to route it
   through another estimator's implementation.
5. Test independent analytic fixtures, special-case reduction, all supported
   weight controls, candidate/scalar equality, missing uncertainty, tabular
   inputs, serialization, frozen geometry and optional dependencies. Exercise
   actual ggplot builds as well as tidy data. Regenerate help and run executable
   guides after documenting the same calls users will make.

This is an incremental alignment rule, not a mandate for one universal fit
schema. PLQR moment evaluation, solved GMM coefficients, set-valued ratio
inversion and NP-only two-scale weights retain their distinct meanings.

## Execution and acceptance

- New summary, weighting, grid, candidate and plotting regressions were observed
  failing before implementation and passing after it. Independent oracles
  include exhaustive exact general-k subsets, literal asymptotic weights and
  paired deletion arrays, and perfect-compliance CATE reduction.
- CATE, NP and ratio/candidate paths now use common rank-weight construction.
  The paired evaluator, shared grid engine and common plotting builders are
  reused rather than copied into an independent CLATE plotting system.
- Tests cover selected methods, normalization, retained model state, optional
  extras, statuses, serialization, compiled/fallback/explicit-deletion parity,
  instrument overlap, paired influence and candidate/scalar covariance equality.
- Logs: ignored `archives/clate-alignment-20261017/work/`. Base `ec987e0`.
  The final full suite has 909 blocks / 6,322 passing assertions / no failures,
  errors or warnings / one opt-in visual skip. Focused plotting tests also pass
  after the last graphics-device cleanup.
- Independent review found no critical or important defects. Three minor
  reporting issues were corrected and checked: misleading unavailable-ratio
  captions, warnings from missing interval/point geometry, and a shadowed
  candidate-table variable. Draw-time checks confirm finite zero-score points
  remain visible; undefined targets remain in the numerical tables.
- Lint, declared-dependency, reference-coverage, repository-hygiene, generated
  namespace/export drift and pinned affine-bundle checks pass. There are still
  131 public reference topics; new methods are grouped with existing generics.
- The local source-package check has zero errors and warnings, with the existing
  11.5 MB installed-size NOTE. The clean checked archive is under
  `work/release/`; 364 source/help/test/vignette files matched its checked
  snapshot. The later article-only refinement prints all five example plots
  explicitly and was rebuilt separately; numerical package inputs are unchanged.
- Full pkgdown build and final article rebuild pass. Across 18 affected pages,
  all 693 local links resolve. Forest estimates, instrument overlap, paired
  influence and candidate/tuning plots were visually inspected. Generated site
  output and verification archives remain ignored; the frozen replication
  source was not edited.
- The documentation-default guard caught an introductory k=2 example. The
  executable introduction retains k=1, while the supported general-k controls
  and their independent numerical tests remain documented.
- Direct `master` publication: `d655282afc84b52ca1449b102b3bdd314421674f`.
  The [R CMD check matrix](https://github.com/JakobJuergens/RdistNN/actions/runs/39070188587)
  completed successfully on Linux devel/release/oldrel-1, Windows x64 release,
  macOS ARM release and macOS Intel release. The
  [lint](https://github.com/JakobJuergens/RdistNN/actions/runs/39070188560) and
  [alignment](https://github.com/JakobJuergens/RdistNN/actions/runs/39070188589)
  workflows also succeeded. These checks cover the final committed article,
  including all explicitly printed plots.
- The [source snapshot workflow](https://github.com/JakobJuergens/RdistNN/actions/runs/39070188584)
  succeeded at its cooldown gate; `zip-and-upload` was skipped by the 12-hour
  policy. No new hosted source archive is claimed. The local checked archive
  and verification receipts remain under the ignored work directory.
- This subsequent acceptance receipt changes only this note. Package source,
  tests, generated help, README, vignettes and site inputs remain identical to
  the verified implementation commit, so successful checks are reused. No tag,
  release, deployment, empirical campaign or Stata edit is included.
