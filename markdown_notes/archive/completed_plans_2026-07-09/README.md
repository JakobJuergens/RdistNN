# Archived completed implementation plans (2026-07-09)

These notes are completed or superseded implementation plans, moved here
during the 2026-07-09 markdown consolidation pass so that the live
`markdown_notes/` surface only carries current material. Each file has a
status header recording when and how it was completed. They are retained
for provenance only and generally use pre-2026-07-08 API naming
(`kappa`, `q_col`, `beta`, `se`/`vcov`/`omega`, `register_rank`, ...).

The live roadmap is `markdown_notes/rdistnn_roadmap.md` (the successor
roadmap `rdistnn_next_overhaul_remaining_tasks.md` was implemented 2026-07-10
and is archived in `../completed_plans_2026-07-11/`).
The canonical design policy corpus remains `markdown_notes/design/`.

Contents:

- `rdistnn_public_release_cleanup_note.md` — pre-public breaking cleanup
  roadmap; implemented 2026-07-08/09, remaining items extracted into the
  live roadmap note.
- `plqr_maturity_parity_implementation_plan.md` — PLQR/CATE maturity
  parity plan; merged 2026-07-08 (see design note 15).
- `plotting_layer_implementation_plan.md` — plotting layer plan; landed
  as `R/plot-*.R`.
- `diagnostics_reporting_implementation_plan.md` — diagnostics layer
  plan; landed 2026-04-26, follow-up closed by design note 14.
- `iv_sim_integration_note.md` — IV simulation layer design; implemented
  as experimental `sim_iv_*`.
- `local_gmm_rdistnn_implementation_note.md` — local GMM estimator-layer
  design; implemented as experimental `gmm_dnn()`.
- `tdnn_rdistnn_implementation_plan.md` — TDNN rank-weight and fast
  jackknife plan; implemented.
