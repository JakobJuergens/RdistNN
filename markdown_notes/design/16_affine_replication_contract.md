# Affine paper replication contract

The installed `paper-replication/affine` bundle represents the current affine
submission candidate, not the independent coupled-score manuscript. The
manifest pins the paper commit, manuscript source hashes, accepted production
configuration, SIPP receipt, and every vendored source/evidence hash.

Scientific modules are copied unchanged from the paper repository. Package
adapters resolve installed paths, require explicit external output directories,
validate dependencies and hashes, and separate smoke validation from production.
The unused legacy CATE shared bootstrap is replaced by an explicit-path shim;
the manifest records this adapter and its upstream identity.
The pinned `source/` snapshot is excluded from package style linting to preserve
its byte identity. Package launchers remain linted; integrity and behavior checks
cover the bundled sources.
No package estimator defaults or exported APIs change for this study.

The bundle supplies three complete workflows: regenerate current manuscript
displays from accepted portable evidence; run all four affine production
profiles with the frozen seed schedule and measured Linux process receipts;
and fit/diagnose/export the raw-dollar SIPP application with literal full-fit,
OOB and five-fold prediction mechanisms plus the eight-order sensitivity grid.
Fresh runs record the actual executable and require explicit review before
manuscript promotion. Saved historical transformed SIPP estimates remain only
because the current supplement reports that distinct target; their deprecated
fit launchers are not bundled. Superseded CATE and IV study scripts are removed.

Refresh with `tools/replication/sync_affine.py` after an accepted paper evidence
update. It checks the frozen source and original configuration hashes before
copying, preserves evidence bytes in a deterministic compressed archive, and
updates the provenance manifest and 21 reference table hashes. Update the
replication article in the same pass if the protocol changes.

Verification covers installed entry points, checksum rejection, production
configuration/row accounting, both affine applications and all prediction
paths, SIPP saved-fit diagnostics, and regeneration of all 21 current table
fragments and six PDF figures. Do not rerun the 10,000-draw campaign merely to
validate packaging of unchanged scientific code.

## Maintenance checks

Use `inst/paper-replication/affine/manifest.json` for the exact pinned paper
commit, source identities, accepted evidence and reference-display hashes.
Run `Rscript inst/paper-replication/affine/verify.R` to verify that snapshot.
Integrity corruption, interrupted-stage reuse, serialization differences,
stale audit binding, dirty-source identity and portable output paths have
regression coverage in `tests/testthat/test-paper-replication-affine.R`.

Package smoke validation, accepted production evidence and a new scientific
study are distinct activities. Preserve this distinction in documentation
and output promotion; a successful packaging check does not constitute a
new production run or empirical finding.
