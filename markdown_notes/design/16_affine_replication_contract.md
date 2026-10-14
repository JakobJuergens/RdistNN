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

## October 14, 2026 refresh verification

The implemented snapshot is paper commit
`683bf69b2f9d096b90011739b32ec3e6e16b0168`. All current numerical/SIPP source
inputs and 21 reference table fragments match that checkout. All nine frozen
affine scientific files retain the accepted production hashes; the numerical
package engine remains unchanged from production pin `40060f1`.

Installed launchers passed four-profile validation, the measured native runner
with two workers, SIPP fitting/diagnostics/export (21 estimates and 168
localization rows), and regeneration/checking of all 21 tables and six PDFs.
Output paths containing spaces, brackets and parentheses passed. All six PDFs
passed qpdf structural checks. Integrity corruption, interrupted-stage reuse,
RDS serialization differences, stale audit binding and dirty-source identity
have regression coverage. Lint, generated-file, dependency/reference and
workflow checks passed.

The final offline R CMD check used installed dependency metadata and finished
with zero errors, zero warnings and one installed-size NOTE (11.1 MB, including
the compressed accepted evidence). All vignettes and package tests passed.
The 10,000-draw production campaign and full-population SIPP fit were not rerun
for this packaging change.
