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

### Accepted executable identities

The accepted simulation freezes identify source commit
`40060f155cc4d84d59ee5af119f6a2534c569bab`, all 87 executable source/metadata
hashes, an explicit `SIM_RDISTNN_SOURCE_DIR` override, and native-library SHA-256
`757e97ce54b39ce0f3f42133be7125acfe3ec6178ea11c904b0e25d06d198c94`.
The study builds its own k=1 full/deletion weight tables and calls the public
full-rank deletion engine. Its four profiles contain 2,516,400 main/control
rows, including 480,000 DNN rows copied across profiles; 10,000 independent
draw families are shared across them.

The current raw-dollar SIPP packet `rdistnn-0.1.4-20261024` records released
RdistNN 0.1.4, source `4acade538f7ac8a88fe7f837b1673ebc17daa17a`, and its
installed-file identity. Its accepted fingerprint is
`5f7203f70e52646d134089a09564a91a2aea45391df8bf1ef2dc8b93bfc836a2`.
Acceptance and localization-derivation receipts verify exact represented
distance averaging and its recomputation after each deletion, including zero
rank weights in a tied block. A stale preparatory tie label remains part of the
accepted fingerprint and is explained separately; future preparation declares
the actual policy. Current fits require explicit aligned retained distances
and policy, and incomplete state requires refitting. Persistent ties and the
practical forest application remain outside the bounded-series primitive theorem.

Historical SIPP evidence instead records installed RdistNN 0.1.0 files, including
native-library SHA-256
`ed485a99bef23ebcaaed1ac1211a9b7a33f79fe3088a682074a4b7558fd724e3`
and R lazy-load database hashes in `fit_manifest.json`. Its clean sibling
source checkout at `40060f1` does not certify the installed binary's origin.
Preserve this installed-file identity rather than substituting the simulation
pin or a current package. Historical-source rebuilds must be labeled as such
unless the executable identity also matches. Rendering accepted evidence does
not require either historical estimator binary; recomputing it does.

Production resumes require matching source, dependency, seed and configuration
identities. R-version serialization differences are reconciled only after exact
equality of settings. Audit/analysis stages require their bound completion
receipts and hashes; interrupted stages are retained and rerun. The proposal's
`not_started` field is historical; completed run manifests determine status.
Use verified retained SIPP arrays locally for dependency comparisons rather
than refitting forests solely for comparison. Those arrays inherit the data's
sharing restrictions.

### Verify and refresh

Use `inst/paper-replication/affine/manifest.json` for the exact pinned paper
commit, source identities, accepted evidence and reference-display hashes.
Run `Rscript inst/paper-replication/affine/verify.R` to verify that snapshot.
Update package-owned README and launchers before synchronization; regenerate
the manifest through the sync tool if any hashed file changes afterward.
The vendored paper commit, manuscript-content snapshot, accepted historical
executable and executing validation package have separate identities.

Refresh scientific snapshots from the package checkout with:

```sh
python3 tools/replication/sync_affine.py --paper PAPER_CHECKOUT --campaign ACCEPTED_CAMPAIGN_ROOT
```

The synchronizer selects the active SIPP packet from the parent of
`application.run_receipt` in committed paper status. Before replacing source or
configurations, it checks every receipt-bound dependency (including PDFs),
nested dependencies, run/fingerprint/policy consistency and committed input
bytes. Packet and nested-fit run IDs have distinct meanings. Simulation roots
and the historical transformed-outcome inputs remain included. A missing or
inconsistent active dependency leaves the existing bundle intact. The
synchronizer needs only the four accepted configuration files, not a simulation
run cache; byte-for-byte copies verified against the accepted proposal can
supply an external campaign layout.

The sync adapter replaces the unused older CATE bootstrap with an explicit-path
shim; all nine frozen scientific files retain their accepted hashes. Rendering
supplies the pinned source identity and uses forward slashes and literal path
prefixes for portable output directories.
Standalone renderers and their default nested `source()` calls share a private
environment for each invocation. The package-owned launcher supplies that
environment while preserving explicit `local` environments; pinned upstream
renderer bytes remain unchanged. Argument or drawing errors must not replace
caller helpers or leave a graphics device open.
Integrity corruption, interrupted-stage reuse, serialization differences,
stale audit binding, dirty-source identity and portable output paths have
regression coverage in `tests/testthat/test-paper-replication-affine.R`.

Package smoke validation, accepted production evidence and a new scientific
study are distinct activities. Preserve this distinction in documentation
and output promotion; a successful packaging check does not constitute a
new production run or empirical finding.
