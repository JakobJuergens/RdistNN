# Releasing RdistNN

Supported workflows include CATE, nonparametric regression, supplied-prediction
CLATE and ordinary paired-signal ratios, with DNN weights, pointwise jackknife
uncertainty, tuning, diagnostics and plots. Ratio confidence-set inversion,
PLQR and local IV-GMM remain experimental.

Version `0.1.2` entered production `master` through the reviewed squash merge
of [PR #20](https://github.com/JakobJuergens/RdistNN/pull/20) on 2026-10-22
at [`c96878562f105a22290edf989631f197d9735e1c`](https://github.com/JakobJuergens/RdistNN/commit/c96878562f105a22290edf989631f197d9735e1c).
It adds duplicate-covariate
block averaging and deletion-aware survivor weights, with a theory warning;
see [NEWS.md](NEWS.md). Production source publication does not itself create
a version tag or replace GitHub release assets.

The `0.1.1` baseline combines numerical fixes, canonical covariate tie ordering,
shared implementation modules and a current-only interface, aligned with the
Stata companion. Earlier development fit schemas and obsolete aliases require
refitting with current syntax.
The previous version `0.1.0.9004` entered production `master` through the reviewed squash
merge of [PR #16](https://github.com/JakobJuergens/RdistNN/pull/16) at
[`3394ae55383ffa4cb6bfad909730378dd8c7d4cb`](https://github.com/JakobJuergens/RdistNN/commit/3394ae55383ffa4cb6bfad909730378dd8c7d4cb),
aligned with the Stata companion. Version `0.1.0.9003` identifies the previous
production snapshot.
The exact-table recurrence, compensated reductions, numerical origins, shared
covariance and fitted whitening policy first shipped in the merged 0.1.0.9003
snapshot, not 0.1.0.9002. See [NEWS.md](NEWS.md). A verified candidate requires
the final-source checks below. Public tags and assets require separate
authorization. Retain the intended release version during candidate preparation.
Matching R/Stata versions do not imply identical feature coverage.

Review the known asymptotic-tail limitation in the
[numerical reference](https://github.com/JakobJuergens/RdistNN/blob/master/markdown_notes/design/19_weighted_accumulation.md) when
assessing a candidate. Engine agreement on fixed tables does not establish
correct construction or cross-language support parity at subnormal scales.
Check that generated `Rplots.pdf` files, including nested outputs, remain
excluded from the source package.

## Candidate identity and scientific reproduction

The next planned version-pinned tag is `v0.1.2` in both packages. The production
PRs are merged; create the tags and publish their verified assets only after
the final production commits pass the required checks and tag publication is
authorized. The dated changelog records production source publication.
After the tagged release is published, install with:

```r
pak::pak("JakobJuergens/RdistNN@v0.1.2")
```

This planned tag is not an installation source yet. Keep the ordinary README
installation on production `master` and retain historical tags, assets and
paper dependency pins. Record the paired final commits and asset checksums in
the release-set manifest below.

Commit candidates locally. Publish review branches/PRs, tags, downloadable
assets and websites only when authorized; do not push directly to `master`.
After both package commits are final, bind them in one external release-set
manifest under ignored outputs. Include tested trees, artifact filenames,
sizes and SHA-256, commands, runtime/compiler information, hosted run IDs/tested
SHAs and explicit PASS/FAIL/BLOCKED/NOT_RUN gates. Do not create circular
committed HEAD references or embed an archive's own hash inside it.

For a historical reference, install the previous immutable `0.1.0.9004` source with:

```r
pak::pak("JakobJuergens/RdistNN@3394ae55383ffa4cb6bfad909730378dd8c7d4cb")
```

The unpinned README command follows production `master`; `@v0.1.0` selects the
existing older public tag. The `0.1.0.9004` source update and pkgdown deployment
did not create a new tag or replace the `v0.1.0` release assets. Different source
snapshots can share a development version, so retain the full commit or checked
archive identity.

For an unpublished future candidate, check out its full recorded commit in a
separate local checkout and use `pak::local_install(".")`, or install its
checked tarball. Put exact candidate commands in the external handoff once
both final identities are known.

The paper-evidence renderer verifies accepted stored bytes and reconstructs 21 tables
and six figures. It does not run the candidate estimator. Simulation freezes
bind source commit `40060f155cc4d84d59ee5af119f6a2534c569bab`, a source override
and a native-library hash. SIPP binds installed files and its DLL; its recorded
sibling checkout alone does not identify the loaded package. Preserve both
historical records. Historical recomputation requires the matching dependency.
Fresh launchers use and record actual installed software; outputs require
review before promotion. Candidate agreement does not change the historical
reproduction reference. See the installed bundle README.

Use `Rscript --vanilla tools/replication/compare_paper_dependency.R --help`
for bounded validation with explicit installations, file-hashed identities and
a frozen case manifest declaring quantity-specific budgets. Fresh processes
compare fixed-table arithmetic, exact-table construction and geometry
separately. Retained SIPP arrays stay private; dependency comparison requires
neither forest refitting nor production simulation.

## Prepare and verify the final source

Run from the repository root with package dependencies and development tools
installed as described in [CONTRIBUTING.md](CONTRIBUTING.md).

1. Confirm that `DESCRIPTION` and the current `NEWS.md` entry agree. Check the
   chosen tag only for a separately authorized public release. Synchronize README and
   website draft references, while preserving historical release notes.
   Regenerate help and README from their sources:

   ```r
   devtools::document()
   devtools::build_readme()
   ```

2. Run the alignment checks in `CONTRIBUTING.md`, then build the installable
   source package with all vignettes and check that exact archive:

   ```r
   dir.create("archives", showWarnings = FALSE)
   src <- devtools::build(path = "archives", vignettes = TRUE, manual = FALSE)
   rcmdcheck::rcmdcheck(
     src, args = "--no-manual", error_on = "warning",
     check_dir = "archives/check"
   )
   ```

   Investigate all errors and warnings, and review any notes. Inspect the
   archive's file list: it should contain package code, generated help,
   tests, and built vignettes, while excluding private notes, agent settings,
   CI logs, prior archives, and repository-only publication instructions.
   `archives/` is ignored and excluded from future source builds.

   Run the full **source** suite with `devtools::test()` as well (including
   checkout-only documentation contracts), mirrored-header checks and focused
   compiled/forced-R numerical/consumer regressions from `tests/README.md`.
   Install this exact archive in a fresh library, execute installed examples,
   and run its paper-replication `verify.R` and `render_evidence.R --output-dir=NEW_DIR`.
   Verify exports, the installed bundle, accepted table reconstruction and
   unchanged installed files afterward. Record rendering separately from
   estimator recomputation.

3. Build the documentation site, including examples and articles:

   ```r
   pkgdown::build_site()
   ```

   Review `docs/index.html`, `docs/articles/getting-started.html`, and
   `docs/reference/index.html`. Verify that CATE and regression lead the
   onboarding, supported CLATE is described under its inference conditions,
   and confidence-set inversion, PLQR and IV-GMM retain experimental labels.
   The generated `docs/` directory is ignored. The canonical site URL is
   <https://JakobJuergens.github.io/RdistNN/>.

4. Review the final diff. Require passing R-CMD-check on Linux, Windows,
   macOS arm64, and macOS x86-64 for the release commit, plus the alignment
   checks. CI also covers R release, devel, and oldrel on Linux. Local checks
   do not substitute for this platform matrix. Reuse successful validation
   while the source snapshot is unchanged.

## Publish on GitHub

Create the chosen `v<version>` tag and GitHub release from the verified commit. Use the
release-scope section of `NEWS.md` for the introduction, followed by a link
to the getting-started guide and the explicit PLQR experimental status.

Attach **`RdistNN_<version>.tar.gz`**, built above, as the installable R source
asset. GitHub's automatically generated repository ZIP/tarball and the CI
source-snapshot ZIP are repository snapshots; the R source package includes
built vignettes and applies `.Rbuildignore` exclusions. Include a checksum:

```bash
cd archives
rdistnn_release_version=$(Rscript -e 'cat(read.dcf("../DESCRIPTION")[1L, "Version"])')
sha256sum "RdistNN_${rdistnn_release_version}.tar.gz" > SHA256SUMS
```

After publication, download both assets without credentials and run
`sha256sum --check SHA256SUMS` from their download directory. Install the
downloaded tarball in a fresh library and execute the getting-started examples.
For version-pinned GitHub installation, use
`pak::pak("JakobJuergens/RdistNN@v<version>")`, replacing `<version>` with the
published version; the default branch remains available
for development installations.

Users can install a downloaded source asset with the actual version in its filename:

```r
install.packages("Rcpp")  # required dependency, if not already installed
install.packages("RdistNN_<version>.tar.gz", repos = NULL, type = "source")
```

Source installation needs C++17. Optional visualization requires ggplot2
(and patchwork for dashboards); other suggested dependencies support
development, replication, or experimental workflows.

Publishing a GitHub release triggers `.github/workflows/pkgdown.yaml`, which
builds the site and deploys it to `gh-pages`. Confirm that this workflow
succeeds and GitHub Pages serves that branch before announcing the public
documentation links. A manual run of the same workflow can refresh the site
without a release. Neither local site building nor pushing package sources
alone deploys the site.

## Introduce the package on jakobjuergens.com

[WEBSITE.md](WEBSITE.md) contains a short software entry and release
announcement ready to copy into the website's existing format. Link to the
GitHub release and canonical documentation site after both are live. Keep
the package documentation on its canonical URL rather than duplicating its
reference pages on the personal site. Review both links after publication.

This file records publication steps; preparing files or building the local
site does not itself tag, publish, or deploy a release.
