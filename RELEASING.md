# Releasing RdistNN

Supported workflows include CATE, nonparametric regression, supplied-prediction
CLATE and ordinary paired-signal ratios, with DNN weights, pointwise jackknife
uncertainty, tuning, diagnostics and plots. Ratio confidence-set inversion,
PLQR and local IV-GMM remain experimental.

The local `development` branch carries version `0.1.0.9003`, aligned with the
Stata companion. Its latest maintenance changes preserve estimator interfaces
and inference limits. Choose a tagged release version after reviewing the
scope in [NEWS.md](NEWS.md); a development-version bump does not publish a
release. Matching R/Stata versions do not imply identical feature coverage.

## Prepare and verify the final source

Run from the repository root with package dependencies and development tools
installed as described in [CONTRIBUTING.md](CONTRIBUTING.md).

1. Confirm that `DESCRIPTION`, the current `NEWS.md` entry, and the intended
   `v<version>` tag agree on the chosen version. Synchronize current README and
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
