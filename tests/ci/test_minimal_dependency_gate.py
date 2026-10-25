"""Exercise the real minimal-library gate with tiny installed source packages."""

import io
from pathlib import Path
import shutil
import subprocess
import tarfile
import tempfile
import unittest


ROOT = Path(__file__).resolve().parents[2]
RSCRIPT = shutil.which("Rscript")


@unittest.skipUnless(RSCRIPT, "Rscript is required for the installed minimal-library gate")
class MinimalDependencyGateTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        probe = subprocess.run(
            [RSCRIPT, "--vanilla", "-e", (
                "packages <- c('Rcpp', 'testthat'); "
                "missing <- packages[!vapply(packages, requireNamespace, logical(1), quietly=TRUE)]; "
                "if (length(missing)) { cat('Missing:', missing, '\\n'); quit(status=77L) }"
            )], capture_output=True, text=True, timeout=30,
        )
        if probe.returncode == 77:
            raise unittest.SkipTest(probe.stdout.strip())
        if probe.returncode:
            raise RuntimeError(probe.stdout + probe.stderr)

    def fixture_archive(self, path, statements):
        files = {
            "DESCRIPTION": (
                "Package: RdistNN\nVersion: 0.0.0\nTitle: Minimal Gate Fixture\n"
                "Description: A disposable installed-test fixture.\nLicense: GPL-3\n"
                "Author: Test Fixture\nMaintainer: Test Fixture <test@example.invalid>\n"
                "Encoding: UTF-8\nSuggests: testthat\nConfig/testthat/edition: 3\n"
            ),
            "NAMESPACE": "export(gate_fixture)\n",
            "R/fixture.R": "gate_fixture <- function() 'installed-fixture'\n",
        }
        if statements is not None:
            files["tests/testthat.R"] = (
                "library(testthat)\nlibrary(RdistNN)\ntest_check('RdistNN')\n"
            )
            files["tests/testthat/test-fixture.R"] = "".join(
                "test_that('installed fixture " + str(index) + "', { " + statement + " })\n"
                for index, statement in enumerate(statements)
            )
        with tarfile.open(path, "w:gz") as archive:
            for name, content in files.items():
                data = content.encode("utf-8")
                entry = tarfile.TarInfo("RdistNN/" + name)
                entry.size = len(data)
                archive.addfile(entry, io.BytesIO(data))

    def test_real_gate_requires_passed_installed_assertions(self):
        passing = "expect_identical(gate_fixture(), 'installed-fixture')"
        skipped = "skip('optional fixture feature unavailable')"
        cases = (
            ("passing", [passing], True),
            ("empty", None, False),
            ("all_skipped", [skipped], False),
            ("passing_with_optional_skip", [passing, skipped], True),
            ("failing", ["expect_identical(gate_fixture(), 'wrong-value')"], False),
        )
        with tempfile.TemporaryDirectory(prefix="rdistnn minimal gate ") as directory:
            for label, statements, should_pass in cases:
                with self.subTest(case=label):
                    archive = Path(directory) / (label + ".tar.gz")
                    self.fixture_archive(archive, statements)
                    result = subprocess.run(
                        [RSCRIPT, "--vanilla", str(ROOT / "tools/ci/check_minimal_dependencies.R"),
                         str(archive)], capture_output=True, text=True, timeout=60,
                    )
                    output = result.stdout + result.stderr
                    self.assertEqual(result.returncode == 0, should_pass, output)
                    self.assertEqual("Minimal-dependency installed tests passed." in output,
                                     should_pass, output)
                    if label == "empty":
                        self.assertIn("Minimal-dependency check executed no installed tests.", output)
                    if label == "all_skipped":
                        self.assertIn("Minimal-dependency check executed no passing assertions.", output)


if __name__ == "__main__":
    unittest.main()
