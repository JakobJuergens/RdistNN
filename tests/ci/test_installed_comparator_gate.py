"""Exercise the package-path and no-skip guards from the actual CI step."""
import os
from pathlib import Path
import re
import shutil
import subprocess
import sys
import tempfile
import unittest


ROOT = Path(__file__).resolve().parents[2]
WORKFLOW = ROOT / ".github/workflows/R-CMD-check.yaml"


class InstalledComparatorGateTests(unittest.TestCase):
    def gate_script(self, delimiter):
        text = WORKFLOW.read_text(encoding="utf-8")
        marker = "      - name: Test installed dependency-comparison controls\n"
        self.assertIn(marker, text, "installed comparator gate is missing")
        step = text.split(marker, 1)[1].split("\n      - ", 1)[0]
        match = re.search(
            r"<<'" + delimiter + r"'\n(.*?)^          " + delimiter + r"$",
            step, re.M | re.S,
        )
        self.assertIsNotNone(match, "gate script is missing: " + delimiter)
        return "\n".join(line[10:] for line in match.group(1).splitlines())

    def test_comparator_runner_rejects_skips_failures_and_empty_suites(self):
        script = self.gate_script("PYTHON")
        with tempfile.TemporaryDirectory(prefix="rdistnn-comparator-gate-") as directory:
            root = Path(directory)
            (root / "tests/ci").mkdir(parents=True)
            (root / "tests/__init__.py").touch()
            (root / "tests/ci/__init__.py").touch()
            module = root / "tests/ci/test_affine_dependency_comparison.py"
            for statement, should_pass in (
                ("self.assertEqual(2 + 2, 4)", True),
                ("self.skipTest('installed package unavailable')", False),
                ("self.fail('control failed')", False),
                (None, False),
            ):
                with self.subTest(statement=statement):
                    source = "import unittest\n"
                    if statement is not None:
                        source += ("class Control(unittest.TestCase):\n"
                                   "    def test_control(self):\n        " + statement + "\n")
                    module.write_text(source)
                    shutil.rmtree(module.parent / "__pycache__", ignore_errors=True)
                    result = subprocess.run(
                        [sys.executable, "-c", script], cwd=root,
                        capture_output=True, text=True, timeout=10,
                    )
                    self.assertEqual(result.returncode == 0, should_pass,
                                     result.stdout + result.stderr)

    @unittest.skipUnless(shutil.which("Rscript"), "Rscript is required")
    def test_package_probe_accepts_only_the_checked_library(self):
        script = self.gate_script("RSCRIPT")
        with tempfile.TemporaryDirectory(prefix="rdistnn-check-library-") as directory:
            root = Path(directory)
            checked, other = root / "checked", root / "other"
            for library in (checked, other):
                package = library / "RdistNN"
                package.mkdir(parents=True)
                (package / "DESCRIPTION").write_text(
                    "Package: RdistNN\nVersion: 0.0.0\nTitle: Path Fixture\n"
                    "Description: A package search-path fixture.\nLicense: GPL-3\n"
                )
            for library, should_pass in ((checked, True), (other, False)):
                with self.subTest(resolved_library=library.name):
                    env = os.environ.copy()
                    env.update(R_LIBS=str(library), RDISTNN_CHECK_LIBRARY=str(checked))
                    result = subprocess.run(
                        [shutil.which("Rscript"), "--vanilla", "-"], input=script,
                        env=env, capture_output=True, text=True, timeout=10,
                    )
                    self.assertEqual(result.returncode == 0, should_pass,
                                     result.stdout + result.stderr)


if __name__ == "__main__":
    unittest.main()
