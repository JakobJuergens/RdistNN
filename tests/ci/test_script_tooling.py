"""Exercise real Ruff diagnostics in disposable independent checkouts."""

from contextlib import redirect_stderr
import importlib.util
import io
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
import unittest
from unittest.mock import patch


ROOT = Path(__file__).resolve().parents[2]


class ScriptToolingTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory(prefix="rdistnn-script-tools-")
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        (self.root / "tools/ci").mkdir(parents=True)
        (self.root / "tests").mkdir()
        for source in (".ruff.toml", "tools/ci/check_scripts.py"):
            shutil.copy2(ROOT / source, self.root / source)
        self.runner = self.root / "tools/ci/check_scripts.py"

    def run_check(self, *args):
        return subprocess.run([sys.executable, str(self.runner), *args],
                              cwd=self.root.parent, capture_output=True, text=True)

    def fixture(self, name, content):
        path = self.root / name
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(content, encoding="utf-8")
        return path

    def test_real_lint_and_selected_formatting_outside_checkout(self):
        selected = self.fixture("tools/selected.py", "value= 1\nprint(value)\n")
        other = self.fixture("tools/other.py", "other= 2\nprint(other)\n")
        self.assertEqual(self.run_check().returncode, 0)
        result = self.run_check("--format-check", "tools/selected.py")
        self.assertNotEqual(result.returncode, 0, result.stdout + result.stderr)
        self.assertEqual(self.run_check("--format", "tools/selected.py").returncode, 0)
        self.assertEqual(selected.read_text(), "value = 1\nprint(value)\n")
        self.assertEqual(other.read_text(), "other= 2\nprint(other)\n")
        self.assertNotEqual(self.run_check("--format").returncode, 0)
        selected.write_text("import math\n", encoding="utf-8")
        result = self.run_check("tools/selected.py")
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("F401", result.stdout + result.stderr)
        self.assertFalse((self.root / ".ruff_cache").exists())

    def test_excluded_sources_and_empty_selection_fail(self):
        for name in ("tools/vendor/broken.py", "tools/generated/broken.py",
                     "tools/build/broken.py", "tools/archives/broken.py",
                     "inst/paper-replication/affine/source/tools/broken.py"):
            path = self.fixture(name, "invalid python (\n")
            # Direct Ruff invocation must also respect force-exclude.
            result = subprocess.run(["ruff", "check", "--no-cache", "--config",
                                     str(self.root / ".ruff.toml"), str(path)],
                                    cwd=self.root.parent, capture_output=True, text=True)
            self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
            if name.startswith("tools/"):
                result = self.run_check(name)
                self.assertNotEqual(result.returncode, 0)
                self.assertIn("zero eligible", result.stderr)
        self.assertEqual(self.run_check().returncode, 0)
        empty = self.root / "tests/empty"
        empty.mkdir()
        self.assertNotEqual(self.run_check("tests/empty").returncode, 0)
        self.assertNotEqual(self.run_check("tools/missing.py").returncode, 0)
        self.fixture("tools/not_python.ps1", "Write-Output 'valid'\n")
        self.assertNotEqual(self.run_check("tools/not_python.ps1").returncode, 0)
        self.assertNotEqual(self.run_check("../outside.py").returncode, 0)

    def test_format_directory_selection_cannot_rewrite_files(self):
        selected = self.fixture("tools/selected.py", "value= 1\nprint(value)\n")
        other = self.fixture("tests/other.py", "other= 2\nprint(other)\n")
        original = {path: path.read_bytes() for path in (selected, other)}
        for option in ("--format", "--format-check"):
            for directory in ("tools", "tests"):
                result = self.run_check(option, "tools/selected.py", directory)
                self.assertNotEqual(result.returncode, 0)
                self.assertIn("not directories", result.stderr)
                for path, content in original.items():
                    self.assertEqual(path.read_bytes(), content)

    def test_missing_and_wrong_ruff_fail_before_analysis(self):
        spec = importlib.util.spec_from_file_location("script_tooling_runner", self.runner)
        module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(module)
        with patch.object(module.shutil, "which", return_value=None), redirect_stderr(io.StringIO()):
            with self.assertRaises(SystemExit) as failure:
                module.main([])
            self.assertEqual(failure.exception.code, 2)
        wrong = subprocess.CompletedProcess([], 0, "ruff 0.0.1\n", "")
        with patch.object(module.shutil, "which", return_value="ruff"):
            with patch.object(module.subprocess, "run", return_value=wrong), redirect_stderr(io.StringIO()):
                with self.assertRaises(SystemExit) as failure:
                    module.main([])
                self.assertEqual(failure.exception.code, 2)


if __name__ == "__main__":
    unittest.main()
