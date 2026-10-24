"""Coverage runner contracts; integration controls opt in to native compilation.

RDISTNN_COVERAGE_INTEGRATION=1 python3 -m unittest discover -s tests/ci \
    -p test_coverage_tooling.py
"""

import hashlib
import json
import os
import shutil
import subprocess
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
SCRIPT = ROOT / "tools/ci/run_coverage.R"


class CoverageToolingTests(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory(prefix="rdistnn-coverage-contract-")
        self.addCleanup(self.temporary.cleanup)
        self.base = Path(self.temporary.name)
        self.root = self.base / "checkout"
        (self.root / "tools/ci").mkdir(parents=True)
        (self.root / "tests/testthat").mkdir(parents=True)
        (self.root / "R").mkdir()
        (self.root / "src").mkdir()
        shutil.copy2(SCRIPT, self.root / "tools/ci/run_coverage.R")
        (self.root / "DESCRIPTION").write_text(
            "Package: RdistNN\nVersion: 0.0.0.9000\nTitle: Coverage Control\n"
            "Description: A small installed coverage control package.\n"
            "Authors@R: person('Test', 'Control', email='test@example.org', role=c('aut','cre'))\n"
            "License: MIT\nEncoding: UTF-8\nSuggests: testthat (>= 3.2.0)\n"
            "Config/testthat/edition: 3\n"
        )
        (self.root / "NAMESPACE").write_text("useDynLib(RdistNN)\nexport(control)\n")
        (self.root / "R/control.R").write_text(
            "control <- function(x) {\n  .Call('coverage_value', as.integer(x), PACKAGE='RdistNN')\n}\n"
        )
        (self.root / "src/control.cpp").write_text(
            '#include <R.h>\n#include <Rinternals.h>\n'
            'extern "C" SEXP coverage_value(SEXP x) {\n'
            '  return Rf_ScalarInteger(INTEGER(x)[0] + 1);\n}\n'
        )
        (self.root / "tests/testthat.R").write_text(
            'library(testthat)\nlibrary(RdistNN)\ntest_check("RdistNN")\n'
        )
        (self.root / "tests/testthat/test-control.R").write_text(
            'test_that("instrumented package computes", {\n'
            '  expect_equal(control(4), 5L)\n})\n'
        )

    def run_cli(self, *args):
        return subprocess.run(
            ["Rscript", str(self.root / "tools/ci/run_coverage.R"), *map(str, args)],
            cwd=self.base, capture_output=True, text=True,
        )

    def source_command(self, expression):
        script = self.root / "tools/ci/run_coverage.R"
        return subprocess.run(
            ["Rscript", "-e", f"source({json.dumps(str(script))}); {expression}"],
            cwd=self.base, capture_output=True, text=True,
        )

    def test_invalid_empty_and_nonmatching_filters_fail_without_outputs(self):
        output = self.base / "report"
        for value, message in (("[", "Invalid --filter"), ("", "nonempty value"),
                               ("^no-such-suite$", "selected no test files")):
            with self.subTest(value=value):
                result = self.run_cli("--filter", value, "--output-dir", output)
                self.assertNotEqual(result.returncode, 0, result.stdout + result.stderr)
                self.assertIn(message, result.stderr)
                self.assertFalse(output.exists())

    def test_existing_source_and_overlapping_outputs_are_rejected(self):
        marker = self.root / "R/control.R"
        before = marker.read_bytes()
        report = self.base / "report"
        cases = (
            ("--output-dir", self.root / "R"),
            ("--output-dir", self.root / "R/new-report"),
            ("--output-dir", report, "--workspace", report, "--install-path", report / "lib"),
            ("--output-dir", report, "--cobertura", report),
            ("--output-dir", report, "--cobertura", report / "coverage.rds"),
            ("--output-dir", report, "--workspace", report / "provenance.json"),
            ("--output-dir", report, "--workspace", report / "summary.json/nested"),
        )
        for args in cases:
            with self.subTest(args=args):
                result = self.run_cli(*args)
                self.assertNotEqual(result.returncode, 0, result.stdout + result.stderr)
                self.assertRegex(result.stderr, "destination|separate directories|reports must|Cobertura")
                self.assertFalse(report.exists())
                self.assertEqual(marker.read_bytes(), before)

    def test_cobertura_cannot_be_a_directory_or_ancestor_of_any_output_directory(self):
        # A dependency sentinel ensures a regression cannot begin compilation.
        with (self.root / "DESCRIPTION").open("a") as description:
            description.write("Imports: Rcpp (>= 999.0.0)\n")
        report = self.base / "report"
        cobertura = self.base / "xml/coverage.xml"
        for flag in ("--output-dir", "--workspace", "--install-path"):
            for directory in (cobertura, cobertura / "nested"):
                with self.subTest(flag=flag, directory=directory):
                    args = ["--output-dir", report, "--cobertura", cobertura, flag, directory]
                    result = self.run_cli(*args)
                    self.assertNotEqual(result.returncode, 0)
                    self.assertRegex(result.stderr, "Cobertura|reports must be outside")
                    self.assertFalse(report.exists())
                    self.assertFalse(cobertura.parent.exists())

    def test_symlink_cannot_hide_source_output(self):
        alias = self.base / "alias"
        alias.symlink_to(self.root / "R", target_is_directory=True)
        result = self.run_cli("--output-dir", alias / "new-report")
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("Unsafe coverage destination", result.stderr)
        self.assertFalse((alias / "new-report").exists())

    def test_wrong_dependency_version_fails_before_snapshot_or_compilation(self):
        with (self.root / "DESCRIPTION").open("a") as description:
            description.write("Imports: Rcpp (>= 999.0.0)\n")
        output = self.base / "report"
        result = self.run_cli("--output-dir", output)
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("requires Rcpp (>= 999.0.0)", result.stderr)
        self.assertFalse(output.exists())

    def test_selected_test_excluded_from_snapshot_fails_before_compilation(self):
        (self.root / ".Rbuildignore").write_text("^tests/testthat/test-control[.]R$\n")
        output = self.base / "report"
        result = self.run_cli("--filter", "^control$", "--output-dir", output)
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("Selected coverage tests are excluded", result.stderr)
        self.assertFalse(output.exists())

    def test_requested_missing_html_dependency_fails_before_snapshot(self):
        output = self.base / "report"
        expression = (
            "requireNamespace <- function(package, ...) package != 'DT'; "
            "Sys.which <- function(x) '/controlled/tool'; "
            f"coverage_dependencies({json.dumps(str(self.root))}, TRUE)"
        )
        result = self.source_command(expression)
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("--html requires package 'DT'", result.stderr)
        self.assertFalse(output.exists())

    def test_missing_native_coverage_tool_fails_preflight(self):
        result = self.source_command(
            "Sys.which <- function(x) if (identical(x, 'gcov')) '' else '/controlled/tool'; "
            f"coverage_dependencies({json.dumps(str(self.root))}, FALSE)"
        )
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("requires 'gcov' on PATH", result.stderr)

    def test_snapshot_copies_current_uncommitted_and_frozen_inputs_excluding_outputs(self):
        header = self.root / "inst/include/owned.hpp"
        frozen = self.root / "inst/paper-replication/affine/source/frozen.R"
        for path in (header, frozen):
            path.parent.mkdir(parents=True, exist_ok=True)
            path.write_bytes(b"unchanged scientific bytes\n")
        (self.root / "src/control.gcda").write_bytes(b"old instrumentation")
        (self.root / "src/control.o").write_bytes(b"old compiler output")
        (self.root / "archives").mkdir()
        (self.root / "archives/private.md").write_text("private output")
        (self.root / "docs").mkdir()
        (self.root / "docs/generated.html").write_text("generated documentation")
        (self.root / ".Rbuildignore").write_text("^docs$\n^tools$\n")
        snapshot = self.base / "snapshot"
        result = self.source_command(
            f"root <- {json.dumps(str(self.root))}; "
            f"coverage_snapshot(root, {json.dumps(str(snapshot))}, coverage_source_files(root))"
        )
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        for relative in ("R/control.R", "inst/include/owned.hpp",
                         "inst/paper-replication/affine/source/frozen.R"):
            self.assertEqual((snapshot / relative).read_bytes(), (self.root / relative).read_bytes())
        self.assertFalse((snapshot / "src/control.gcda").exists())
        self.assertFalse((snapshot / "src/control.o").exists())
        self.assertFalse((snapshot / "archives").exists())
        self.assertFalse((snapshot / "docs").exists())
        self.assertFalse((snapshot / "tools").exists())

    @unittest.skipUnless(os.environ.get("RDISTNN_COVERAGE_INTEGRATION") == "1",
                         "set RDISTNN_COVERAGE_INTEGRATION=1 to compile tiny installed controls")
    def test_actual_instrumented_install_reports_native_hits_and_propagates_failure(self):
        source = self.root / "src/control.cpp"
        initial = hashlib.sha256(source.read_bytes()).hexdigest()
        output = self.base / "focused"
        result = self.run_cli("--filter", "^control$", "--require-native", "--output-dir", output)
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        provenance = json.loads((output / "provenance.json").read_text())
        summary = json.loads((output / "summary.json").read_text())
        self.assertEqual(provenance["status"], "passed")
        self.assertEqual(provenance["makeflags"], "-j2")
        self.assertGreater(summary["R_hit_expressions"], 0)
        self.assertGreater(summary["native_hit_expressions"], 0)
        self.assertEqual(provenance["source_sha256"]["src/control.cpp"], initial)
        self.assertTrue(provenance["installed_dll_sha256"])
        self.assertEqual(Path(provenance["execution"]["package_path"]), output / "library/RdistNN")
        for name in ("coverage.rds", "coverage.csv", "summary.json", "cobertura.xml"):
            self.assertGreater((output / name).stat().st_size, 0)
        self.assertEqual(hashlib.sha256(source.read_bytes()).hexdigest(), initial)
        self.assertEqual(list((self.root / "src").glob("*.gc*")), [])
        (self.root / "tests/testthat/test-control.R").write_text(
            'test_that("deliberate installed failure", { expect_equal(control(4), -1L) })\n'
        )
        failure_output = self.base / "failure"
        result = self.run_cli("--filter", "^control$", "--output-dir", failure_output)
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("deliberate installed failure", result.stdout + result.stderr)
        self.assertEqual(json.loads((failure_output / "provenance.json").read_text())["status"], "failed")
        self.assertFalse((failure_output / "cobertura.xml").exists())

    @unittest.skipUnless(os.environ.get("RDISTNN_COVERAGE_INTEGRATION") == "1",
                         "set RDISTNN_COVERAGE_INTEGRATION=1 to install the pure-R negative control")
    def test_required_native_absence_fails_without_publishable_coverage_reports(self):
        shutil.rmtree(self.root / "src")
        (self.root / "NAMESPACE").write_text("export(control)\n")
        (self.root / "R/control.R").write_text("control <- function(x) {\n  as.integer(x) + 1L\n}\n")
        output = self.base / "no-native"
        result = self.run_cli("--filter", "^control$", "--require-native", "--output-dir", output)
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("no positive owned C++ source hits", result.stderr)
        provenance = json.loads((output / "provenance.json").read_text())
        self.assertEqual(provenance["status"], "failed")
        self.assertGreater(provenance["summary"]["R_hit_expressions"], 0)
        self.assertEqual(provenance["summary"]["native_hit_expressions"], 0)
        for name in ("coverage.rds", "coverage.csv", "summary.json", "cobertura.xml", "coverage.html"):
            self.assertFalse((output / name).exists(), name)


if __name__ == "__main__":
    unittest.main()
