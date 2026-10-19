import shutil
import subprocess
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
CHECKER = ROOT / 'tools/ci/check_generated_file_drift.R'


@unittest.skipUnless(shutil.which('Rscript'), 'Rscript is required')
class GeneratedRdDriftTests(unittest.TestCase):
    def test_stale_missing_and_orphan_generated_rd_are_rejected(self):
        with tempfile.TemporaryDirectory(prefix='rdistnn-generated-docs-') as directory:
            root = Path(directory)
            (root / 'R').mkdir()
            (root / 'DESCRIPTION').write_text(
                'Package: driftfixture\nVersion: 0.0.1\nTitle: Drift Fixture\n'
                'Description: A small fixture for generated documentation checks.\n'
                'License: GPL-3\nEncoding: UTF-8\n'
            )
            (root / 'R/answer.R').write_text(
                "#' A documented answer\n#' @return The number 42.\n#' @export\n"
                'answer <- function() 42\n'
            )
            subprocess.run(['Rscript', '-e',
                'roxygen2::roxygenise(commandArgs(TRUE)[1], roclets=c("namespace","rd"), load="source")',
                str(root)], check=True, capture_output=True, text=True)
            def check():
                return subprocess.run(['Rscript', str(CHECKER), '--package', str(root)],
                    cwd=ROOT, capture_output=True, text=True)
            self.assertEqual(check().returncode, 0)
            rd = root / 'man/answer.Rd'
            original = rd.read_bytes()
            rd.write_bytes(original + b'\nStale documentation.\n')
            result = check()
            self.assertNotEqual(result.returncode, 0)
            self.assertIn('man/answer.Rd', result.stdout + result.stderr)
            rd.unlink()
            self.assertNotEqual(check().returncode, 0)
            rd.write_bytes(original)
            (root / 'man/orphan.Rd').write_bytes(original.replace(b'answer', b'orphan'))
            result = check()
            self.assertNotEqual(result.returncode, 0)
            self.assertIn('man/orphan.Rd', result.stdout + result.stderr)
