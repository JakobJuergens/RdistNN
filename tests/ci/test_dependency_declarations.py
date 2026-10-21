import shutil
import subprocess
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]


@unittest.skipUnless(shutil.which('Rscript'), 'Rscript is required')
class TestDependencyDeclarations(unittest.TestCase):
    def setUp(self):
        directory = tempfile.TemporaryDirectory(prefix='rdistnn-declarations-')
        self.addCleanup(directory.cleanup)
        self.root = Path(directory.name)
        (self.root / 'tests/testthat').mkdir(parents=True)
        (self.root / 'DESCRIPTION').write_text(
            'Package: RdistNN\nVersion: 0.0.0\nImports: stats\nSuggests: testthat\n'
        )

    def check(self, code):
        (self.root / 'tests/testthat/test-example.R').write_text(code)
        return subprocess.run(
            ['Rscript', str(ROOT / 'tools/ci/check_tests_declared_deps.R')],
            cwd=self.root, capture_output=True, text=True
        )

    def test_declared_public_and_private_namespace_references_are_accepted(self):
        result = self.check('stats::median(1:3)\ntestthat:::local_test_directory()\n')
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)

    def test_undeclared_public_namespace_reference_is_rejected(self):
        result = self.check('undeclared.pkg::example()\n')
        self.assertNotEqual(result.returncode, 0, result.stdout + result.stderr)
        self.assertIn('undeclared.pkg', result.stdout)
        self.assertIn('tests/testthat/test-example.R:1', result.stdout)

    def test_undeclared_private_namespace_reference_is_rejected(self):
        result = self.check('undeclared.pkg:::example()\n')
        self.assertNotEqual(result.returncode, 0, result.stdout + result.stderr)
        self.assertIn('undeclared.pkg', result.stdout)
        self.assertIn('tests/testthat/test-example.R:1', result.stdout)


if __name__ == '__main__':
    unittest.main()
