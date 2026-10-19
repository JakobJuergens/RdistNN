import shutil
import subprocess
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]


@unittest.skipUnless(shutil.which('Rscript') and shutil.which('git'), 'Rscript and git are required')
class PackageHygieneTests(unittest.TestCase):
    def setUp(self):
        directory = tempfile.TemporaryDirectory(prefix='rdistnn-hygiene-')
        self.addCleanup(directory.cleanup)
        self.root = Path(directory.name)
        (self.root / 'tools/ci').mkdir(parents=True)
        (self.root / '.github/workflows').mkdir(parents=True)
        (self.root / 'R').mkdir()
        self.script = self.root / 'tools/ci/check_repository_hygiene.R'
        shutil.copy2(ROOT / 'tools/ci/check_repository_hygiene.R', self.script)
        shutil.copy2(ROOT / '.gitignore', self.root / '.gitignore')
        subprocess.run(['git', 'init', '--quiet', str(self.root)], check=True,
                       capture_output=True, text=True)

    def check(self):
        return subprocess.run(['Rscript', str(self.script)], capture_output=True, text=True)

    def test_candidate_files_are_checked_and_ignored_builds_are_allowed(self):
        self.assertEqual(self.check().returncode, 0)
        artifact = self.root / 'debug.pyc'
        artifact.write_bytes(b'generated')
        result = self.check()
        self.assertNotEqual(result.returncode, 0)
        self.assertIn('debug.pyc', result.stdout + result.stderr)
        artifact.unlink()
        (self.root / 'build').mkdir()
        (self.root / 'build/check.log').write_text('local verification output')
        self.assertEqual(self.check().returncode, 0)

    def test_top_level_functions_have_one_owner_but_local_closures_can_share_names(self):
        (self.root / 'R/a.R').write_text('shared_helper <- function(x) x\n')
        other = self.root / 'R/b.R'
        for assignment in ('<-', '='):
            other.write_text(f'shared_helper {assignment} function(x) x + 1\n')
            result = self.check()
            self.assertNotEqual(result.returncode, 0)
            self.assertIn('shared_helper', result.stdout + result.stderr)
            self.assertIn('R/a.R', result.stdout + result.stderr)
            self.assertIn('R/b.R', result.stdout + result.stderr)
        other.write_text('other_helper <- function(x) {\n'
                         '  shared_helper <- function(y) y\n'
                         '  shared_helper(x)\n}\n')
        self.assertEqual(self.check().returncode, 0)
