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
        (self.root / 'DESCRIPTION').write_text('Package: RdistNN\nVersion: 0.1.3.9001\n')
        (self.root / 'NEWS.md').write_text('# RdistNN 0.1.3.9001\n\n# RdistNN 0.1.3\n')
        for name in ('README.Rmd', 'README.md'):
            (self.root / name).write_text('Version **0.1.3.9001** (development).\n')
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

    def test_current_package_version_owners_must_agree(self):
        for owner in ('DESCRIPTION', 'NEWS.md', 'README.Rmd', 'README.md'):
            with self.subTest(owner=owner):
                path = self.root / owner
                original = path.read_text()
                try:
                    path.write_text(original.replace('0.1.3.9001', '0.1.3.9000', 1))
                    result = self.check()
                    self.assertNotEqual(result.returncode, 0, result.stdout + result.stderr)
                    self.assertIn('version', (result.stdout + result.stderr).lower())
                finally:
                    path.write_text(original)

    def test_current_news_heading_accepts_undated_and_iso_dated_versions(self):
        for version in ('0.1.3.9001', '0.1.4'):
            (self.root / 'DESCRIPTION').write_text(f'Package: RdistNN\nVersion: {version}\n')
            for name in ('README.Rmd', 'README.md'):
                (self.root / name).write_text(f'Version **{version}**.\n')
            for suffix in ('', ' (2026-10-23)', ' (2024-02-29)'):
                with self.subTest(version=version, suffix=suffix):
                    (self.root / 'NEWS.md').write_text(
                        f'# RdistNN {version}{suffix}\n\n# RdistNN 0.1.3 (2026-10-22)\n')
                    result = self.check()
                    self.assertEqual(result.returncode, 0, result.stdout + result.stderr)

    def test_version_declarations_are_unambiguous_and_present(self):
        cases = (
            ('DESCRIPTION', 'Package: RdistNN\nVersion: 0.1.3.9001\nVersion: invalid\n'),
            ('DESCRIPTION', 'Package: RdistNN\nVersion: invalid\n'),
            ('NEWS.md', '# Upcoming\n\n# RdistNN 0.1.3.9001\n'),
            ('NEWS.md', '# RdistNN 0.1.3.9000 (2026-10-23)\n'),
            ('NEWS.md', '# RdistNN 0.1.3.9001.1 (2026-10-23)\n'),
            ('NEWS.md', '# RdistNN 0.1.3.9001 (2026/10/23)\n'),
            ('NEWS.md', '# RdistNN 0.1.3.9001 (2026-1-23)\n'),
            ('NEWS.md', '# RdistNN 0.1.3.9001 (2026-02-30)\n'),
            ('NEWS.md', '# RdistNN 0.1.3.9001 (2026-10-23) extra\n'),
            ('README.md', 'Version **0.1.3.9001**\nVersion **0.1.3.9000**\n'),
            ('README.Rmd', 'See DESCRIPTION for the version.\n'),
        )
        for owner, content in cases:
            with self.subTest(owner=owner, content=content):
                path = self.root / owner
                original = path.read_text()
                try:
                    path.write_text(content)
                    result = self.check()
                    self.assertNotEqual(result.returncode, 0, result.stdout + result.stderr)
                    self.assertIn(owner, result.stdout + result.stderr)
                finally:
                    path.write_text(original)

    def test_version_hygiene_preserves_history_pins_and_policy_examples(self):
        with (self.root / 'NEWS.md').open('a') as news:
            news.write('\nPrevious production version 0.1.2; old version 0.1.1.9000.\n')
        for owner in ('README.Rmd', 'README.md'):
            with (self.root / owner).open('a') as readme:
                readme.write('Install the pinned production release RdistNN@v0.1.3.\n')
        (self.root / 'RELEASING.md').write_text('Example: production0.1.2 starts development0.1.2.9000.\n')
        self.assertEqual(self.check().returncode, 0)
