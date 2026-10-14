import importlib.util
from pathlib import Path
import subprocess
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[2]
spec = importlib.util.spec_from_file_location("affine_sync", ROOT / "tools/replication/sync_affine.py")
module = importlib.util.module_from_spec(spec)
spec.loader.exec_module(module)


class SnapshotIdentityTests(unittest.TestCase):
    def test_dirty_or_uncommitted_input_cannot_claim_head_identity(self):
        with tempfile.TemporaryDirectory() as directory:
            repo = Path(directory)
            subprocess.run(["git", "init", "-q", str(repo)], check=True)
            path = repo / "study.R"
            path.write_text("current <- TRUE\n")
            subprocess.run(["git", "-C", str(repo), "add", "study.R"], check=True)
            subprocess.run(["git", "-C", str(repo), "-c", "user.name=Test", "-c",
                            "user.email=test@example.invalid", "commit", "-qm", "fixture"], check=True)
            module.verify_snapshot(repo, ["study.R"])
            path.write_text("current <- FALSE\n")
            with self.assertRaisesRegex(SystemExit, "committed paper snapshot"):
                module.verify_snapshot(repo, ["study.R"])
            (repo / "new.R").write_text("new <- TRUE\n")
            with self.assertRaisesRegex(SystemExit, "committed paper snapshot"):
                module.verify_snapshot(repo, ["new.R"])


if __name__ == "__main__":
    unittest.main()
