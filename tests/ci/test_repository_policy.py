import re
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
POLICY = ROOT / ".github" / "workflows" / "repository-policy.yml"
ALIGNMENT = ROOT / ".github" / "workflows" / "alignment-guardrails.yaml"
DESCRIPTION = ROOT / "DESCRIPTION"


class RepositoryPolicyTests(unittest.TestCase):
    def test_standard_policy_replaces_minimal_policy(self):
        self.assertTrue(POLICY.is_file())
        self.assertFalse((POLICY.parent / "action-policy.yml").exists())

    def test_policy_enforces_shared_core(self):
        text = POLICY.read_text(encoding="utf-8")
        self.assertIn("Reject mutable external action references", text)
        self.assertIn("Check finite workflow artifact retention", text)
        self.assertIn("Lint GitHub Actions workflows", text)
        self.assertIn("persist-credentials: false", text)
        self.assertRegex(text, r"reviewdog/action-actionlint@[0-9a-f]{40}")

    def test_every_external_action_reference_is_immutable(self):
        for workflow in sorted((ROOT / ".github" / "workflows").glob("*.y*ml")):
            for line in workflow.read_text(encoding="utf-8").splitlines():
                match = re.search(r"\buses:\s*([^\s#]+)", line)
                if match and not match.group(1).startswith(("./", "docker://")):
                    self.assertRegex(match.group(1), r"@[0-9a-fA-F]{40}$")

    def test_every_uploaded_artifact_has_finite_retention(self):
        for workflow in sorted((ROOT / ".github" / "workflows").glob("*.y*ml")):
            lines = workflow.read_text(encoding="utf-8").splitlines()
            for index, line in enumerate(lines):
                if "uses: actions/upload-artifact@" not in line:
                    continue
                uses_indent = len(line) - len(line.lstrip())
                block = []
                for following in lines[index + 1:]:
                    indent = len(following) - len(following.lstrip())
                    if following.lstrip().startswith("-") and indent <= uses_indent:
                        break
                    block.append(following)
                self.assertTrue(
                    any(re.match(r"^\s*retention-days:\s*[0-9]+\s*$", item) for item in block),
                    "{}:{} upload-artifact has no finite retention".format(workflow, index + 1),
                )

    def test_alignment_uses_declared_roxygen_version(self):
        description = DESCRIPTION.read_text(encoding="utf-8")
        match = re.search(
            r"^Config/roxygen2/version:\s*(\S+)\s*$", description, re.MULTILINE
        )
        self.assertIsNotNone(match, "DESCRIPTION must declare a roxygen2 version")

        alignment = ALIGNMENT.read_text(encoding="utf-8")
        self.assertIn("roxygen2@{}".format(match.group(1)), alignment)
