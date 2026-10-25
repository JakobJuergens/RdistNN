import importlib.util
import hashlib
import json
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


class ActivePacketTests(unittest.TestCase):
    def setUp(self):
        self.directory = tempfile.TemporaryDirectory()
        self.addCleanup(self.directory.cleanup)
        root = Path(self.directory.name)
        self.paper = root / "paper"
        self.packet = self.paper / "markdown_notes/affine_scores/sipp_development/selected"
        self.packet.mkdir(parents=True)
        self.application = {"run_receipt": str((self.packet / "README.md").relative_to(self.paper)),
                            "run_id": "accepted-packet", "fit_fingerprint": "fingerprint"}
        (self.packet / "README.md").write_text("Selected by metadata, not directory date.\n")
        evidence = self.packet / "evidence"
        evidence.mkdir()
        (evidence / "raw_dollars_intervals.pdf").write_bytes(b"accepted PDF")
        for name in ("raw_dollars_concentration.pdf", "raw_dollars_localization.pdf"):
            (evidence / name).write_bytes(b"accepted PDF")
        self.documents = {
            "acceptance_verification.json": dict(self.application, status="complete", n=9915,
                localization_policy="exact_distance_block_averaging_v1", averaged_weights_reconstructed=True,
                deletions_reconstructed=True, represented_distances_verified=True),
            "localization_derivation.json": {"policy": "exact_distance_block_averaging_v1",
                "represented_distance_dimensions": [9915, 7], "original_fit_sha256": "original"},
            "fit_manifest.json": {"run_id": "paper-current", "status": "complete",
                "outputs": {"fit": {"sha256": "original"}}},
            "selection.json": {"run_id": "paper-current", "fit_fingerprint": "fingerprint",
                "purpose": "empirical", "files": {}},
            "sensitivity.json": {"fit_fingerprint": "fingerprint", "purpose": "empirical", "files": {}},
        }
        self.write_packet()
        metadata = self.paper / "Paper/metadata/affine_status.json"
        metadata.parent.mkdir(parents=True)
        metadata.write_text(json.dumps({"application": self.application}))
        subprocess.run(["git", "init", "-q", str(self.paper)], check=True)
        subprocess.run(["git", "-C", str(self.paper), "add", "."], check=True)
        subprocess.run(["git", "-C", str(self.paper), "-c", "user.name=Test", "-c",
                        "user.email=test@example.invalid", "commit", "-qm", "fixture"], check=True)
        self.package = root / "package"
        self.bundle = self.package / "inst/paper-replication/affine"
        for name in ("source/sipp/fit.R", "configs/dnn/config.rds", "evidence.zip", "manifest.json"):
            path = self.bundle / name
            path.parent.mkdir(parents=True, exist_ok=True)
            path.write_bytes(b"existing bundle: " + name.encode())
        self.before = self.bundle_bytes()

    def write_packet(self):
        for name, document in self.documents.items():
            (self.packet / "evidence" / name).write_text(json.dumps(document))
        files = {str(path.relative_to(self.packet)): hashlib.sha256(path.read_bytes()).hexdigest()
                 for path in (self.packet / "evidence").iterdir()}
        (self.packet / "receipt.json").write_text(json.dumps({
            "run_id": "accepted-packet", "fit_fingerprint": "fingerprint", "files": files}))

    def bundle_bytes(self):
        return {str(path.relative_to(self.bundle)): path.read_bytes()
                for path in self.bundle.rglob("*") if path.is_file()}

    def assert_preflight_failure(self, message):
        with self.assertRaisesRegex(SystemExit, message):
            module.sync(self.paper, Path(self.directory.name) / "campaign", self.package)
        self.assertEqual(self.bundle_bytes(), self.before)

    def complete_snapshot(self):
        """Small committed source/config/evidence graph for real sync calls."""
        sim = "Codes/Simulations/affine_scores"
        sipp = "Codes/Applications/SIPP_401k"
        source = ["Codes/Simulations/bootstrap.R", f"{sim}/fixture.R", f"{sipp}/R/helper.R"]
        source += [f"{sim}/campaign_tools/{name}" for name in (
            "production_spec.json", "measure_process.py", "audit_production.R", "analyze_production.R",
            "render_production.R", "verify_production_analysis.R", "render_manuscript_simulations.R")]
        source += [f"{sipp}/{name}" for name in (
            "01_prepare_pension.R", "02_fit_cate_profiles.R", "03_profile_diagnostics.R",
            "05_export_affine_application.R", "06_localization_sensitivity.R", "07_render_affine_manuscript.R")]
        historical = "markdown_notes/affine_scores/sipp_development/2026-10-12_overhaul/evidence/accepted_dnn_estimates.csv"
        references = ["Tables/Simulations/affine_production/simulation.tex",
                      "Tables/Applications/SIPP_401k/affine_overhaul_20261013/application.tex"]
        inputs = source + references + [historical,
            "markdown_notes/affine_scores/sipp_development/2026-10-12_code_overhaul/receipt.json",
            "markdown_notes/affine_scores/simulation_program/production_analysis/data/summary.csv"]
        for relative in inputs:
            path = self.paper / relative
            path.parent.mkdir(parents=True, exist_ok=True)
            path.write_text("Committed fixture input: " + relative + "\n")
        config = Path(self.directory.name) / "campaign/configs/dnn/config.rds"
        config.parent.mkdir(parents=True)
        config.write_bytes(b"frozen configuration bytes")
        proposal = self.paper / "markdown_notes/affine_scores/simulation_program/full_10000_20261012/provenance/proposal.json"
        proposal.parent.mkdir(parents=True)
        proposal.write_text(json.dumps({"profiles": ["dnn"],
            "source_hashes": {"fixture.R": module.sha(self.paper / sim / "fixture.R")},
            "config_hashes": {"dnn": module.sha(config)}}))
        metadata = self.paper / "Paper/metadata/affine_status.json"
        metadata.write_text(json.dumps({"application": self.application, "production": {},
            "source_snapshot_commit": "manuscript-content-snapshot", "as_of": "2026-10-25",
            "manuscript": {"affine_content_hashes": {}}, "supplement": {"source_hashes": {}}}))
        adapter = self.package / "tools/replication/affine_bootstrap.R"
        adapter.parent.mkdir(parents=True)
        adapter.write_text("# Package-owned bootstrap\n")
        (self.bundle / "README.md").write_text("Package-owned replication guide\n")
        self.before = self.bundle_bytes()
        subprocess.run(["git", "-C", str(self.paper), "add", "."], check=True)
        subprocess.run(["git", "-C", str(self.paper), "-c", "user.name=Test", "-c",
                        "user.email=test@example.invalid", "commit", "-qm", "complete fixture"], check=True)
        return historical, f"{sipp}/R/helper.R", references[1]

    def test_complete_committed_snapshot_syncs_source_and_reference_evidence(self):
        historical, helper, reference = self.complete_snapshot()
        module.sync(self.paper, Path(self.directory.name) / "campaign", self.package)
        self.assertEqual((self.bundle / "source/sipp/R/helper.R").read_bytes(), (self.paper / helper).read_bytes())
        manifest = json.loads((self.bundle / "manifest.json").read_text())
        self.assertIn(historical, manifest["evidence_files"])
        self.assertIn(reference, manifest["reference_tables"])
        self.assertNotEqual(self.bundle_bytes(), self.before)

    def test_deleted_committed_inputs_leave_existing_bundle_intact(self):
        for relative in self.complete_snapshot():
            with self.subTest(relative=relative):
                path = self.paper / relative
                data = path.read_bytes()
                path.unlink()
                try:
                    self.assert_preflight_failure("Committed replication input is missing")
                finally:
                    path.write_bytes(data)

    def test_active_selection_preserves_packet_and_nested_fit_identities(self):
        later = self.packet.parent / "2099-12-31"
        later.mkdir()
        packet, identity = module.preflight_sipp(self.paper, self.application)
        self.assertEqual(packet, self.packet)
        self.assertEqual(identity["nested_fit_run_id"], "paper-current")

    def test_absent_active_packet_leaves_bundle_intact(self):
        (self.packet / "README.md").unlink()
        self.assert_preflight_failure("Required replication input")

    def test_missing_receipt_leaves_bundle_intact(self):
        (self.packet / "receipt.json").unlink()
        self.assert_preflight_failure("Required replication input")

    def test_missing_pdf_or_json_dependency_leaves_bundle_intact(self):
        for name in ("raw_dollars_intervals.pdf", "sensitivity.json"):
            with self.subTest(name=name):
                path = self.packet / "evidence" / name
                data = path.read_bytes()
                path.unlink()
                self.assert_preflight_failure("Required replication input")
                path.write_bytes(data)

    def test_corrupt_dependency_leaves_bundle_intact(self):
        (self.packet / "evidence/raw_dollars_localization.pdf").write_bytes(b"changed PDF")
        self.assert_preflight_failure("receipt hash mismatch")

    def test_rebound_inconsistent_identity_or_policy_leaves_bundle_intact(self):
        for document, field in (("selection.json", "run_id"), ("sensitivity.json", "fit_fingerprint"),
                                ("localization_derivation.json", "policy")):
            with self.subTest(document=document, field=field):
                original = self.documents[document][field]
                self.documents[document][field] = "inconsistent"
                self.write_packet()
                self.assert_preflight_failure("inconsistent")
                self.documents[document][field] = original


if __name__ == "__main__":
    unittest.main()
