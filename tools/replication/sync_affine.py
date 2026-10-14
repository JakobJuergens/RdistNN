#!/usr/bin/env python3
"""Vendor the current affine replication source and portable evidence.

Read the paper checkout; write only the package's affine bundle. Evidence is
compressed without changing its bytes. Run after a paper evidence update.
"""
import argparse
import hashlib
import json
from pathlib import Path
import shutil
import subprocess
import zipfile


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def verify_snapshot(paper, paths, commit="HEAD"):
    for relative in sorted(set(paths)):
        try:
            committed = subprocess.check_output(
                ["git", "-C", str(paper), "show", f"{commit}:{relative}"], stderr=subprocess.PIPE)
        except subprocess.CalledProcessError:
            raise SystemExit(f"Replication input is absent from the committed paper snapshot: {relative}")
        if (paper / relative).read_bytes() != committed:
            raise SystemExit(f"Replication input differs from the committed paper snapshot: {relative}")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--paper", type=Path, required=True)
    parser.add_argument("--campaign", type=Path, required=True)
    args = parser.parse_args()
    paper = args.paper.resolve()
    campaign = args.campaign.resolve()
    package = Path(__file__).resolve().parents[2]
    bundle = package / "inst/paper-replication/affine"
    commit = subprocess.check_output(["git", "-C", str(paper), "rev-parse", "HEAD"], text=True).strip()
    status = json.loads((paper / "Paper/metadata/affine_status.json").read_text())
    proposal_path = paper / "markdown_notes/affine_scores/simulation_program/full_10000_20261012/provenance/proposal.json"
    proposal = json.loads(proposal_path.read_text())
    source_paths = ["Codes/Simulations/bootstrap.R"]
    sim = Path("Codes/Simulations/affine_scores")
    source_paths += [str(sim / name) for name in proposal["source_hashes"]]
    source_paths += [str(sim / "campaign_tools" / name) for name in (
        "production_spec.json", "measure_process.py", "audit_production.R",
        "analyze_production.R", "render_production.R", "verify_production_analysis.R",
        "render_manuscript_simulations.R")]
    sipp = Path("Codes/Applications/SIPP_401k")
    source_paths += [str(sipp / name) for name in (
        "01_prepare_pension.R", "02_fit_cate_profiles.R", "03_profile_diagnostics.R",
        "05_export_affine_application.R", "06_localization_sensitivity.R",
        "07_render_affine_manuscript.R")]
    source_paths += [str(path.relative_to(paper)) for path in sorted((paper / sipp / "R").glob("*.R"))]
    for name, expected in proposal["source_hashes"].items():
        if sha(paper / sim / name) != expected:
            raise SystemExit(f"Frozen affine source changed: {name}")
    for profile, expected in proposal["config_hashes"].items():
        if sha(campaign / "configs" / profile / "config.rds") != expected:
            raise SystemExit(f"Accepted configuration changed: {profile}")
    files = {}
    source_mapping = {}
    # Validate all inputs before replacing any bundled files. The checks below
    # also protect SIPP/render/evidence provenance, not just the frozen simulator.
    roots = [
        "markdown_notes/affine_scores/simulation_program/full_10000_20261012",
        "markdown_notes/affine_scores/simulation_program/production_analysis/data",
        "markdown_notes/affine_scores/sipp_development/2026-10-12_code_overhaul",
        "markdown_notes/affine_scores/sipp_development/2026-10-12_overhaul/evidence",
    ]
    evidence_paths = {str(path.relative_to(paper)): path
                      for root in roots for path in sorted((paper / root).rglob("*")) if path.is_file()}
    reference_tables = {}
    for root in ("Tables/Simulations/affine_production", "Tables/Applications/SIPP_401k/affine_overhaul_20261013"):
        for path in sorted((paper / root).glob("*.tex")):
            relative = str(path.relative_to(paper))
            evidence_paths["reference/" + relative] = path
            reference_tables[relative] = sha(path)
    verify_snapshot(paper, source_paths + [str(path.relative_to(paper)) for path in evidence_paths.values()]
                    + ["Paper/metadata/affine_status.json"], commit=commit)
    for child in ("source", "configs"):
        if (bundle / child).exists():
            shutil.rmtree(bundle / child)
    for relative in source_paths:
        if relative == "Codes/Simulations/bootstrap.R":
            stored = "source/bootstrap.R"
        elif relative.startswith(str(sim) + "/"):
            stored = "source/sim/" + relative[len(str(sim)) + 1:]
        else:
            stored = "source/sipp/" + relative[len(str(sipp)) + 1:]
        destination = bundle / stored
        destination.parent.mkdir(parents=True, exist_ok=True)
        origin = package / "tools/replication/affine_bootstrap.R" if relative == "Codes/Simulations/bootstrap.R" else paper / relative
        shutil.copyfile(origin, destination)
        files[str(destination.relative_to(bundle))] = sha(destination)
        source_mapping[stored] = relative
    for profile in proposal["profiles"]:
        destination = bundle / "configs" / profile / "config.rds"
        destination.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(campaign / "configs" / profile / "config.rds", destination)
        files[str(destination.relative_to(bundle))] = sha(destination)
    destination = bundle / "configs/proposal.json"
    shutil.copyfile(proposal_path, destination)
    files[str(destination.relative_to(bundle))] = sha(destination)
    archive = bundle / "evidence.zip"
    with zipfile.ZipFile(archive, "w", compression=zipfile.ZIP_DEFLATED, compresslevel=9) as output:
        for relative, path in sorted(evidence_paths.items()):
            info = zipfile.ZipInfo(relative, date_time=(1980, 1, 1, 0, 0, 0))
            info.compress_type = zipfile.ZIP_DEFLATED
            info.external_attr = 0o644 << 16
            output.writestr(info, path.read_bytes(), compresslevel=9)
    files["evidence.zip"] = sha(archive)
    for path in sorted(bundle.glob("*.R")) + [bundle / "README.md"]:
        files[path.name] = sha(path)
    snapshot = {
        "schema_version": 1,
        "paper_repository": "JakobJuergens/Cond_Z_Estimation",
        "paper_commit": commit,
        "manuscript_source_snapshot_commit": status["source_snapshot_commit"],
        "paper_variant": "affine_submission_candidate",
        "paper_title": "Orthogonal Inference for Conditional Z-Estimation without Cross-Fitting: A Distributional Nearest-Neighbor Approach",
        "snapshot_date": status["as_of"],
        "production": status["production"],
        "application": status["application"],
        "manuscript_source_hashes": {name: value for name, value in
                                     (status["manuscript"]["affine_content_hashes"] | status["supplement"]["source_hashes"]).items()
                                     if name.endswith(".tex") and "affine_scores" in name},
        "files": dict(sorted(files.items())),
        "source_mapping": source_mapping,
        "evidence_files": {name: sha(path) for name, path in sorted(evidence_paths.items())},
        "reference_tables": reference_tables,
        "scope": "Current affine simulation and SIPP replication. Historical transformed SIPP estimates are saved evidence used by the current supplement; no deprecated fit or simulation launchers are bundled.",
        "source_adapters": {
            "source/bootstrap.R": {
                "upstream_sha256": sha(paper / "Codes/Simulations/bootstrap.R"),
                "reason": "Replace unused legacy CATE orchestration with an explicit-path bootstrap; the nine frozen affine scientific files are unchanged."
            }
        },
    }
    (bundle / "manifest.json").write_text(json.dumps(snapshot, indent=2) + "\n")
    print(f"Synced {len(source_paths)} source files, four accepted configurations, {len(evidence_paths)} evidence files and {len(reference_tables)} reference tables from {commit}.")
    print(f"Compressed evidence: {archive.stat().st_size / 1024**2:.2f} MiB")


if __name__ == "__main__":
    main()
