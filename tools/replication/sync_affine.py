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


def required_file(root, relative):
    path = root / relative
    if Path(relative).is_absolute() or not path.resolve().is_relative_to(root.resolve()) or not path.is_file():
        raise SystemExit(f"Required replication input is missing or outside its root: {relative}")
    return path


def read_json(path):
    try:
        return json.loads(path.read_text())
    except (OSError, ValueError) as error:
        raise SystemExit(f"Cannot read replication metadata {path}: {error}") from error


def preflight_sipp(paper, application):
    """Validate the selected packet and receipt chain without changing the bundle.

    Packet run IDs identify accepted evidence; nested run IDs identify its fit.
    Recorded historical source/executable hashes are provenance, not hashes of
    the scientific source being vendored today.
    """
    receipt_reference = required_file(paper, application["run_receipt"])
    packet = receipt_reference.parent
    if not packet.is_relative_to(paper / "markdown_notes/affine_scores/sipp_development"):
        raise SystemExit("Active SIPP receipt must select a paper SIPP evidence packet.")
    receipt = read_json(required_file(packet, "receipt.json"))
    files = receipt.get("files")
    required = {"evidence/" + name for name in (
        "acceptance_verification.json", "localization_derivation.json", "fit_manifest.json",
        "selection.json", "sensitivity.json", "raw_dollars_intervals.pdf",
        "raw_dollars_concentration.pdf", "raw_dollars_localization.pdf")}
    if not isinstance(files, dict) or not required.issubset(files):
        raise SystemExit("Active SIPP receipt is missing required evidence dependencies.")
    for relative, expected in files.items():
        if sha(required_file(packet, relative)) != expected:
            raise SystemExit(f"Active SIPP receipt hash mismatch: {relative}")
    evidence = packet / "evidence"
    accepted = read_json(evidence / "acceptance_verification.json")
    derivation = read_json(evidence / "localization_derivation.json")
    fit = read_json(evidence / "fit_manifest.json")
    selection = read_json(evidence / "selection.json")
    sensitivity = read_json(evidence / "sensitivity.json")
    policy = "exact_distance_block_averaging_v1"
    coherent = (
        receipt.get("run_id") == accepted.get("run_id") == application["run_id"]
        and receipt.get("fit_fingerprint") == accepted.get("fit_fingerprint")
        == selection.get("fit_fingerprint") == sensitivity.get("fit_fingerprint") == application["fit_fingerprint"]
        and bool(fit.get("run_id")) and fit["run_id"] == selection.get("run_id")
        and fit.get("status") == accepted.get("status") == "complete"
        and selection.get("purpose") == sensitivity.get("purpose") == "empirical"
        and accepted.get("n") == 9915
        and accepted.get("localization_policy") == derivation.get("policy") == policy
        and derivation.get("represented_distance_dimensions") == [9915, 7]
        and derivation.get("original_fit_sha256") == fit.get("outputs", {}).get("fit", {}).get("sha256")
        and all(accepted.get(key) is True for key in (
            "averaged_weights_reconstructed", "deletions_reconstructed", "represented_distances_verified"))
    )
    if not coherent:
        raise SystemExit("Active SIPP run, fingerprint, retained geometry or policy is inconsistent.")
    for document in (selection, sensitivity):
        for relative, expected in document["files"].items():
            if files.get("evidence/" + relative) != expected or sha(required_file(evidence, relative)) != expected:
                raise SystemExit(f"Active SIPP nested dependency hash mismatch: {relative}")
    return packet, {"packet": str(packet.relative_to(paper)), "nested_fit_run_id": fit["run_id"],
                    "receipt_sha256": sha(packet / "receipt.json"), "localization_policy": policy}


def verify_snapshot(paper, paths, commit="HEAD"):
    for relative in sorted(set(paths)):
        try:
            committed = subprocess.check_output(
                ["git", "-C", str(paper), "show", f"{commit}:{relative}"], stderr=subprocess.PIPE)
        except subprocess.CalledProcessError:
            raise SystemExit(f"Replication input is absent from the committed paper snapshot: {relative}")
        if required_file(paper, relative).read_bytes() != committed:
            raise SystemExit(f"Replication input differs from the committed paper snapshot: {relative}")


def verify_inventory(paper, paths, roots, commit="HEAD", suffix=None, recursive=True):
    """Include committed-but-deleted files in preflight, before replacing outputs.

    Present-file globs alone cannot detect a missing helper, historical input or
    reference table. Uncommitted additions still go through verify_snapshot.
    """
    roots = [Path(root).as_posix() for root in roots]
    present = {Path(path).as_posix() for path in paths}
    names = subprocess.check_output(
        ["git", "-C", str(paper), "ls-tree", "-r", "--name-only", "-z", commit, "--", *roots]
    ).decode().split("\0")
    expected = {name for name in names if name and (suffix is None or Path(name).suffix == suffix)
                and (recursive or Path(name).parent.as_posix() in roots)}
    for root in roots:
        if not any(name.startswith(root + "/") for name in expected):
            raise SystemExit(f"Required committed replication inventory is empty: {root}")
    missing = sorted(expected - present)
    if missing:
        raise SystemExit(f"Committed replication input is missing from the paper worktree: {missing[0]}")


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--paper", type=Path, required=True)
    parser.add_argument("--campaign", type=Path, required=True)
    args = parser.parse_args(argv)
    sync(args.paper.resolve(), args.campaign.resolve(), Path(__file__).resolve().parents[2])


def sync(paper, campaign, package):
    bundle = package / "inst/paper-replication/affine"
    commit = subprocess.check_output(["git", "-C", str(paper), "rev-parse", "HEAD"], text=True).strip()
    status = read_json(required_file(paper, "Paper/metadata/affine_status.json"))
    active_packet, accepted_sipp = preflight_sipp(paper, status["application"])
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
    helper_paths = [str(path.relative_to(paper)) for path in sorted((paper / sipp / "R").glob("*.R"))]
    verify_inventory(paper, helper_paths, [str(sipp / "R")], commit, suffix=".R", recursive=False)
    source_paths += helper_paths
    for name, expected in proposal["source_hashes"].items():
        if sha(paper / sim / name) != expected:
            raise SystemExit(f"Frozen affine source changed: {name}")
    for profile, expected in proposal["config_hashes"].items():
        if sha(required_file(campaign, f"configs/{profile}/config.rds")) != expected:
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
        str(active_packet.relative_to(paper)),
    ]
    evidence_paths = {str(path.relative_to(paper)): path
                      for root in roots for path in sorted((paper / root).rglob("*")) if path.is_file()}
    verify_inventory(paper, evidence_paths, roots, commit)
    reference_tables = {}
    reference_roots = ["Tables/Simulations/affine_production", "Tables/Applications/SIPP_401k/affine_overhaul_20261013"]
    reference_paths = {str(path.relative_to(paper)): path
                       for root in reference_roots for path in sorted((paper / root).glob("*.tex"))}
    verify_inventory(paper, reference_paths, reference_roots, commit, suffix=".tex", recursive=False)
    for relative, path in reference_paths.items():
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
        "accepted_sipp": accepted_sipp,
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
