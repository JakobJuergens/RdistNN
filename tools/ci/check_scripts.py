#!/usr/bin/env python3
"""Run pinned Ruff on maintained Python; never install tools or fix lint."""

import argparse
from pathlib import Path
import shutil
import subprocess
import sys


ROOT = Path(__file__).resolve().parents[2]
RUFF_VERSION = "0.17.2"


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    mode = parser.add_mutually_exclusive_group()
    mode.add_argument("--format", action="store_true", help="format selected files in place")
    mode.add_argument("--format-check", action="store_true", help="check selected formatting")
    parser.add_argument("paths", nargs="*", help="paths relative to the repository root")
    args = parser.parse_args(argv)
    if (args.format or args.format_check) and not args.paths:
        parser.error("formatting requires explicit selected paths")
    paths = []
    for item in args.paths or ["tools", "tests"]:
        path = (ROOT / item).resolve()
        try:
            relative = path.relative_to(ROOT)
        except ValueError:
            parser.error("paths must stay inside this repository: " + item)
        if not relative.parts or relative.parts[0] not in ("tools", "tests"):
            parser.error("select maintained paths under tools/ or tests/: " + item)
        if not path.exists():
            parser.error("selected path does not exist: " + item)
        if (args.format or args.format_check) and not path.is_file():
            parser.error("formatting requires selected .py files, not directories: " + item)
        if path.is_file() and path.suffix != ".py":
            parser.error("select Python .py files or directories: " + item)
        paths.append(str(path))
    ruff = shutil.which("ruff")
    hint = "Install with: python3 -m pip install ruff==" + RUFF_VERSION
    if not ruff:
        parser.error("Ruff " + RUFF_VERSION + " is required. " + hint)
    version = subprocess.run([ruff, "--version"], capture_output=True, text=True)
    if version.returncode or version.stdout.strip() != "ruff " + RUFF_VERSION:
        parser.error("Ruff " + RUFF_VERSION + " is required; got "
                     + (version.stdout + version.stderr).strip() + ". " + hint)
    common = ["--config", str(ROOT / ".ruff.toml"), "--no-cache"]
    selected = subprocess.run([ruff, "check", *common, "--show-files", *paths],
                              capture_output=True, text=True, cwd=ROOT)
    if selected.returncode:
        sys.stderr.write(selected.stdout + selected.stderr)
        return selected.returncode
    if not selected.stdout.strip():
        parser.error("zero eligible maintained Python files selected (check exclusions)")
    command = [ruff, "format" if args.format or args.format_check else "check", *common]
    if args.format_check:
        command.append("--check")
    return subprocess.run([*command, *paths], cwd=ROOT).returncode


if __name__ == "__main__":
    sys.exit(main())
