#!/usr/bin/env python3
"""Linux process-group measurement; sampled RSS includes shared pages per PID."""
import argparse
import csv
import hashlib
import json
import math
import os
from pathlib import Path
import resource
import signal
import subprocess
import time


class MeasurementInterrupted(Exception):
    def __init__(self, signum):
        self.signum = signum
        super().__init__(f"measurement interrupted by signal {signum}")


def stop_group(child):
    """Terminate our isolated group, including descendants after leader exit."""
    try:
        os.killpg(child.pid, signal.SIGTERM)
    except ProcessLookupError:
        pass
    try:
        child.wait(timeout=1)
    except subprocess.TimeoutExpired:
        pass
    # A reaped leader does not imply that its worker processes have stopped.
    try:
        os.killpg(child.pid, signal.SIGKILL)
    except ProcessLookupError:
        pass
    return child.wait()


def group_sample(group, coordinator):
    result = subprocess.run(["ps", "-e", "-o", "pid=,pgid=,rss="], check=True,
                            capture_output=True, text=True)
    members = []
    for line in result.stdout.splitlines():
        pid, pgid, rss = map(int, line.split())
        if pgid != group:
            continue
        try:
            stat = Path(f"/proc/{pid}/stat").read_text().rsplit(")", 1)[1].split()
            cpu = (int(stat[11]) + int(stat[12])) / os.sysconf("SC_CLK_TCK")
        except (FileNotFoundError, ProcessLookupError):
            cpu = 0
        members.append({"pid": pid, "rss_kib": rss, "cpu_seconds": cpu,
                        "coordinator": pid == coordinator})
    return members


def measure(command, output_dir, cwd=None, env=None, interval=.2, wall_limit=3600,
            rss_limit_mib=None):
    limits = [interval, wall_limit] + ([] if rss_limit_mib is None else [rss_limit_mib])
    if any(not math.isfinite(v) or v <= 0 for v in limits):
        raise ValueError("measurement intervals and limits must be positive")
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=False)
    usage0 = resource.getrusage(resource.RUSAGE_CHILDREN)
    start = time.monotonic()
    seen = {}
    combined_peak = coordinator_peak = worker_peak = 0
    deadline = False
    memory_limit = False
    reason = None
    error = None
    child = None
    returncode = 1
    state = {"status": "running", "command": command, "cwd": str(cwd),
             "started_unix": time.time()}
    old_handlers = {}
    def interrupted(signum, frame):
        raise MeasurementInterrupted(signum)
    for signum in (signal.SIGINT, signal.SIGTERM):
        old_handlers[signum] = signal.signal(signum, interrupted)
    with (output_dir / "stdout.log").open("w") as stdout, (output_dir / "stderr.log").open("w") as stderr:
      try:
        child = subprocess.Popen(command, cwd=cwd, env=env, stdout=stdout, stderr=stderr,
                                 start_new_session=True)
        state.update(pid=child.pid, process_group=child.pid)
        (output_dir / "process.json").write_text(json.dumps(state, indent=2))
        with (output_dir / "process_samples.csv").open("w", newline="") as sample_file:
            writer = csv.writer(sample_file)
            writer.writerow(["wall_seconds", "process_count", "combined_rss_kib", "coordinator_rss_kib", "largest_worker_rss_kib"])
            while child.poll() is None:
                members = group_sample(child.pid, child.pid)
                total = sum(m["rss_kib"] for m in members)
                main = max((m["rss_kib"] for m in members if m["coordinator"]), default=0)
                worker = max((m["rss_kib"] for m in members if not m["coordinator"]), default=0)
                combined_peak = max(combined_peak, total)
                coordinator_peak = max(coordinator_peak, main)
                worker_peak = max(worker_peak, worker)
                for m in members:
                    seen[m["pid"]] = max(seen.get(m["pid"], 0), m["cpu_seconds"])
                elapsed = time.monotonic() - start
                writer.writerow([elapsed, len(members), total, main, worker])
                sample_file.flush()
                if elapsed > wall_limit:
                    deadline = True
                    reason = "wall_limit"
                    stop_group(child)
                    break
                if rss_limit_mib is not None and total / 1024 > rss_limit_mib:
                    memory_limit = True
                    reason = "rss_limit"
                    stop_group(child)
                    break
                time.sleep(interval)
        returncode = 124 if deadline else 125 if memory_limit else child.wait()
      except MeasurementInterrupted as exc:
        reason = f"signal_{exc.signum}"
        returncode = 128 + exc.signum
        if child is not None:
            stop_group(child)
      except BaseException as exc:
        reason = "exception"
        error = repr(exc)
        if child is not None:
            stop_group(child)
      finally:
        # A failed or completed coordinator may leave descendants behind.
        # This group is owned solely by this measurement invocation.
        if child is not None:
            stop_group(child)
        for signum, handler in old_handlers.items():
            signal.signal(signum, handler)
    usage1 = resource.getrusage(resource.RUSAGE_CHILDREN)
    result = {**state, "status": "interrupted" if reason and reason.startswith("signal_") else
              ("completed" if returncode == 0 else "failed"),
              "returncode": returncode, "wall_seconds": time.monotonic() - start,
              "child_returncode": child.returncode if child is not None else None,
              "termination_reason": reason, "error": error,
              "coordinator_rss_peak_mib": coordinator_peak / 1024,
              "worker_rss_peak_mib": worker_peak / 1024,
              "combined_rss_peak_mib": combined_peak / 1024,
              "sampled_group_cpu_seconds": sum(seen.values()),
              "reaped_child_cpu_seconds": usage1.ru_utime + usage1.ru_stime - usage0.ru_utime - usage0.ru_stime,
              "observed_processes": len(seen), "sample_interval_seconds": interval,
              "wall_deadline_reached": deadline,
              "rss_limit_mib": rss_limit_mib, "rss_limit_reached": memory_limit,
              "memory_definition": "sampled sum RSS; shared pages counted for each process; brief peaks may be missed",
              "helper_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
    (output_dir / "process.json").write_text(json.dumps(result, indent=2))
    return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", required=True)
    parser.add_argument("--cwd", default=os.getcwd())
    parser.add_argument("--wall-limit", type=float, default=3600)
    parser.add_argument("--rss-limit-mib", type=float)
    parser.add_argument("command", nargs=argparse.REMAINDER)
    args = parser.parse_args()
    command = args.command[1:] if args.command[:1] == ["--"] else args.command
    if not command:
        parser.error("supply a command after --")
    result = measure(command, args.output, args.cwd, os.environ.copy(), wall_limit=args.wall_limit,
                     rss_limit_mib=args.rss_limit_mib)
    print(json.dumps(result, sort_keys=True))
    raise SystemExit(result["returncode"])
