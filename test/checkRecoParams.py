#!/usr/bin/env python3
"""Run a complete one-event reconstruction for every params/**/*.ini file."""

import argparse
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile


def run(command, directory, log_name):
    log = directory / log_name
    print(f"Running {command!r} in {directory}", flush=True)
    with log.open("w") as output:
        try:
            result = subprocess.run(command, cwd=directory, stdout=output,
                                    stderr=subprocess.STDOUT, timeout=180)
        except subprocess.TimeoutExpired as error:
            raise RuntimeError(f"Command timed out; see {log}") from error
    if result.returncode:
        print(log.read_text(errors="replace"), file=sys.stderr)
        raise RuntimeError(f"Command exited with {result.returncode}; see {log}")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("sim", "rec", "params", "generator", "check-output", "work-dir"):
        parser.add_argument(f"--{name}", required=True, type=Path)
    args = parser.parse_args()
    configs = sorted(args.params.resolve().rglob("*.ini"))
    if not configs:
        parser.error(f"No .ini files found in {args.params}")
    if not os.environ.get("NA6PROOT_ROOT"):
        parser.error("Source the installed init.sh before running this test")

    args.work_dir.mkdir(parents=True, exist_ok=True)
    workspace = Path(tempfile.mkdtemp(prefix="run-", dir=args.work_dir.resolve()))
    print(f"Test logs and outputs: {workspace}", flush=True)
    simulation = workspace / "simulation"
    simulation.mkdir()
    run([str(args.sim.resolve()), "-n", "1", "--rnd-seed", "12345",
         "--doDigitization", "false", "-g", str(args.generator.resolve())],
        simulation, "simulation.log")

    failed = []
    for index, config in enumerate(configs):
        directory = workspace / f"config-{index}"
        directory.mkdir()
        for name in ("geometry.root", "MCKine.root", "HitsVerTel.root",
                     "HitsMuonSpecModular.root"):
            (directory / name).symlink_to(simulation / name)
        try:
            run([str(args.rec.resolve()), "--load-recoparam", str(config),
                 "--firstevent", "0", "--lastevent", "0"], directory, "reconstruction.log")
            run([str(args.check_output.resolve())],
                directory, "check-output.log")
        except RuntimeError as error:
            print(f"FAIL {config}: {error}", file=sys.stderr)
            failed.append(config)
        else:
            print(f"PASS {config}", flush=True)
    if failed:
        raise RuntimeError(f"{len(failed)}/{len(configs)} parameter files failed; see {workspace}")
    shutil.rmtree(workspace)
    print(f"All {len(configs)} parameter files passed.")


if __name__ == "__main__":
    try:
        main()
    except (RuntimeError, OSError) as error:
        sys.exit(str(error))
