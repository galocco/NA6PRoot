#!/usr/bin/env python3
"""Embed fresh thermal dimuons into saved backgrounds and reconstruct the mixed events."""

import argparse
from concurrent.futures import ThreadPoolExecutor
import json
import os
from pathlib import Path
import shutil
import subprocess


HERE = Path(__file__).resolve().parent
REPO = HERE.parent
SIM_TREES = {
    "MCKine.root": "mckine",
    "HitsVerTel.root": "hitsVerTel",
    "HitsMuonSpecModular.root": "hitsMuonSpecModular",
}
MATCH_INPUTS = (
    "MCKine.root", "TracksVerTel.root", "TracksMuonSpec.root",
    "ClustersVerTel.root", "ClustersMuonSpec.root",
)


def call(command, directory, log, env=None):
    with log.open("w") as stream:
        try:
            subprocess.run(command, cwd=directory, env=env, stdout=stream,
                           stderr=subprocess.STDOUT, check=True)
        except subprocess.CalledProcessError as error:
            raise RuntimeError(f"Command failed ({error.returncode}), see {log}") from error


def root_macro(path, *args):
    arguments = ",".join(json.dumps(arg) for arg in args)
    return ["root", "-l", "-b", "-q", f"{path}({arguments})"]


def verify_trees(directory, trees, expected=None):
    tests = [f'TFile file{i}({json.dumps(str(directory / filename))}); '
             f'auto* t{i}=file{i}.Get<TTree>("{tree}"); '
             f'if (!t{i} || t{i}->GetEntries()!=n) gSystem->Exit({i + 2});'
             for i, (filename, tree) in enumerate(trees.items())]
    filename, tree = next(iter(trees.items()))
    number = str(expected) if expected is not None else "tree->GetEntries()"
    code = (f'TFile input({json.dumps(str(directory / filename))}); '
            f'auto* tree=input.Get<TTree>("{tree}"); if (!tree) gSystem->Exit(1); '
            f'auto n={number}; {" ".join(tests)} '
            'std::cout << "ENTRIES=" << n << std::endl; gSystem->Exit(0);')
    result = subprocess.run(["root", "-l", "-b", "-q", "-e", code], check=True,
                            stdout=subprocess.PIPE, text=True)
    return int(next(line.split("=", 1)[1] for line in result.stdout.splitlines() if line.startswith("ENTRIES=")))


def combine(inputs, output, log):
    partial = output.with_suffix(".merging.root")
    call(["hadd", "-f", str(partial), *map(str, inputs)], output.parent, log)
    partial.replace(output)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("output", type=Path)
    parser.add_argument("--background", type=Path, default=HERE / "run-5photons-1000evt" / "sim")
    parser.add_argument("--replicas", type=int, default=100, help="fresh signal events per background entry")
    parser.add_argument("--workers", type=int, default=50)
    parser.add_argument("--events", type=int, help="optional event limit for validation runs")
    args = parser.parse_args()
    if min(args.replicas, args.workers) < 1 or (args.events is not None and args.events < 1):
        parser.error("replicas, workers and the event limit must be positive")

    background = args.background.resolve()
    for filename in (*SIM_TREES, "geometry.root", "na6pLayout.ini"):
        if not (background / filename).is_file():
            parser.error(f"missing background input: {background / filename}")
    n_background = verify_trees(background, SIM_TREES)
    if n_background <= 0:
        parser.error("the background sample is empty")
    n_events = min(args.events, n_background * args.replicas) if args.events else n_background * args.replicas
    workers = min(args.workers, n_events)
    output = args.output.resolve()
    output.mkdir()
    print(f"{n_events} mixed events from {n_background} backgrounds, {workers} workers", flush=True)

    layout = background / "na6pLayout.ini"
    reco_config = output / "recoparams_optmized.ini"
    reco_config.write_text("\n".join(
        line for line in (REPO / "config/reco/recoparams_optmized.ini").read_text().splitlines()
        if not line.startswith(("vtMaxChi2VertexConstrain=", "msMaxChi2VertexConstrain="))
    ) + "\n")
    geometry = background / "geometry.root"
    vertex_macro = output / "embedThermalVertex.C"
    shutil.copy2(HERE / vertex_macro.name, vertex_macro)
    call(["root", "-l", "-b", "-q", "-e",
          f'int error=0; gROOT->LoadMacro("{vertex_macro}+", &error); gSystem->Exit(error != 0);'],
         output, output / "vertexMacro.log")
    signal_generator = (
        f'{REPO}/test/genVirtualPhotonBgEvent.C("{REPO}/data/thermal_radiation_40GeV.root",'
        '"h3MPtY_HG_40GeV",1,0.,5.,0.,6.,false)'
    )
    for executable in ("na6psim", "na6prec", "hadd"):
        if not shutil.which(executable):
            parser.error(f"{executable} was not found in PATH")

    shards = []
    first = 0
    for i in range(workers):
        count = n_events // workers + (i < n_events % workers)
        shard = output / "workers" / f"worker{i:03}"
        shard.mkdir(parents=True)
        shards.append((shard, first, count))
        first += count

    def embed(shard, start, count):
        signal = shard / "signal"
        mixed = shard / "mixed"
        signal.mkdir()
        mixed.mkdir()
        env = dict(os.environ, NA6P_EMBED_BACKGROUND=str(background / "MCKine.root"),
                   NA6P_EMBED_FIRST=str(start), NA6P_EMBED_REPLICAS=str(args.replicas))
        call(["na6psim", "--nevents", str(count), "--event-offset", str(start),
              "--rnd-seed", str(482123 + start), "--doDigitization", "false",
              "--load-ini", str(layout), "--generator", signal_generator,
              "--user-vertex", str(vertex_macro),
              "--configKeyValues", f"keyval.output_dir={signal}"], signal, shard / "sim.log", env)
        call(root_macro(HERE / "embedThermalHits.C", str(background), str(signal), str(mixed), start, args.replicas),
             mixed, shard / "embedding.log")
        verify_trees(mixed, SIM_TREES, count)
        print(f"{shard.name} embedded {count} events", flush=True)

    def run_parallel(function):
        with ThreadPoolExecutor(max_workers=workers) as pool:
            futures = [pool.submit(function, *shard) for shard in shards]
            for future in futures:
                future.result()

    run_parallel(embed)
    sim = output / "sim"
    sim.mkdir()
    for filename in SIM_TREES:
        combine([shard / "mixed" / filename for shard, _, _ in shards],
                sim / filename, output / f"merge-{filename}.log")
    (sim / "geometry.root").symlink_to(geometry)
    (sim / "na6pLayout.ini").symlink_to(layout)
    verify_trees(sim, SIM_TREES, n_events)

    reco = output / "reco"
    settings = "reco.mtMCMatching=false;reco.mtZMatching=45"
    reco_options = ["--load-ini", str(layout), "--load-recoparam", str(reco_config), "--geometry", str(geometry)]
    print("Reconstructing VT and MS tracks", flush=True)
    call(["python3", str(REPO / "scripts/na6prec_parallel.py"), "--workers", str(workers),
          "--nevents", str(n_events), "--input-dir", str(sim), *reco_options,
          "--doMatching", "false", "--configKeyValues", f"{settings};keyval.output_dir={reco}"],
         output, output / "reco.log")
    verify_trees(reco, {"TracksVerTel.root": "tracksVerTel", "TracksMuonSpec.root": "tracksMuonSpec"}, n_events)

    def match(shard, start, count):
        for threshold in (3, 4, 5):
            directory = shard / f"matches{threshold}"
            directory.mkdir()
            for filename in MATCH_INPUTS:
                (directory / filename).symlink_to(reco / filename)
            call(["na6prec", *reco_options, "--firstevent", str(start), "--lastevent", str(start + count - 1),
                  "--doHitsToRecPoints", "false", "--doDigitsToRecPoints", "false",
                  "--doTrackletVertex", "false", "--doVTTracking", "false", "--doMSTracking", "false",
                  "--configKeyValues", f"{settings};reco.mtMinVTHits={threshold};keyval.output_dir={directory}"],
                 directory, shard / f"matching{threshold}.log")
            verify_trees(directory, {"TracksMatching.root": "tracksMatching"}, count)
        print(f"{shard.name} matched {count} events at all three thresholds", flush=True)

    print("Matching VT and MS tracks at thresholds 3, 4 and 5", flush=True)
    run_parallel(match)
    for threshold in (3, 4, 5):
        directory = reco if threshold == 3 else output / f"reco_{threshold}hits"
        directory.mkdir(exist_ok=True)
        if threshold != 3:
            for filename in MATCH_INPUTS:
                (directory / filename).symlink_to(reco / filename)
        combine([shard / f"matches{threshold}/TracksMatching.root" for shard, _, _ in shards],
                directory / "TracksMatching.root", output / f"merge-matching{threshold}.log")
        shutil.copy2(shards[0][0] / f"matches{threshold}/na6pRecoParam.ini", directory / "na6pRecoParam.ini")
        verify_trees(directory, {"TracksMatching.root": "tracksMatching"}, n_events)
    print(f"Complete: {n_events} reconstructed and matched events at each threshold", flush=True)


if __name__ == "__main__":
    main()
