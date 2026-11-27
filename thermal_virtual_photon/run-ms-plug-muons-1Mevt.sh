#!/usr/bin/env bash
set -euo pipefail

# Reproduce the dedicated MS plug-muon sample with 1M events (2 muons/event).
# Usage: ./run-ms-plug-muons-1Mevt.sh [workers=40] [output_dir]
#        ./run-ms-plug-muons-1Mevt.sh --reco-only [workers=40] [run_dir]

script_dir=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd -P)
repo_dir=$(dirname -- "$script_dir")
data_root=${NA6P_DATA_DIR:-/data/galocco}
events=1000000
reco_only=false
if [[ ${1:-} == --reco-only ]]; then
  reco_only=true
  shift
fi
workers=${1:-40}
run_dir=${2:-"$data_root/run-ms-plug-muons-1000000evt"}

if (( $# > 2 )) || [[ ! $workers =~ ^[1-9][0-9]*$ ]]; then
  echo "Uso: $0 [workers=40] [directory_output]" >&2
  echo "     $0 --reco-only [workers=40] [directory_run]" >&2
  exit 2
fi
if "$reco_only"; then
  if [[ ! -f $run_dir/sim/geometry.root || ! -f $run_dir/sim/MCKine.root ||
        ! -f $run_dir/sim/HitsMuonSpecModular.root ]]; then
    echo "Simulazione incompleta o mancante: $run_dir/sim" >&2
    exit 1
  fi
elif [[ -e $run_dir ]]; then
  echo "La directory esiste già: $run_dir" >&2
  echo "Scegli una directory di output nuova; il run non sovrascrive dati." >&2
  exit 1
fi

layout="$repo_dir/config/setup/layout_short.ini"
generator="$repo_dir/test/genMSPlugMuons.C"
reco_params_source="$repo_dir/config/reco/reconstruction_params.ini"
for input in "$layout" "$generator" "$reco_params_source"; do
  if [[ ! -r $input ]]; then
    echo "File mancante o illeggibile: $input" >&2
    exit 1
  fi
done

bin_dir=${NA6P_BIN_DIR:-"$repo_dir/build"}
export PATH="$bin_dir:$PATH"
export NA6PROOT_ROOT=${NA6PROOT_ROOT:-"$repo_dir/install"}
export ROOT_INCLUDE_PATH="$NA6PROOT_ROOT/include${ROOT_INCLUDE_PATH:+:$ROOT_INCLUDE_PATH}"
for command_name in python3 hadd na6psim na6prec root; do
  if ! command -v "$command_name" >/dev/null; then
    echo "Comando mancante: $command_name. Carica l'ambiente ALICE/ROOT e NA6PRoot." >&2
    exit 1
  fi
done

if "$reco_only"; then
  run_dir=$(cd -- "$run_dir" && pwd -P)
  mkdir -p -- "$run_dir/reco"
else
  mkdir -p -- "$run_dir"
  run_dir=$(cd -- "$run_dir" && pwd -P)
  mkdir -- "$run_dir/sim" "$run_dir/reco"
fi
reco_params="$run_dir/reco_params.ini"
awk '!/^[[:space:]]*(vt|ms)MaxChi2VertexConstrain[[:space:]]*=/' \
  "$reco_params_source" > "$reco_params"
echo "Eventi: $events (2 muoni generati per evento)"
echo "Worker: $workers"
echo "Output: $run_dir"

if ! "$reco_only"; then
  PYTHONUNBUFFERED=1 python3 "$repo_dir/scripts/na6psim_parallel.py" \
    --workers "$workers" --nevents "$events" \
    --generator "$generator()" --load-ini "$layout" \
    --record-muon-plug-length \
    --doDigitization false \
    --configKeyValues "keyval.output_dir=$run_dir/sim" \
    2>&1 | tee "$run_dir/sim.log"
fi

# Reconstruct only the MS: no VT hit clustering, vertexing, or VT tracking.
PYTHONUNBUFFERED=1 python3 "$repo_dir/scripts/na6prec_parallel.py" \
  --workers "$workers" --nevents "$events" --input-dir "$run_dir/sim" \
  --load-ini "$layout" --load-recoparam "$reco_params" \
  --geometry "$run_dir/sim/geometry.root" \
  --doDigitsToRecPoints false \
  --doVTHitsToRecPoints false --doTrackletVertex false --doVTTracking false \
  --doMSTracking true --doMatching false \
  --configKeyValues "keyval.output_dir=$run_dir/reco" \
  2>&1 | tee "$run_dir/reco.log"

# The parallel wrapper merges ROOT trees but does not copy these macro inputs.
cp -- "$run_dir/sim/na6pLayout.ini" "$run_dir/reco/na6pLayout.ini"
cp -- "$reco_params" "$run_dir/reco/na6pRecoParam.ini"

# Energy bias/resolution and fitted-versus-true W path-length diagnostics.
root -l -b -q "$repo_dir/test/compareMSMuonEnergy.C(\"$run_dir\",-1,\"\",0.1)" \
  2>&1 | tee "$run_dir/energy_comparison.log"
root -l -b -q "$repo_dir/test/compareMSMuonPlugTruth.C(\"$run_dir\")" \
  2>&1 | tee "$run_dir/plug_truth.log"

echo "Run completato: $run_dir"
