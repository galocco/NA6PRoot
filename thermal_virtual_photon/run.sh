#!/usr/bin/env bash
set -euo pipefail

script_dir=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd -P)
repo_dir=$(dirname -- "$script_dir")

if [[ ${1:-} == -h || ${1:-} == --help ]]; then
  echo "Uso: $0 [eventi=100] [worker=4] [directory_output_nuova]"
  echo "     $0 --reco-only eventi worker directory_simulata"
  exit 0
fi

reco_only=false
if [[ ${1:-} == --reco-only ]]; then
  reco_only=true
  shift
fi
events=${1:-100}
workers=${2:-4}
if (( $# > 3 )) || [[ ! $events =~ ^[1-9][0-9]*$ || ! $workers =~ ^[1-9][0-9]*$ ]] ||
    { "$reco_only" && (( $# != 3 )); }; then
  echo "Argomenti non validi. Vedi $0 --help." >&2
  exit 2
fi

layout="$repo_dir/config/setup/layout_short.ini"
reco_params_source="$repo_dir/config/reco/recoparams_optmized.ini"
generator="$repo_dir/test/genVirtualPhotonBgEvent.C"
histogram="$repo_dir/data/thermal_radiation_40GeV.root"
photons=${THERMAL_VIRTUAL_PHOTONS:-1}
if [[ ! $photons =~ ^[1-9][0-9]*$ ]]; then
  echo "THERMAL_VIRTUAL_PHOTONS deve essere un intero positivo." >&2
  exit 2
fi
generator_call="$generator(\"$histogram\",\"h3MPtY_HG_40GeV\",${photons},0.,5.,0.,6.,true)"
for input in "$layout" "$reco_params_source" "$generator" "$histogram"; do
  if [[ ! -r $input ]]; then
    echo "File mancante o illeggibile: $input" >&2
    exit 1
  fi
done

# Prefer the local build so the simulation and reconstruction use the current code.
bin_dir=${NA6P_BIN_DIR:-"$repo_dir/build"}
export PATH="$bin_dir:$PATH"
export NA6PROOT_ROOT=${NA6PROOT_ROOT:-"$repo_dir/install"}
export ROOT_INCLUDE_PATH="$NA6PROOT_ROOT/include${ROOT_INCLUDE_PATH:+:$ROOT_INCLUDE_PATH}"

for command_name in python3 hadd na6psim na6prec; do
  if ! command -v "$command_name" > /dev/null; then
    echo "Comando mancante: $command_name. Carica l'ambiente ALICE/ROOT e NA6PRoot." >&2
    exit 1
  fi
done

if "$reco_only"; then
  run_dir=$(cd -- "$3" && pwd -P)
  if [[ ! -f $run_dir/sim/geometry.root || -e $run_dir/reco || -e $run_dir/reco_5hits ]]; then
    echo "Serve una directory con simulazione completata e senza output di ricostruzione." >&2
    exit 1
  fi
elif [[ -n ${3:-} ]]; then
  mkdir -- "$3"
  run_dir=$(cd -- "$3" && pwd -P)
else
  run_dir=$(mktemp -d "$script_dir/run-XXXXXXXX")
fi
echo "Output in $run_dir"
cd -- "$run_dir"

# The current reconstruction no longer recognizes these two vertex-fit keys.
reco_params="$run_dir/recoparams_optmized.ini"
awk '!/^[[:space:]]*(vt|ms)MaxChi2VertexConstrain[[:space:]]*=/' "$reco_params_source" > "$reco_params"

# Load the macro in each worker without shared ACLiC build files.
# Reconstruction clusters the simulated hits directly.
if ! "$reco_only"; then
  PYTHONUNBUFFERED=1 python3 "$repo_dir/scripts/na6psim_parallel.py" \
    --workers "$workers" --nevents "$events" \
    --generator "$generator_call" --load-ini "$layout" \
    --doDigitization false \
    --configKeyValues "keyval.output_dir=$run_dir/sim" \
    2>&1 | tee "$run_dir/sim.log"
fi

PYTHONUNBUFFERED=1 python3 "$repo_dir/scripts/na6prec_parallel.py" \
  --workers "$workers" --nevents "$events" --input-dir "$run_dir/sim" \
  --load-ini "$layout" --load-recoparam "$reco_params" \
  --geometry "$run_dir/sim/geometry.root" \
  --configKeyValues "keyval.output_dir=$run_dir/reco;reco.mtMCMatching=false;reco.mtMinVTHits=3" \
  2>&1 | tee "$run_dir/reco.log"

# Rerun matching on the same reconstructed tracks with the tighter VT requirement.
for threshold in 4 5; do
(
  mkdir -- "$run_dir/reco_${threshold}hits"
  cd -- "$run_dir/reco_${threshold}hits"
  for input in TracksVerTel.root TracksMuonSpec.root ClustersVerTel.root ClustersMuonSpec.root MCKine.root; do
    ln -s -- "$run_dir/reco/$input" "$input"
  done
  na6prec --load-ini "$layout" --load-recoparam "$reco_params" \
    --geometry "$run_dir/sim/geometry.root" \
    --doHitsToRecPoints false --doDigitsToRecPoints false \
    --doTrackletVertex false --doVTTracking false --doMSTracking false --doMatching true \
    --configKeyValues "keyval.output_dir=$PWD;reco.mtMCMatching=false;reco.mtMinVTHits=$threshold"
) 2>&1 | tee "$run_dir/reco_${threshold}hits.log"
done

root -l -b -q "$repo_dir/thermal_muon_short/matchingQuality.C(\"$run_dir\")" 2>&1 | tee "$run_dir/matchingQuality.log"
root -l -b -q "$repo_dir/thermal_virtual_photon/plotMatchingPNG.C(\"$run_dir\",$events)" 2>&1 | tee "$run_dir/plotMatchingPNG.log"
