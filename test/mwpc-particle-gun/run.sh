#!/usr/bin/env bash
set -euo pipefail

EVENTS="${1:-100}"
MOMENTUM_GEV="${2:-20}"
PDG="${3:-13}"
SEED="${4:-20260914}"

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO="$(cd "$HERE/../.." && pwd)"
LOCAL_INSTALL="$REPO/install"
BUILD="$HERE/build"
RUNS="$HERE/runs"

# This test is specifically meant to exercise the MWPC code from the current
# checkout.  A previously loaded NA6PRoot installation may still be present in
# LD_LIBRARY_PATH, so source the local install and put its bin/lib directories
# first.  Otherwise the dynamic loader can silently pick an older libsimLib.so.
if [[ ! -f "$LOCAL_INSTALL/init.sh" ]]; then
  cat >&2 <<EOF
$LOCAL_INSTALL/init.sh was not found.

Build and install the current feature/mwpc-dice checkout first:
  cd "$REPO"
  cmake -S . -B build -DCMAKE_INSTALL_PREFIX="$LOCAL_INSTALL"
  cmake --build build --parallel 5 --target install
EOF
  exit 2
fi

# shellcheck disable=SC1091
source "$LOCAL_INSTALL/init.sh"
export NA6PROOT_ROOT="$LOCAL_INSTALL"
export PATH="$LOCAL_INSTALL/bin:${PATH:-}"
export LD_LIBRARY_PATH="$LOCAL_INSTALL/lib:${LD_LIBRARY_PATH:-}"
export ROOT_INCLUDE_PATH="$LOCAL_INSTALL/include:${ROOT_INCLUDE_PATH:-}"

if [[ ! -f "$LOCAL_INSTALL/include/NA6PMWPCChamber.h" ]]; then
  cat >&2 <<EOF
The local NA6PRoot installation does not contain NA6PMWPCChamber.h yet.

Rebuild and install the current branch:
  cd "$REPO"
  cmake -S . -B build -DCMAKE_INSTALL_PREFIX="$LOCAL_INSTALL"
  cmake --build build --parallel 5 --target install

Then rerun this script.
EOF
  exit 2
fi

mkdir -p "$BUILD" "$RUNS"
cmake -S "$HERE" -B "$BUILD" -DCMAKE_BUILD_TYPE=Release -DNA6P_INSTALL="$LOCAL_INSTALL"
cmake --build "$BUILD" --parallel 2

STAMP="$(date +%Y%m%d-%H%M%S)"
RUN="$RUNS/${STAMP}-${EVENTS}ev-${MOMENTUM_GEV}GeV-pdg${PDG}"

"$BUILD/mwpc_particle_gun" "$RUN" "$EVENTS" "$MOMENTUM_GEV" "$PDG" "$SEED"

printf '\nSaved run: %s\n' "$RUN"
printf '\n--- summary.txt ---\n'
cat "$RUN/summary.txt"
printf '\n--- geometry_summary.txt ---\n'
cat "$RUN/geometry_summary.txt"
printf '\nFirst event records:\n'
head -n 6 "$RUN/events.csv"
