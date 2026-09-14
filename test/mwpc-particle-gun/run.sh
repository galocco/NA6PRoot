#!/usr/bin/env bash
set -euo pipefail

EVENTS="${1:-100}"
MOMENTUM_GEV="${2:-20}"
PDG="${3:-13}"
SEED="${4:-20260914}"

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO="$(cd "$HERE/../.." && pwd)"
BUILD="$HERE/build"
RUNS="$HERE/runs"

if [[ -z "${NA6PROOT_ROOT:-}" ]]; then
  if [[ -f "$REPO/install/init.sh" ]]; then
    # Use the local installation when it is available.
    # shellcheck disable=SC1091
    source "$REPO/install/init.sh"
  else
    echo "NA6PROOT_ROOT is not set and $REPO/install/init.sh was not found." >&2
    echo "Load the NA6PRoot/O2 environment and build/install this branch first." >&2
    exit 2
  fi
fi

if [[ ! -f "$NA6PROOT_ROOT/include/NA6PMWPCChamber.h" ]]; then
  cat >&2 <<EOF
The installed NA6PRoot at
  $NA6PROOT_ROOT
does not contain NA6PMWPCChamber.h yet.

Build and install the current feature/mwpc-dice branch first, for example:
  cd "$REPO"
  mkdir -p build && cd build
  cmake -DCMAKE_INSTALL_PREFIX="$REPO/install" ..
  make -j5 install
  source "$REPO/install/init.sh"

Then rerun this script.
EOF
  exit 2
fi

mkdir -p "$BUILD" "$RUNS"
cmake -S "$HERE" -B "$BUILD" -DCMAKE_BUILD_TYPE=Release -DNA6P_INSTALL="$NA6PROOT_ROOT"
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
