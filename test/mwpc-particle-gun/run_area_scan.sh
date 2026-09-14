#!/usr/bin/env bash
set -euo pipefail

# 700 x 700 mm = 70 x 70 cm. With 71 points per axis this gives a 1 cm grid
# from -35 to +35 cm in both x and y, at normal incidence.
GRID_N="${1:-71}"
AREA_CM="${2:-70}"
MOMENTUM_GEV="${3:-20}"
PDG="${4:-13}"
SEED="${5:-20260914}"

if (( GRID_N < 2 )); then
  echo "GRID_N must be at least 2" >&2
  exit 2
fi

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
EVENTS=$((GRID_N * GRID_N))

export MWPC_GUN_GRID_N="$GRID_N"
export MWPC_GUN_AREA_CM="$AREA_CM"

printf 'MWPC active-area scan\n'
printf '  grid: %s x %s = %s muons\n' "$GRID_N" "$GRID_N" "$EVENTS"
printf '  area: %s x %s cm\n' "$AREA_CM" "$AREA_CM"
printf '  momentum: %s GeV/c\n' "$MOMENTUM_GEV"
printf '  direction: normal to chamber (+z)\n\n'

bash "$HERE/run.sh" "$EVENTS" "$MOMENTUM_GEV" "$PDG" "$SEED"

RUN="$(find "$HERE/runs" -mindepth 1 -maxdepth 1 -type d -printf '%T@ %p\n' | sort -nr | head -n 1 | cut -d' ' -f2-)"

printf '\n--- active-area comparison ---\n'
python3 "$HERE/analyze_area_scan.py" "$RUN"

printf '\nDetailed point-by-point table:\n  %s/area_scan.csv\n' "$RUN"
printf 'Summary:\n  %s/area_scan_summary.txt\n' "$RUN"
