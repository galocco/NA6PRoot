#!/usr/bin/env bash
SETTINGS=(
    "layout.useCoolingPlate=false;layout.coolingPlateThickness=0.4;layout.medPlateVerTel=Air"
    "layout.useCoolingPlate=true;layout.coolingPlateThickness=0.4;layout.medPlateVerTel=Air"
    "layout.useCoolingPlate=true;layout.coolingPlateThickness=0.4;layout.medPlateVerTel=Aluminium"
    "layout.useCoolingPlate=true;layout.coolingPlateThickness=0.4;layout.medPlateVerTel=CarbonFoam"
)


run_one() {
    local SETTING="$1"
    local MAT
    MAT=$(echo "$SETTING" | sed -n 's/.*layout.medPlateVerTel=\([^;]*\).*/\1/p')

    local DIR="LAMBDA_C_${MAT}"
    mkdir -p "$DIR"
    cd "$DIR" || exit 1

    echo "Running with setting: $SETTING"

    root -l -q ../plotVTLambdaC.C
}

export -f run_one
export NA6PROOT_ROOT

parallel -j 4 run_one ::: "${SETTINGS[@]}"

#root -l -q overlayMassSigmaVsRapidity.C