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

    local DIR="D0_${MAT}"
    mkdir -p "$DIR"
    cd "$DIR" || exit 1

    echo "Running with setting: $SETTING"

    #na6psim -n50000 \
    #    -g "$NA6PROOT_ROOT/share/test/genOpenCharmBgEvent.C+(1,\"D0\",0.,5.,0,6,false)" \
    #    --load-ini ../na6pLayout_Jelle.ini \
    #    --configKeyValues "$SETTING"

    #na6prec \
    #    --load-recoparam ../na6pRecoParam.ini \
    #    --load-ini ../na6pLayout_Jelle.ini \
    #    --configKeyValues "$SETTING" \
    #    --doMSTracking false \
    #    --doMatching false

    root -l -q ../plotVTD0.C
}

export -f run_one
export NA6PROOT_ROOT

parallel -j 4 run_one ::: "${SETTINGS[@]}"

#root -l -q overlayMassSigmaVsRapidity.C