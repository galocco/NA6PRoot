#!/bin/bash
set -euo pipefail

NTASKS=10     # total job numbers
NWORKERS=10    # number of parallel jobs

#SYSTEM="PbPb.8.77.C0-5"
#SYSTEM="PbPb.17.3.C0-5"
SYSTEM="PbPb.7.5.C0-5"

FIST_BIN="/home/galocco/Thermal-FIST/fist-sampler/build/FISTSampler"
FIST_INPUT="/home/galocco/Thermal-FIST/fist-sampler/input/input.${SYSTEM}"
FIST_ASCII="${SYSTEM}.events.dat"
FIST_ROOT="output.${SYSTEM}.root"
OUTDIR="/data/galocco/TheFIST_${SYSTEM}-Lambda"
RUN_FIST=false
RUN_RECO=true

[ -x "$FIST_BIN" ] || { echo "ERROR: missing executable $FIST_BIN"; exit 1; }
[ -f "$FIST_INPUT" ] || { echo "ERROR: missing input $FIST_INPUT"; exit 1; }
[ -n "${NA6PROOT_ROOT:-}" ] || { echo "ERROR: NA6PROOT_ROOT is not set"; exit 1; }
export SYSTEM OUTDIR FIST_BIN FIST_INPUT FIST_ASCII FIST_ROOT RUN_FIST RUN_RECO NA6PROOT_ROOT
run_job() {
    ID=$1
    WORKDIR="$OUTDIR/job_$ID"
    mkdir -p "$WORKDIR"
    echo "Running job $ID in $WORKDIR"
    cd "$WORKDIR"
    if [ "$RUN_FIST" = true ]; then
        "$FIST_BIN" "$FIST_INPUT" > "output.${SYSTEM}.txt" 2>&1
        [ -s "$FIST_ASCII" ] || { echo "ERROR: $FIST_ASCII not produced in $WORKDIR"; return 1; }
        root -l -q -b "/home/galocco/Thermal-FIST/fist-sampler/hepmc3_to_root.C(\"${FIST_ASCII}\", \"${FIST_ROOT}\")" > convert.log 2>&1
        [ -s "$FIST_ROOT" ] || { echo "ERROR: $FIST_ROOT not produced in $WORKDIR"; return 1; }
        rm ${FIST_ASCII}

    fi

    if [ "$RUN_RECO" = true ]; then
         #na6psim -n1 -g "${NA6PROOT_ROOT}/share/test/genHEPMC.C+(\"${FIST_ROOT}\")" --load-ini /home/galocco/NA6PRoot/na6pLayout_Jelle.ini > sim.log 2>&1
        na6prec --load-ini /home/galocco/NA6PRoot/na6pLayout_Jelle.ini --doMSTracking false --doMatching false > rec.log 2>&1
        root -l -q /home/galocco/NA6PRoot/test/PlotChargeEfficiency.C > build.log 2>&1
        #root -l -q '/home/galocco/NA6PRoot/test/BuildCandidates.C("Lambda0")' > build.log 2>&1
        cd ../..
    fi
}

export -f run_job
export NA6PROOT_ROOT
export OUTDIR FIST_BIN FIST_INPUT FIST_ASCII FIST_ROOT RUN_FIST RUN_RECO

parallel --halt soon,fail=1 -j "$NWORKERS" run_job {} ::: $(seq 1 "$NTASKS")

if [ "$RUN_RECO" = true ]; then
    DIR="$OUTDIR"
    cd "$DIR"
    # merge log files
    cat "$DIR"/job_*/sim.log  > "$DIR"/sim_all.log
    cat "$DIR"/job_*/rec.log  > "$DIR"/rec_all.log
    cat "$DIR"/job_*/build.log > "$DIR"/build_all.log

    # merge the results
    hadd -f Lambda0-Candidates.root "$DIR"/job_*/Lambda0-Candidates.root

    cp "$DIR"/job_1/na6pLayout.ini "$DIR"/
    cp "$DIR"/job_1/na6pRecoParam.ini "$DIR"/
    cp "$DIR"/job_1/geometry.root "$DIR"/

    root -q -b -l '/home/galocco/NA6PRoot/test/PlotLambda0EffPurity.C("Lambda0-Candidates.root","Lambda0-EffPurity","na6pLayout.ini")' > plot.log 2>&1
fi