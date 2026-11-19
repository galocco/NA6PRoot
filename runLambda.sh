#!/bin/bash

NJOB=1000   # numero di job paralleli

run_job() {
    ID=$1
    WORKDIR=FullEv_40AGeV/job_$ID
    mkdir -p $WORKDIR
    cd $WORKDIR
    
    na6psim -n1000 -g "$NA6PROOT_ROOT/share/test/genStrangenessBgEvent.C+"  --load-ini ../../na6pLayout_Jelle.ini > sim.log 2>&1
    na6prec --load-recoparam ../../na6pRecoParam.ini --load-ini ../../na6pLayout_Jelle.ini --doMSTracking false --doMatching false > rec.log 2>&1
    root -l -q '../../test/BuildCandidates.C("Lambda0")' > build.log 2>&1
    cd ../..
}

export -f run_job
export NA6PROOT_ROOT

parallel -j $NJOB run_job {} ::: $(seq 1 $NJOB)
DIR=FullEv_40AGeV

# merge log files
cat $DIR/job_*/sim.log  > $DIR/sim_all.log
cat $DIR/job_*/rec.log  > $DIR/rec_all.log
cat $DIR/job_*/build.log > $DIR/build_all.log

rm $DIR/*root
# merge the results
hadd -f Lambda0-Candidates.root $DIR/job_*/Lambda0-Candidates.root

cp $DIR/job_1/na6pLayout.ini $DIR/
cp $DIR/job_1/na6pRecoParam.ini $DIR/
cp $DIR/job_1/geometry.root $DIR/

cp Lambda0-Candidates.root $DIR/
cd $DIR
root -q -b -l '../test/PlotLambda0EffPurity.C("Lambda0-Candidates.root","Lambda0-EffPurity","na6pLayout.ini")' > plot.log 2>&1