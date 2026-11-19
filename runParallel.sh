#!/bin/bash

NJOB=100   # numero di job paralleli

run_job() {
    ID=$1
    WORKDIR=jobs_d0_158gev/job_$ID
    mkdir -p $WORKDIR
    cd $WORKDIR
    
    na6psim -n1000 -g "$NA6PROOT_ROOT/share/test/genDimuonBgEvent.C+(1,\"Jpsi\",0.,5.,0,6,false)"  --load-ini ../na6pLayout_Jelle.ini

    na6prec --load-recoparam ../../na6pRecoParam.ini \
        --load-ini ../../na6pLayout_Jelle.ini

    cd ..
}

export -f run_job
export NA6PROOT_ROOT

parallel -j $NJOB run_job {} ::: $(seq 1 $NJOB)