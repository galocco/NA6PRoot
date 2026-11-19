#!/usr/bin/env bash
set -euo pipefail



# Generate J/psi events
#na6psim_parallel --workers 10 --nevents 50000 -g "$NA6PROOT_ROOT/share/test/genDimuonBgEvent.C+(1,\"Jpsi\",0.,5.,0,6,false)"  --load-ini na6pLayout_Optimized.ini #../na6pLayout_Jelle.ini

start=$SECONDS
na6prec_parallel --workers 10 \
        --doHitsToRecPoints true \
        --doTrackletVertex false \
        --doVTTracking true  \
        --doMSTracking false \
        --doMatching false \
        --load-ini na6pLayout_Optimized.ini \
        --load-recoparam na6pRecoParam.ini > log.txt 2>&1
echo "Time taken: $(( SECONDS - start )) seconds"


#root -l -b -q '../plotMSDimuons.C(443,4)'
root -l -b -q '../plotVTDimuons.C(443,5)'