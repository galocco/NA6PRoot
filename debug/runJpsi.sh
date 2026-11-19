#!/usr/bin/env bash
set -euo pipefail



# Generate J/psi events
na6psim_parallel --workers 10 --nevents 100 -g "$NA6PROOT_ROOT/share/test/genDimuonBgEvent.C+(1,\"Jpsi\",0.,5.,0,6,true)"  --load-ini na6pLayout_Optimized.ini #../na6pLayout_Jelle.ini

start=$SECONDS
na6prec --lastevent 99 \
        --doHitsToRecPoints true \
        --doTrackletVertex false \
        --doVTTracking true  \
        --doMSTracking true \
        --doMatching false \
        --load-ini na6pLayout_Optimized.ini \
        --load-recoparam na6pRecoParam.ini \
        --configKeyValues "reco.caDoInwardSeedExtension=true" > log.txt 2>&1
echo "na6prec took $(( SECONDS - start )) seconds to run"

mkdir -p caDoInwardSeedExtension
cp *png caDoInwardSeedExtension/
cp *txt caDoInwardSeedExtension/
cp *ini caDoInwardSeedExtension/

start=$SECONDS
na6prec --lastevent 99 \
        --doHitsToRecPoints true \
        --doTrackletVertex false \
        --doVTTracking true  \
        --doMSTracking true \
        --doMatching false \
        --load-ini na6pLayout_Optimized.ini \
        --load-recoparam na6pRecoParam.ini \
        --configKeyValues "reco.caDoInwardSeedExtension=false" > log.txt 2>&1
echo "na6prec took $(( SECONDS - start )) seconds to run"
mkdir -p caDoInwardSeedExtensionFalse
cp *png caDoInwardSeedExtensionFalse/
cp *txt caDoInwardSeedExtensionFalse/
cp *ini caDoInwardSeedExtensionFalse/