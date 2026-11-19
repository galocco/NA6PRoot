#!/usr/bin/env bash


na6psim -n10000 \
        -g "$NA6PROOT_ROOT/share/test/genDimuonBgEvent.C+(1,\"Jpsi\",0.,5.,0,6,false)" \
        --load-ini na6pLayout_optimized.ini 

na6prec \
        --load-recoparam ../na6pRecoParam.ini \
        --load-ini na6pLayout_optimized.ini

root -l -q '../plotMTDimuons.C(443)'
