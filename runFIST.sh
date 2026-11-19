na6psim -n1 \
  -g "$NA6PROOT_ROOT/share/test/genHEPMC.C+(\"/home/galocco/Thermal-FIST/fist-sampler/output.root\")" \
  --load-ini ../na6pLayout_Jelle.ini

na6prec \
    --load-recoparam ../na6pRecoParam.ini \
    --load-ini ../na6pLayout_Jelle.ini \
    --configKeyValues "$SETTING" \
    --doMSTracking false \
    --doMatching false \
    --doTOFMatching false \
    --doTOFHitsToRecPoints false
