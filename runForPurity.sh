LAYOUT=$NA6PROOT_ROOT/../na6pLayout_Jelle.ini
FIELD_VALUES="0.2 0.4 0.6 0.8 1.0 1.2 1.4"

OPTIONS="reco.msNLayers=4;reco.msMinNClusTracksCA[0]=4;reco.msMinNClusTracksCA[1]=4;"

OPTIONAL=(
          "reco.msDoOutwardPropagation=false;reco.msDoInwardRefit=false;"
          "reco.msDoOutwardPropagation=true;reco.msDoInwardRefit=false;"
          "reco.msDoOutwardPropagation=true;reco.msDoInwardRefit=true;"
          )

RESX=(0.01 0.02 0.05 0.005 0.002 0.0005)
RESY=(0.05 0.10 0.25 0.025 0.010 0.0005)

for opt in "${OPTIONAL[@]}"; do
  echo "Running with options: $opt"

  run_field () {
    field=$1
    mkdir -p field_scaled${field}
    cd field_scaled${field} || exit

    echo "Running with field scale factor: $field"

    na6psim --load-ini $LAYOUT -n 10000 \
      --configKeyValues "layout.BScaleFactor=$field" \
      -g $NA6PROOT_ROOT/share/test/genDimuonBgEvent.C+ > output.log

    RESX=$(echo "0.01*$field" | bc)
    RESY=$(echo "0.05*$field" | bc)

    na6prec --load-ini na6pLayout.ini \
      --configKeyValues "reco.msSpatialResolutionX=${RESX};reco.msSpatialResolutionY=${RESY};${OPTIONS};${opt}" \
      > output2.log
    
    root -l -q $NA6PROOT_ROOT/../test/checkMSPurity.C
  }

  export -f run_field
  export LAYOUT NA6PROOT_ROOT

  parallel -j 7 run_field ::: $FIELD_VALUES

  length=${#RESX[@]}
  for ((i=0; i<length; i++)); do
    echo "resX: ${RESX[i]}, resY: ${RESY[i]}"
    rm -rf field_scaled1.0_resX${RESX[i]}_resY${RESY[i]}
    cp -r field_scaled1.0 field_scaled1.0_resX${RESX[i]}_resY${RESY[i]}
    cd field_scaled1.0_resX${RESX[i]}_resY${RESY[i]}
    na6prec --load-ini $LAYOUT --configKeyValues "reco.msSpatialResolutionX=${RESX[i]};reco.msSpatialResolutionY=${RESY[i]};${OPTIONS};${opt}" > output.log
    root -l -q $NA6PROOT_ROOT/../test/checkMSPurity.C
    cd ..
  done

  root -l -q overlayPurity.C 
  mkdir -p $opt
  mv field_scaled* $opt/
  mv resolution* $opt/
done