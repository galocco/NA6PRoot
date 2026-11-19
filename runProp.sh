DIRS=(
  "propStepSize_0.5"
  "propStepSize_0.2"
  "propStepSize_0.1"
)

STEPSIZES=(
  "0.5"
  "0.2"
  "0.1"
)

parallel --link '
  mkdir -p {1} &&
  cp std_abs_newbb/*.root {1}/ &&
  cd {1} &&
  na6prec --load-recoparam ../na6pRecoParam.ini --load-ini ../na6pLayout_Jelle.ini --configKeyValues "reco.maxStepForProp={2};" --doMatching false --doVTTracking false --doTrackletVertex false --doHitsToRecPoints false &&
  root -l -q ../plotMSDimuons.C
' ::: "${DIRS[@]}" ::: "${STEPSIZES[@]}"

