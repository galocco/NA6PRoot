DIRS=(
  "matStepSize_1"
  "matStepSize_10"
  "matStepSize_100"
  "matStepSize_1000"
)

STEPSIZES=(
  "1"
  "10"
  "100"
  "1000"
)

parallel --link '
  mkdir -p {1} &&
  cp std_abs_newbb/*.root {1}/ &&
  cd {1} &&
  na6prec --load-recoparam ../na6pRecoParam.ini --load-ini ../na6pLayout_Jelle.ini --configKeyValues "reco.stepSizeForMat={2};" --doMatching false --doVTTracking false --doTrackletVertex false --doHitsToRecPoints false &&
  root -l -q ../plotMSDimuons.C
' ::: "${DIRS[@]}" ::: "${STEPSIZES[@]}"

