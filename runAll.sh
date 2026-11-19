INI_FILES=(na6pLayout_Jelle_SiOnly.ini na6pLayout_Jelle_SiOnly.ini  na6pLayout_Jelle_Air.ini na6pLayout_Jelle.ini)
DIRS=(si_all si_abs air_abs std_abs)

INI_FILES=(na6pLayout_Jelle.ini)
DIRS=(std_abs_newbb)

parallel --link '
  mkdir -p {2} &&
  cd {2} &&
  na6psim -n10000 -g '$NA6PROOT_ROOT/share/test/genDimuonBgEvent.C+' --load-ini ../na6pLayout_Jelle.ini &&
  na6prec --load-recoparam ../na6pRecoParam.ini &&
  root -l -q ../plotMSDimuons.C &&
  root -l -q ../plotVTDimuons.C
' ::: "${INI_FILES[@]}" ::: "${DIRS[@]}"

genDimuonBgEvent(1, "Jpsi", 0., 5., 0, 5, true)
genDimuonBgEvent(1, "Omega", 0., 5., 0, 5, true)