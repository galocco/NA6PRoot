#if !defined(__CINT__) || defined(__MAKECINT__)
#include "NA6PGenParamTH3.h"
#endif

// Example generator for D0 from a TH2(pt,y) stored in data/
NA6PGenerator* genD0()
{
  // Update these two to match your ROOT input in data/
  const std::string sourceFile = "data/pp0_frag-PtSpectra-Boost.root";
  const std::string sourceHisto = "hPtY_D0";

  // name, PDG, multiplicity, source ROOT file, TH2 name, poisson multiplicity
  auto gen = new NA6PGenParamTH3("D0_th2", 421, 1.0f, sourceFile, sourceHisto, true);

  return gen;
}
