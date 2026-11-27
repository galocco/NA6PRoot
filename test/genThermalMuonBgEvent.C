// NA6PCCopyright

#include "genBgEvent.C"

#include <TDatabasePDG.h>
#include <cmath>
#include <memory>

// Phenomenological single-muon spectrum: pT * exp(-mT / inverseSlope)
// and a Gaussian in lab rapidity centered at the beam's yCM.
// The same pT (GeV/c) and lab-rapidity ranges apply to the hadronic background.
NA6PGenerator* genThermalMuonBgEvent(float ptMin = 0.f, float ptMax = 5.f,
                                     float yMin = 0.f, float yMax = 6.f, bool bg = true,
                                     int nMuons = 2, float inverseSlope = 0.25f, float sigmaY = 1.f)
{
  if (nMuons < 1 || !std::isfinite(ptMin) || !std::isfinite(ptMax) || ptMin < 0.f || ptMax <= ptMin ||
      !std::isfinite(yMin) || !std::isfinite(yMax) || yMax <= yMin ||
      !std::isfinite(inverseSlope) || inverseSlope <= 0.f || !std::isfinite(sigmaY) || sigmaY <= 0.f) {
    LOGP(error, "Invalid thermal muon generator parameters: nMuons={}, pT=[{}, {}], y=[{}, {}], inverseSlope={}, sigmaY={}",
         nMuons, ptMin, ptMax, yMin, yMax, inverseSlope, sigmaY);
    return nullptr;
  }

  auto cocktail = std::make_unique<NA6PGenCocktail>("thermalMuonBgCocktail");
  cocktail->setOwner(true);
  if (bg) {
    addBgEventGenerator(cocktail.get(), ptMin, ptMax, yMin, yMax);
  }

  const float muonMass = TDatabasePDG::Instance()->GetParticle(13)->Mass();
  const float yCM = NA6PBeamParam::Instance().getYCM();
  const int nMinus = nMuons / 2;
  const int multiplicities[] = {nMuons - nMinus, nMinus};
  const int pdgCodes[] = {-13, 13};
  const char* names[] = {"thermalMuPlus", "thermalMuMinus"};
  for (int charge = 0; charge < 2; ++charge) {
    if (multiplicities[charge] == 0) {
      continue;
    }
    auto muons = std::make_unique<NA6PGenParam>(names[charge], pdgCodes[charge], multiplicities[charge],
                                              "x*exp(-sqrt(x*x+[0]*[0])/[1])", "exp(-0.5*pow((x-[0])/[1],2))",
                                              ptMin, ptMax, yMin, yMax, true, true, false);
    muons->setParametersTrans({muonMass, inverseSlope});
    muons->setParametersLong({yCM, sigmaY});
    cocktail->addGenerator(muons.release());
  }
  return cocktail.release();
}
