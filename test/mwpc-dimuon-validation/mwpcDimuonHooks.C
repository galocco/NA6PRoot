#include <cmath>
#include <fairlogger/Logger.h>
#include <TVirtualMC.h>

#include "NA6PMCStack.h"
#include "NA6PSimMisc.h"

namespace
{
void forceDimuonDecay(int pdg, const char* name)
{
  int decayMode[6][3] = {{0}};
  float branching[6] = {100.f, 0.f, 0.f, 0.f, 0.f, 0.f};
  decayMode[0][0] = 13;  // mu-
  decayMode[0][1] = -13; // mu+
  TVirtualMC::GetMC()->SetDecayMode(pdg, branching, decayMode);
  LOGP(debug, "Forced {} (PDG {}) to mu+mu- for MWPC validation", name, pdg);
}

bool hasSignalAncestor(NA6PMCStack* stack, const TParticle* particle)
{
  int motherID = particle->GetFirstMother();
  while (motherID >= 0) {
    const auto* mother = stack->GetParticle(motherID);
    if (!mother) {
      break;
    }
    const int pdg = std::abs(mother->GetPdgCode());
    if (pdg == 223 || pdg == 333 || pdg == 443) {
      return true;
    }
    motherID = mother->GetFirstMother();
  }
  return false;
}
} // namespace

int mwpcDimuonHooks(int arg, bool inout)
{
  if (arg == UserHook::ADDParticles) {
    // NA6PMC applies some internal decay overrides between the entry and exit
    // calls of this hook.  Set the modes on both calls so the final VMC table
    // is guaranteed to contain the dimuon channels used in this validation.
    forceDimuonDecay(443, "J/psi");
    forceDimuonDecay(333, "phi");
    forceDimuonDecay(223, "omega");
    return 0;
  }

  if (arg == UserHook::SelectParticles) {
    if (!inout) {
      return 0;
    }

    auto* stack = static_cast<NA6PMCStack*>(TVirtualMC::GetMC()->GetStack());
    if (!stack) {
      LOGP(error, "MWPC dimuon validation hook could not access NA6PMCStack");
      return -1;
    }

    // Keep every generated signal daughter muon in MCKine.root, including
    // muons which miss all sensitive detector volumes.  This makes the parent
    // and daughter momentum spectra unbiased by detector acceptance.
    const int nTracks = stack->GetNtrack();
    const int nPrimaries = stack->GetNprimary();
    for (int i = nPrimaries; i < nTracks; ++i) {
      auto* particle = stack->GetParticle(i);
      if (!particle || std::abs(particle->GetPdgCode()) != 13) {
        continue;
      }
      if (hasSignalAncestor(stack, particle)) {
        particle->SetBit(UserHook::KeepParticleBit);
      }
    }
    return 0;
  }

  LOGP(error, "Unknown MWPC dimuon validation hook ID {}", arg);
  return -1;
}
