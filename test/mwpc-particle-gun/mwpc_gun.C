#if !defined(__CINT__) || defined(__MAKECINT__)
#include "NA6PGenerator.h"
#include "NA6PMCStack.h"
#include <TDatabasePDG.h>
#include <TMCProcess.h>
#include <cmath>
#include <stdexcept>
#endif

// One forward muon per event. The source is 10 cm upstream of the chamber
// centre, so a zero-angle track passes through the middle of the sensitive gas.
class MWPCGun final : public NA6PGenerator
{
 public:
  MWPCGun(double momentumGeV, int pdg)
    : NA6PGenerator("MWPCGun"), mMomentum(momentumGeV), mPdg(pdg)
  {
    if (mMomentum <= 0.) {
      throw std::runtime_error("MWPCGun momentum must be positive");
    }
    if (mPdg != 13 && mPdg != -13) {
      throw std::runtime_error("MWPCGun supports mu- (13) and mu+ (-13) only");
    }
  }

  void generate() override
  {
    constexpr double vx = 0.;
    constexpr double vy = 0.;
    constexpr double vz = -10.; // cm
    const double mass = TDatabasePDG::Instance()->GetParticle(mPdg)->Mass();
    const double energy = std::sqrt(mMomentum * mMomentum + mass * mass);

    auto* stack = getStack();
    auto* header = stack->getEventHeader();
    header->setVX(vx);
    header->setVY(vy);
    header->setVZ(vz);
    header->setEventID(mEvent);
    stack->setPVGenerated(true);

    int track = -1;
    stack->PushTrack(true, -1, mPdg,
                     0., 0., mMomentum, energy,
                     vx, vy, vz, 0.,
                     0., 0., 0., TMCProcess::kPPrimary,
                     track, 1., 0);

    header->getGenHeaders().emplace_back(1, 0, header->getNPrimaries(), 0, "MWPCGun");
    header->incNPrimaries(1);
    ++mEvent;
  }

 private:
  double mMomentum = 20.;
  int mPdg = 13;
  int mEvent = 0;
};

NA6PGenerator* mwpc_gun(double momentumGeV = 20., int pdg = 13)
{
  return new MWPCGun(momentumGeV, pdg);
}
