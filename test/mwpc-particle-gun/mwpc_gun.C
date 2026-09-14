#if !defined(__CINT__) || defined(__MAKECINT__)
#include "NA6PGenerator.h"
#include "NA6PMCStack.h"
#include <TDatabasePDG.h>
#include <TMCProcess.h>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <stdexcept>
#include <string>
#endif

// One forward muon per event. By default the gun is centred on the chamber.
// If MWPC_GUN_GRID_N and MWPC_GUN_AREA_CM are set, the source walks over a
// regular square grid while keeping every muon normal to the chamber plane.
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

    const char* gridN = std::getenv("MWPC_GUN_GRID_N");
    const char* area = std::getenv("MWPC_GUN_AREA_CM");
    if ((gridN && !area) || (!gridN && area)) {
      throw std::runtime_error("Set both MWPC_GUN_GRID_N and MWPC_GUN_AREA_CM, or neither");
    }
    if (gridN && area) {
      mGridN = std::stoi(gridN);
      mAreaCM = std::stod(area);
      if (mGridN < 2 || mAreaCM <= 0.) {
        throw std::runtime_error("MWPC grid scan requires GRID_N >= 2 and AREA_CM > 0");
      }
      mPositions.open("gun_xy.csv");
      if (!mPositions) {
        throw std::runtime_error("Cannot create gun_xy.csv");
      }
      mPositions << "event,x_cm,y_cm\n" << std::setprecision(17);
    }
  }

  void generate() override
  {
    double vx = 0.;
    double vy = 0.;
    constexpr double vz = -10.; // cm, safely upstream of the chamber

    if (mGridN > 1) {
      const int total = mGridN * mGridN;
      if (mEvent >= total) {
        throw std::runtime_error("MWPCGun received more events than grid points");
      }
      const int ix = mEvent % mGridN;
      const int iy = mEvent / mGridN;
      vx = -0.5 * mAreaCM + mAreaCM * static_cast<double>(ix) / (mGridN - 1);
      vy = -0.5 * mAreaCM + mAreaCM * static_cast<double>(iy) / (mGridN - 1);
      mPositions << mEvent << ',' << vx << ',' << vy << '\n';
      mPositions.flush();
    }

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
  int mGridN = 1;
  double mAreaCM = 0.;
  std::ofstream mPositions;
};

NA6PGenerator* mwpc_gun(double momentumGeV = 20., int pdg = 13)
{
  return new MWPCGun(momentumGeV, pdg);
}
