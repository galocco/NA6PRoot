// NA6PCCopyright

#include "Propagator.h"
#include "fairlogger/Logger.h"
#include <TGeoGlobalMagField.h>
#include <TGeoManager.h>
#include <TFile.h>
#include <algorithm>

namespace
{
// Groom, Mokhov and Striganov, PDG Atomic and Nuclear Properties of Materials (MUB).
// Total b(E) in 10^-6 cm^2/g, including bremsstrahlung, pair production and photonuclear loss.
constexpr float muonEnergyGrid[] = {2.f, 5.f, 10.f, 20.f, 50.f, 100.f, 1000.f, 10000.f, 100000.f};
struct MuonRadiativeRow {
  float z;
  float b[std::size(muonEnergyGrid)];
};
constexpr MuonRadiativeRow muonRadiativeTable[] = {
  {4, {0.6925, 0.8662, 0.9822, 1.1030, 1.2787, 1.4179, 1.8043, 2.0432, 2.1973}}, // Be
  {6, {0.8175, 1.0848, 1.2750, 1.4747, 1.7590, 1.9711, 2.5470, 2.8583, 3.0292}}, // C
  {8, {0.9120, 1.2531, 1.5007, 1.7603, 2.1273, 2.3955, 3.1067, 3.4669, 3.6488}}, // O
  {13, {1.1189, 1.6230, 1.9974, 2.3869, 2.9339, 3.3256, 4.3192, 4.7837, 4.9873}}, // Al
  {14, {1.1915, 1.7485, 2.1643, 2.5962, 3.2023, 3.6349, 4.7250, 5.2270, 5.4403}}, // Si
  {26, {1.6337, 2.5534, 3.2517, 3.9583, 4.9622, 5.6621, 7.3528, 8.0708, 8.3323}}, // Fe
  {74, {2.6518, 4.7694, 6.4130, 7.9914, 10.2825, 11.8361, 15.4003, 16.7697, 17.1793}}, // W
  {82, {2.7173, 5.0077, 6.7899, 8.4967, 10.9682, 12.6440, 16.4724, 17.9313, 18.3613}}, // Pb
};

float muonRadiativeRate(float z, size_t low, float fraction)
{
  // Interpolate in Z for elements without an entry and clamp to the table endpoints.
  z = std::clamp(z, muonRadiativeTable[0].z, muonRadiativeTable[std::size(muonRadiativeTable) - 1].z);
  size_t upper = 1;
  while (z > muonRadiativeTable[upper].z) {
    ++upper;
  }
  const auto& a = muonRadiativeTable[upper - 1];
  const auto& b = muonRadiativeTable[upper];
  const float zFraction = (z - a.z) / (b.z - a.z);
  return std::lerp(std::lerp(a.b[low], a.b[low + 1], fraction),
                   std::lerp(b.b[low], b.b[low + 1], fraction), zFraction);
}
}

Propagator::Propagator(bool uninitialized)
{
  if (uninitialized) {
    return;
  }
  ///< construct checking if needed components were initialized
  updateField();
}

//____________________________________________________________
bool Propagator::loadField()
{
  if (dynamic_cast<MagneticField*>(TGeoGlobalMagField::Instance()->GetField())) {
    LOGP(warn, "Global magnetic field was already loaded");
    return true;
  }
  auto prop = Instance(true);
  auto f = new MagneticField();
  f->loadField();
  f->setAsGlobalField();
  prop->updateField();
  return true;
}

//____________________________________________________________
bool Propagator::loadGeometry(const std::string& path, const std::string geomName)
{
  if (gGeoManager) {
    LOGP(warn, "Geometry was already loaded");
    return true;
  }
  auto f = TFile::Open(path.c_str());
  if (!f || f->IsZombie()) {
    LOGP(error, "Could not open geometry file {}", path);
    return false;
  }
  if (!f->Get(geomName.c_str())) {
    LOGP(error, "Could not load geometry {} from file {}", geomName, path);
  }
  f->Close();
  delete f;
  return gGeoManager != nullptr;
}

//____________________________________________________________
void Propagator::updateField()
{
  if (!mField) {
    mField = static_cast<MagneticField*>(TGeoGlobalMagField::Instance()->GetField());
    if (!mField) {
      LOG(fatal) << "Magnetic field is not initialized!";
    }
  }
  mNDOF = getBy({0.f, 0.f, 0.f}) != 0 ? 5 : 4;
}

bool Propagator::hasGeometryLoaded() const
{
  return gGeoManager != nullptr;
}

template <typename T>
MatBudget Propagator::getMeanMaterialBudgetFromGeom(const T* startF, const T* endF, float muonEnergy, float* radiativeDEDX) const
{
  // "mparam" - parameters used for the energy and multiple scattering
  //  corrections:
  MatBudget mbTot;
  size_t eBin = 0;
  float eFraction = 0.f, eScale = 0.f, radiativeSum = 0.f;
  if (radiativeDEDX) {
    *radiativeDEDX = 0.f;
  }
  const bool addRadiative = radiativeDEDX && std::isfinite(muonEnergy);
  if (addRadiative) {
    const float e = std::clamp(muonEnergy, muonEnergyGrid[0], muonEnergyGrid[std::size(muonEnergyGrid) - 1]);
    while (e > muonEnergyGrid[eBin + 1]) {
      ++eBin;
    }
    eFraction = std::log(e / muonEnergyGrid[eBin]) / std::log(muonEnergyGrid[eBin + 1] / muonEnergyGrid[eBin]);
    // Below the table, taper b to zero at rest using the kinetic energy.
    const float m = PID::getMass(PID::Muon);
    eScale = 1.e-6f * muonEnergy * std::clamp((muonEnergy - m) / (muonEnergyGrid[0] - m), 0.f, 1.f);
  }
  if (gGeoManager == nullptr) {
    LOGP(fatal, "No geometry was loaded, cannot extract  material corrections");
  }
  float tolerance = 1e-9f;
  float dx = endF[0] - startF[0], dy = endF[1] - startF[1], dz = endF[2] - startF[2];
  float length = std::sqrt(dx * dx + dy * dy + dz * dz);
  if (length < tolerance) {
    mbTot.meanL = length;
    return mbTot;
  }
  float invlen = 1.f / length;
  double dir[3] = {dx * invlen, dy * invlen, dz * invlen}, start[3] = {startF[0], startF[1], startF[2]};
  TGeoNode* currNode = gGeoManager->InitTrack(start, dir);
  if (!currNode) {
    mbTot.meanL = 0.f;
    return mbTot;
  }
  float sumSteps = 0.0f, minStep = 1e-4f; // 1 micron

  while (sumSteps < length) {
    float stepMax = length - sumSteps;

    // Find next boundary (or max step to end point)
    TGeoNode* nextNode = gGeoManager->FindNextBoundaryAndStep(stepMax, false);
    if (!nextNode) {
      break;
    }
    float snext = std::min<float>(gGeoManager->GetStep(), stepMax);
    // Handle numerical zero-step
    if (snext < tolerance) {
      nextNode = gGeoManager->Step(std::min(minStep, stepMax));
      snext = std::min<float>(gGeoManager->GetStep(), stepMax);
      if (snext < tolerance) {
        break; // protection in case still zero
      }
    }
    TGeoNode* nodeForMaterial = currNode;
    currNode = nextNode ? nextNode : currNode;

    if (!nodeForMaterial) {
      sumSteps += snext;
      continue;
    }

    TGeoMedium* med = nodeForMaterial->GetVolume()->GetMedium();
    if (!med) {
      sumSteps += snext;
      continue;
    }
    TGeoMaterial* mat = med->GetMaterial();
    if (!mat) {
      sumSteps += snext;
      continue;
    }
    const float radLen = mat->GetRadLen();
    if (radLen > 0.f) {
      mbTot.meanX2X0 += snext / radLen;
    }
    const float gram = snext * mat->GetDensity();
    mbTot.meanRho += gram;
    mbTot.meanA += snext * mat->GetA();
    mbTot.meanZ += snext * mat->GetZ();
    if (mat->IsMixture()) {
      TGeoMixture* mixture = (TGeoMixture*)mat;
      float meanZ2A = 0.f, sum = 0.f, b = 0.f;
      for (int iel = 0; iel < mixture->GetNelements(); iel++) {
        const float a = mixture->GetAmixt()[iel], w = mixture->GetWmixt()[iel];
        if (a > 0.f) {
          sum += w;
          meanZ2A += mixture->GetZmixt()[iel] * w / a;
          if (addRadiative) {
            b += w * muonRadiativeRate(mixture->GetZmixt()[iel], eBin, eFraction);
          }
        }
      }
      if (sum > 0.f) {
        mbTot.meanZ2A += snext * meanZ2A / sum;
        if (addRadiative) {
          radiativeSum += gram * b / sum;
        }
      }
    } else if (mat->GetA() > 0.f) {
      mbTot.meanZ2A += snext * mat->GetZ() / mat->GetA();
      if (addRadiative) {
        radiativeSum += gram * muonRadiativeRate(mat->GetZ(), eBin, eFraction);
      }
    }
    sumSteps += snext;
    mbTot.nBoundaries += 1;
    currNode = nextNode;
  }

  if (addRadiative && mbTot.meanRho > 0.f) {
    *radiativeDEDX = eScale * radiativeSum / mbTot.meanRho;
  }
  if (sumSteps == 0.f) {
    return mbTot;
  }
  float norm = 1.f / sumSteps;
  mbTot.meanRho *= norm;
  mbTot.meanA *= norm;
  mbTot.meanZ *= norm;
  mbTot.meanZ2A *= norm;
  mbTot.meanL = sumSteps;
  return mbTot;
}

bool Propagator::propagateToZ(NA6PTrackParCov& track, float zToGo, const Propagator::PropOpt& opt) const
{
  //----------------------------------------------------------------
  //
  // Propagates the track to the plane at z (cm)
  // and correcting for the crossed material.
  // opt.maxStep  - maximal step for propagation
  // opt.byOnly   - pure dipole propagation
  // opt.matCorr  - material correction type, it is up to the user to make sure the pointer is attached (if LUT is requested)
  //----------------------------------------------------------------
  auto dz = zToGo - track.getZ();
  int dir = dz > 0.f ? 1 : -1;
  NA6PTrackPar& linRef = opt.linRef ? *opt.linRef : track;
  while (std::abs(dz) > Epsilon) {
    auto step = std::min(std::abs(dz), opt.maxStep);
    if (dir < 0) {
      step = -step;
    }
    auto z = track.getZ() + step;
    auto xyz0 = linRef.getXYZ();
    auto correct = [&track, &linRef, &xyz0, dir, &opt, this]() {
      bool res = true;
      if (opt.matCorr != MatCorrType::USEMatCorrNONE) {
        auto xyz1 = linRef.getXYZ();
        float radiativeDEDX = 0.f;
        const auto mb = opt.muonRadiativeLoss && track.getPID() == PID::Muon
                          ? getMeanMaterialBudgetFromGeom(xyz0.data(), xyz1.data(), std::hypot(linRef.getP(), PID::getMass(PID::Muon)), &radiativeDEDX)
                          : getMeanMaterial(xyz0, xyz1);
        if (!track.correctForMaterial(mb.meanX2X0, (dir > 0 ? -mb.meanRho : mb.meanRho) * mb.meanL, mb.meanRho, mb.meanZ, mb.meanZ2A, linRef, false, radiativeDEDX)) {
          res = false;
        }
      }
      return res;
    };
    if (!(opt.byOnly ? track.propagateToZ(z, getBy(xyz0), opt.linRef) : track.propagateToZ(z, getFieldXYZ(xyz0), opt.linRef))) {
      return false;
    }
    if (!correct()) {
      return false;
    }
    dz = zToGo - track.getZ();
  }
  track.setZ(zToGo);
  if (opt.linRef) {
    opt.linRef->setZ(zToGo);
  }
  return true;
}

bool Propagator::propagateToZ(NA6PTrackPar& track, float zToGo, const Propagator::PropOpt& opt) const
{
  //----------------------------------------------------------------
  //
  // Propagates the track to the plane at z (cm)
  // and correcting for the crossed material.
  // opt.maxStep  - maximal step for propagation
  // opt.byOnly   - pure dipole propagation
  // opt.matCorr  - material correction type, it is up to the user to make sure the pointer is attached (if LUT is requested)
  //----------------------------------------------------------------
  auto dz = zToGo - track.getZ();
  int dir = dz > 0.f ? 1 : -1;
  while (std::abs(dz) > Epsilon) {
    auto step = std::min(std::abs(dz), opt.maxStep);
    if (dir < 0) {
      step = -step;
    }
    auto z = track.getZ() + step;
    auto xyz0 = track.getXYZ();
    auto correct = [&track, &xyz0, dir, &opt, this]() {
      bool res = true;
      if (opt.matCorr != MatCorrType::USEMatCorrNONE) {
        auto xyz1 = track.getXYZ();
        float radiativeDEDX = 0.f;
        const auto mb = opt.muonRadiativeLoss && track.getPID() == PID::Muon
                          ? getMeanMaterialBudgetFromGeom(xyz0.data(), xyz1.data(), std::hypot(track.getP(), PID::getMass(PID::Muon)), &radiativeDEDX)
                          : getMeanMaterial(xyz0, xyz1);
        if (!track.correctForELoss((dir > 0 ? -mb.meanRho : mb.meanRho) * mb.meanL, mb.meanRho, mb.meanZ, mb.meanZ2A, false, radiativeDEDX)) {
          res = false;
        }
      }
      return res;
    };
    if (!(opt.byOnly ? track.propagateParamToZ(z, getBy(xyz0)) : track.propagateParamToZ(z, getFieldXYZ(xyz0)))) {
      return false;
    }
    if (!correct()) {
      return false;
    }
    dz = zToGo - track.getZ();
  }
  track.setZ(zToGo);
  return true;
}

template <class T>
bool Propagator::propagatePCAToLine(T& track, float x0, float y0, float tolerance, const Propagator::PropOpt& opt) const
{
  // propagate the track to the point of closest aproach to the straight line x(z) = x0, y(z) = y0
  float zca = 0;
  int iterLeft = 2;
  while (iterLeft--) {
    if (!track.getZPCAToLine(x0, y0, getBy(track.getXYZ()), zca) || !propagateToZ(track, zca, opt)) {
      return false;
    }
    float dx = x0 - track.getX(), dy = y0 - track.getY();
    if (dx * dx + dy * dy < tolerance * tolerance) {
      break;
    }
  }
  return true;
}

template MatBudget Propagator::getMeanMaterialBudgetFromGeom<float>(const float* start, const float* end, float muonEnergy, float* radiativeDEDX) const;
template MatBudget Propagator::getMeanMaterialBudgetFromGeom<double>(const double* start, const double* end, float muonEnergy, float* radiativeDEDX) const;
template bool Propagator::propagatePCAToLine(NA6PTrackPar& track, float x0, float y0, float tolerance, const Propagator::PropOpt& opt) const;
template bool Propagator::propagatePCAToLine(NA6PTrackParCov& track, float x0, float y0, float tolerance, const Propagator::PropOpt& opt) const;
