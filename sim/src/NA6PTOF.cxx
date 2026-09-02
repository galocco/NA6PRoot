// NA6PCCopyright

#include "NA6PTOF.h"

#include "NA6PLayoutParam.h"
#include "NA6PMCStack.h"
#include "NA6PTGeoHelper.h"

#include <TColor.h>
#include <TFile.h>
#include <TGeoBBox.h>
#include <TGeoManager.h>
#include <TGeoMatrix.h>
#include <TGeoVolume.h>
#include <TTree.h>
#include <TVirtualMC.h>
#include <fairlogger/Logger.h>

#include <array>

void NA6PTOF::createMaterials()
{
  auto& matPool = NA6PTGeoHelper::instance().getMatPool();
  const auto name = addName("Silicon");
  if (matPool.find(name) == matPool.end()) {
    matPool[name] = new TGeoMaterial(name.c_str(), 28.09, 14, 2.33);
    NA6PTGeoHelper::instance().addMedium(name, "", kCyan + 1);
  }
}

void NA6PTOF::createGeometry(TGeoVolume* world)
{
  const auto& param = NA6PLayoutParam::Instance();
  createMaterials();

  const float tileHalfSide = 0.5f * param.tofTileSide;
  const float tileHalfThickness = 0.5f * param.tofTileThickness;
  const float holeHalfSize = 0.5f * param.tofHoleSize;

  auto* tileShape = new TGeoBBox("TOFTileShape", tileHalfSide, tileHalfSide, tileHalfThickness);
  auto* tile = new TGeoVolume("TOFTile", tileShape, NA6PTGeoHelper::instance().getMedium(addName("Silicon")));
  tile->SetLineColor(NA6PTGeoHelper::instance().getMediumColor(addName("Silicon")));

  auto* layer = new TGeoVolumeAssembly("TOFLayer");
  constexpr std::array<float, NTiles> signsX = {1.f, -1.f, -1.f, 1.f};
  constexpr std::array<float, NTiles> signsY = {1.f, 1.f, -1.f, -1.f};
  for (int tileID = 0; tileID < NTiles; ++tileID) {
    // Arrange the tiles as a pinwheel so that only the central square is open.
    // A purely radial offset would instead leave a cross-shaped gap through the layer.
    const float tileX = signsX[tileID] * tileHalfSide + signsY[tileID] * holeHalfSize;
    const float tileY = signsY[tileID] * tileHalfSide - signsX[tileID] * holeHalfSize;
    layer->AddNode(tile, composeSensorVolID(tileID),
                   new TGeoTranslation(tileX, tileY, 0.f));
  }

  world->AddNode(layer, composeNonSensorVolID(0),
                 new TGeoTranslation(param.posTOF[0], param.posTOF[1], param.posTOF[2]));
}

void NA6PTOF::setAlignableEntries()
{
  const std::string topNodeName = gGeoManager->GetTopNode()->GetName();
  for (int tileID = 0; tileID < NTiles; ++tileID) {
    const int id = getActiveID() * 100 + tileID;
    const std::string name = fmt::format("TOF_Tile{}", tileID);
    const std::string path = fmt::format("/{}/TOFLayer_{}/TOFTile_{}", topNodeName,
                                         composeNonSensorVolID(0), composeSensorVolID(tileID));
    auto* entry = gGeoManager->SetAlignableEntry(name.c_str(), path.c_str(), id);
    if (entry) {
      LOGP(info, "Successfully added {} {} as alignable sensor {}", name, path, id);
    } else {
      LOGP(error, "FAILED to add alignable entry {} {}", name, path);
    }
  }
}

bool NA6PTOF::stepManager(int volID)
{
  const int sensorID = volID2SensID(volID);
  if (sensorID < 0) {
    LOGP(fatal, "Non-sensor volID={} was provided to stepManager of {}", volID, getName());
  }

  auto* mc = TVirtualMC::GetMC();
  bool startHit = false;
  bool stopHit = false;
  unsigned char status = 0;
  if (mc->IsTrackEntering()) {
    status |= NA6PBaseHit::kTrackEntering;
  }
  if (mc->IsTrackInside()) {
    status |= NA6PBaseHit::kTrackInside;
  }
  if (mc->IsTrackExiting()) {
    status |= NA6PBaseHit::kTrackExiting;
  }
  if (mc->IsTrackOut()) {
    status |= NA6PBaseHit::kTrackOut;
  }
  if (mc->IsTrackStop()) {
    status |= NA6PBaseHit::kTrackStopped;
  }
  if (mc->IsTrackAlive()) {
    status |= NA6PBaseHit::kTrackAlive;
  }

  if ((status & NA6PBaseHit::kTrackEntering) ||
      (status & NA6PBaseHit::kTrackInside && !mTrackData.mHitStarted)) {
    startHit = true;
  } else if (status & (NA6PBaseHit::kTrackExiting | NA6PBaseHit::kTrackOut | NA6PBaseHit::kTrackStopped)) {
    stopHit = true;
  }
  if (!startHit) {
    mTrackData.mEnergyLoss += mc->Edep();
  }
  if (!(startHit || stopHit)) {
    return false;
  }

  auto* stack = static_cast<NA6PMCStack*>(mc->GetStack());
  if (startHit) {
    mTrackData.mEnergyLoss = 0.;
    mc->TrackMomentum(mTrackData.mMomentumStart);
    mc->TrackPosition(mTrackData.mPositionStart);
    mTrackData.mTrkStatusStart = status;
    mTrackData.mHitStarted = true;
  }
  if (stopHit) {
    TLorentzVector positionStop, momentumStop;
    mc->TrackMomentum(momentumStop);
    mc->TrackPosition(positionStop);
    auto* hit = addHit(stack->GetCurrentTrackNumber(), sensorID, mTrackData.mPositionStart.Vect(),
                       positionStop.Vect(), mTrackData.mMomentumStart.Vect(), momentumStop.Vect(),
                       positionStop.T(), mTrackData.mEnergyLoss, mTrackData.mTrkStatusStart, status);
    mTrackData.mHitStarted = false;
    if (mVerbosity > 0) {
      LOGP(info, "{} Tr{} {}", getName(), stack->GetCurrentTrackNumber(), hit->asString());
    }
    stack->addHit(getActiveIDBit());
    return true;
  }
  return false;
}

NA6PTOFHit* NA6PTOF::addHit(int trackID, int detID, const TVector3& startPos, const TVector3& endPos,
                            const TVector3& startMom, const TVector3& endMom, float endTime,
                            float eLoss, unsigned char startStatus, unsigned char endStatus)
{
  mHits.emplace_back(trackID, detID, startPos, endPos, startMom, endMom, endTime, eLoss,
                     startStatus, endStatus);
  return &mHits.back();
}

void NA6PTOF::createHitsOutput(const std::string& outDir)
{
  const auto name = fmt::format("{}Hits{}.root", outDir, getName());
  mHitsFile = TFile::Open(name.c_str(), "recreate");
  mHitsTree = new TTree(fmt::format("hits{}", getName()).c_str(), fmt::format("{} Hits", getName()).c_str());
  mHitsTree->Branch(getName().c_str(), &hHitsPtr);
  LOGP(info, "Will store {} hits in {}", getName(), name);
}

void NA6PTOF::closeHitsOutput()
{
  if (mHitsTree && mHitsFile) {
    mHitsFile->cd();
    mHitsTree->Write();
    delete mHitsTree;
    mHitsTree = nullptr;
    mHitsFile->Close();
    delete mHitsFile;
    mHitsFile = nullptr;
  }
}

void NA6PTOF::writeHits(const std::vector<int>& remapping)
{
  for (size_t i = 0; i < mHits.size(); ++i) {
    auto& hit = mHits[i];
    if (remapping[hit.getTrackID()] < 0) {
      LOGP(error, "Track {} hit {} in {} was not remapped!", hit.getTrackID(), i, getName());
    }
    hit.setTrackID(remapping[hit.getTrackID()]);
  }
  if (mHitsTree) {
    mHitsTree->Fill();
  }
}
