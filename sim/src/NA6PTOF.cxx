// NA6PCCopyright

#include "NA6PTOF.h"
#include "NA6PDetector.h"
#include "NA6PTGeoHelper.h"
#include "NA6PLayoutParam.h"
#include "NA6PMCStack.h"

#include <TVirtualMC.h>
#include <TGeoVolume.h>
#include <TGeoNode.h>
#include <TGeoBBox.h>
#include <TGeoManager.h>
#include <TGeoCompositeShape.h>
#include <TGeoBoolNode.h>
#include <TColor.h>
#include <fairlogger/Logger.h>
#include <TFile.h>
#include <TTree.h>
#include <TMath.h>

// ============================================================
//  NA6PTOF::createMaterials
// ============================================================

void NA6PTOF::createMaterials()
{
  auto& matPool = NA6PTGeoHelper::instance().getMatPool();
  std::string nameM;

  nameM = addName("Silicon");
  if (matPool.find(nameM) == matPool.end()) {
    matPool[nameM] = new TGeoMaterial(nameM.c_str(), 28.09, 14, 2.33);
    NA6PTGeoHelper::instance().addMedium(nameM, "", kCyan + 1);
  }
  nameM = addName("CarbonFiber");
  if (matPool.find(nameM) == matPool.end()) {
    auto mixt = new TGeoMixture(nameM.c_str(), 1, 1.91);
    mixt->AddElement(12.01, 6, 1.0);
    matPool[nameM] = mixt;
    NA6PTGeoHelper::instance().addMedium(nameM, "", kGray + 1);
  }
  nameM = addName("Air");
  if (matPool.find(nameM) == matPool.end()) {
    auto mixt = new TGeoMixture(nameM.c_str(), 2, 0.001);
    mixt->AddElement(new TGeoElement("N", "Nitrogen", 7, 14.01), 0.78);
    mixt->AddElement(new TGeoElement("O", "Oxygen", 8, 16.00), 0.22);
    matPool[nameM] = mixt;
    NA6PTGeoHelper::instance().addMedium(nameM);
  }
}

// ============================================================
//  NA6PTOF::createGeometry
// ============================================================

void NA6PTOF::createGeometry(TGeoVolume* world)
{
  const auto& param = NA6PLayoutParam::Instance();

  createMaterials();

  // ── Pixel station dimensions ──────────────────────────────────────
  float pixChipContainerDX = 30.0f;
  float pixChipContainerDY = 30.0f;
  float pixChipContainerDz = 0.1f;
  float pixChipDz = 50e-4f;
  float pixChipOffsX = 0.29f;
  float pixChipOffsY = 0.31f;
  float carbonPlateDz = 400e-4f;

  float pixChipDX = 15;
  float pixChipDY = 15;
  float carbonPlateDX = 2 * pixChipDX + 2 * pixChipOffsX;
  float carbonPlateDY = 2 * pixChipDY + 2 * pixChipOffsY;

  TGeoMedium* medAir = NA6PTGeoHelper::instance().getMedium(addName("Air"));

  // ── TOF container ────────────────────────────────────────────────
  float boxDZMargin = pixChipContainerDz + 0.5f;
  float boxDZ = 2 * boxDZMargin; // single plane

  auto* tofShape = new TGeoBBox("TOFContainer",
                                pixChipContainerDX / 2.f + 1.f,
                                pixChipContainerDY / 2.f + 1.f,
                                boxDZ / 2);
  TGeoVolume* tofContainer = new TGeoVolume("TOFContainer", tofShape, medAir);

  auto* tofTransform = new TGeoCombiTrans(
    param.shiftTOF[0],
    param.shiftTOF[1],
    param.shiftTOF[2] + param.posTOFPlaneZ,
    NA6PTGeoHelper::rotAroundVector(0.0, 0.0, 0.0, 0.0));
  world->AddNode(tofContainer, composeNonSensorVolID(0), tofTransform);

  // ── Pixel sensor station ──────────────────────────────────────────
  auto* pixelStationShape = new TGeoBBox("TOFPixelStationShape",
                                         pixChipContainerDX / 2,
                                         pixChipContainerDY / 2,
                                         pixChipContainerDz / 2);
  auto* sensorShape = new TGeoBBox("TOFSensorShape",
                                   pixChipDX / 2, pixChipDY / 2, pixChipDz / 2);

  auto* pixelStationVol = new TGeoVolume("TOFPixelStationVol", pixelStationShape, medAir);
  TGeoVolume* pixelSensor = new TGeoVolume("TOFPixelSensor", sensorShape,
                                           NA6PTGeoHelper::instance().getMedium(addName("Silicon")));
  pixelSensor->SetLineColor(NA6PTGeoHelper::instance().getMediumColor(addName("Silicon")));

  // ── Carbon-fiber plate with central beam hole ─────────────────────
  auto* carbonplateFullShape = new TGeoBBox("TOFCarbonPlateFullShape",
                                            carbonPlateDX / 2, carbonPlateDY / 2, carbonPlateDz / 2);
  auto* beamHole = new TGeoBBox("TOFCarbonPlateBeamHole", pixChipOffsX, pixChipOffsY, carbonPlateDz);
  auto* holeRemoval = new TGeoSubtraction(carbonplateFullShape, beamHole);
  auto* cbPlateWithHoleShape = new TGeoCompositeShape("TOFCarbonPlateWithHoleShape", holeRemoval);
  TGeoVolume* cbPlate = new TGeoVolume("TOFCarbonPlateWithHole", cbPlateWithHoleShape,
                                       NA6PTGeoHelper::instance().getMedium(addName("CarbonFiber")));
  cbPlate->SetLineColor(NA6PTGeoHelper::instance().getMediumColor(addName("CarbonFiber")));

  // ── Place sensors + carbon plate into the station volume ──────────
  std::vector<float> alpdx = {pixChipDX / 2 + pixChipOffsX, -pixChipDX / 2 + pixChipOffsX,
                              -pixChipDX / 2 - pixChipOffsX, pixChipDX / 2 - pixChipOffsX};
  std::vector<float> alpdy = {pixChipDY / 2 - pixChipOffsY, pixChipDY / 2 + pixChipOffsY,
                              -pixChipDY / 2 + pixChipOffsY, -pixChipDY / 2 - pixChipOffsY};
  for (size_t ii = 0; ii < alpdx.size(); ++ii) {
    pixelStationVol->AddNode(pixelSensor, composeSensorVolID(ii),
                             new TGeoTranslation(alpdx[ii], alpdy[ii], 0));
  }
  pixelStationVol->AddNode(cbPlate, composeNonSensorVolID(20),
                           new TGeoCombiTrans(0., 0.,
                                              pixChipDz / 2 + carbonPlateDz / 2,
                                              NA6PTGeoHelper::rotAroundVector(0, 0.0, 0.0, 0.0)));

  // ── Place the single station in the container ─────────────────────
  tofContainer->AddNode(pixelStationVol, composeNonSensorVolID(0),
                        new TGeoCombiTrans(0., 0., 0.,
                                           NA6PTGeoHelper::rotAroundVector(0.0, 0.0, 0.0, 0.0)));
}

// ============================================================
//  NA6PTOF::setAlignableEntries
// ============================================================

void NA6PTOF::setAlignableEntries()
{
  int svolCnt = 0;
  for (int ii = 0; ii < 4; ++ii) {
    int id = getActiveID() * 100 + svolCnt;
    std::string nm = fmt::format("TOF_Lr0_Sens{}", ii);
    std::string path = fmt::format("/World_1/TOFContainer_{}/TOFPixelStationVol_{}/TOFPixelSensor_{}",
                                   composeNonSensorVolID(0),
                                   composeNonSensorVolID(0),
                                   composeSensorVolID(ii));
    gGeoManager->SetAlignableEntry(nm.c_str(), path.c_str(), id);
    LOGP(info, "Adding {} {} as alignable sensor {}", nm, path, id);
    svolCnt++;
  }
}

// ============================================================
//  NA6PTOF::stepManager
// ============================================================

bool NA6PTOF::stepManager(int volID)
{
  int sensID = NA6PModule::volID2SensID(volID);
  if (sensID < 0) {
    LOGP(fatal, "Non-sensor volID={} was provided to stepManager of {}", volID, getName());
  }
  auto mc = TVirtualMC::GetMC();
  bool startHit = false, stopHit = false;
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

  if ((status & NA6PBaseHit::kTrackEntering) || (status & NA6PBaseHit::kTrackInside && !mTrackData.mHitStarted)) {
    startHit = true;
  } else if ((status & (NA6PBaseHit::kTrackExiting | NA6PBaseHit::kTrackOut | NA6PBaseHit::kTrackStopped))) {
    stopHit = true;
  }
  if (!startHit) {
    mTrackData.mEnergyLoss += mc->Edep();
  }
  if (!(startHit | stopHit)) {
    return false;
  }
  auto stack = (NA6PMCStack*)mc->GetStack();
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
    int stationID(-1);
    mc->CurrentVolOffID(1, stationID);
    stationID = volID2NonSensID(stationID);

    int chipindex = NChipsPerStation * stationID + sensID;
    auto* p = addHit(stack->GetCurrentTrackNumber(), chipindex, mTrackData.mPositionStart.Vect(), positionStop.Vect(),
                     mTrackData.mMomentumStart.Vect(), momentumStop.Vect(), positionStop.T(),
                     mTrackData.mEnergyLoss, mTrackData.mTrkStatusStart, status);
    if (mVerbosity > 0) {
      LOGP(info, "{} Tr{} {}", getName(), stack->GetCurrentTrackNumber(), p->asString());
    }
    stack->addHit(getActiveIDBit());
    return true;
  }
  return false;
}

// ============================================================
//  NA6PTOF::addHit
// ============================================================

NA6PTOFHit* NA6PTOF::addHit(int trackID, int detID, const TVector3& startPos, const TVector3& endPos, const TVector3& startMom, const TVector3& endMom,
                             float endTime, float eLoss, unsigned char startStatus, unsigned char endStatus)
{
  mHits.emplace_back(trackID, detID, startPos, endPos, startMom, endMom, endTime, eLoss, startStatus, endStatus);
  return &(mHits.back());
}

// ============================================================
//  NA6PTOF::createHitsOutput
// ============================================================

void NA6PTOF::createHitsOutput(const std::string& outDir)
{
  auto nm = fmt::format("{}Hits{}.root", outDir, getName());
  mHitsFile = TFile::Open(nm.c_str(), "recreate");
  mHitsTree = new TTree(fmt::format("hits{}", getName()).c_str(), fmt::format("{} Hits", getName()).c_str());
  mHitsTree->Branch(getName().c_str(), &hHitsPtr);
  LOGP(info, "Will store {} hits in {}", getName(), nm);
}

// ============================================================
//  NA6PTOF::closeHitsOutput
// ============================================================

void NA6PTOF::closeHitsOutput()
{
  if (mHitsTree && mHitsFile) {
    mHitsFile->cd();
    mHitsTree->Write();
    delete mHitsTree;
    mHitsTree = 0;
    mHitsFile->Close();
    delete mHitsFile;
    mHitsFile = 0;
  }
}

// ============================================================
//  NA6PTOF::writeHits
// ============================================================

void NA6PTOF::writeHits(const std::vector<int>& remapping)
{
  int nh = mHits.size();
  for (int i = 0; i < nh; i++) {
    auto& h = mHits[i];
    if (remapping[h.getTrackID()] < 0) {
      LOGP(error, "Track {} hit {} in {} was not remapped!", h.getTrackID(), i, getName());
    }
    h.setTrackID(remapping[h.getTrackID()]);
  }
  if (mHitsTree) {
    mHitsTree->Fill();
  }
}
