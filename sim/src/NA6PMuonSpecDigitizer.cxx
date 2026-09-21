// NA6PCCopyright

#include <ranges>
#include <numeric>
#include <fairlogger/Logger.h>
#include <TGeoNode.h>
#include <TGeoBBox.h>
#include <TGeoManager.h>
#include <TFile.h>
#include <TTree.h>

#include "NA6PLayoutParam.h"
#include "ConfigurableParam.h"
#include "StringUtils.h"
#include "NA6PMuonSpecDigitizer.h"

void NA6PMuonSpecDigitizer::init(const char* filename, const char* geoname)
{
  const auto& param = NA6PLayoutParam::Instance();
  const auto& mwpc = NA6PMWPCParam::Instance();
  mNumberOfModules = mwpc.getNModules(param.nMSPlanes);
  mModules.resize(mNumberOfModules);
  mThresholds.assign(mNumberOfModules, kDefaultThresholdEl);
  mSegmentation.configure(mwpc);
  initGeometry(filename, geoname);
  const auto outDir = na6p::utils::Str::rectifyDirectory(na6p::conf::ConfigurableParam::getOutputDir());
  createDigitsOutput(outDir);
}

void NA6PMuonSpecDigitizer::initGeometry(const char* filename, const char* geoname)
{
  if (!mGeoManager.loadGeometry(filename, geoname)) {
    LOGP(fatal, "Load of geometry not successful");
  }
}

void NA6PMuonSpecDigitizer::createDigitsOutput(const std::string& outDir)
{
  auto nm = fmt::format("{}Digits{}.root", outDir, getName());
  mDigitsFile = TFile::Open(nm.c_str(), "recreate");
  mDigitsTree = new TTree(fmt::format("digits{}", getName()).c_str(), fmt::format("{} Digits", getName()).c_str());
  mDigitsTree->Branch(getName().c_str(), &hDigitsPtr);
  mDigitsTree->Branch(fmt::format("{}MCTruth", getName()).c_str(), &hMCLabelsPtr);
  LOGP(info, "Will store {} digits in {}", getName(), nm);
}

void NA6PMuonSpecDigitizer::closeDigitsOutput()
{
  if (mDigitsTree && mDigitsFile) {
    mDigitsFile->cd();
    mDigitsTree->Write();
    delete mDigitsTree;
    mDigitsTree = nullptr;
    mDigitsFile->Close();
    delete mDigitsFile;
    mDigitsFile = nullptr;
  }
}

void NA6PMuonSpecDigitizer::writeDigits()
{
  if (mDigitsTree) {
    mDigitsTree->Fill();
  }
  LOGP(info, "Saved {} digits in tree with {} entries", mDigits.size(), mDigitsTree->GetEntries());
}

void NA6PMuonSpecDigitizer::process(const std::vector<NA6PMuonSpecHit>& hits, int layer)
{
  clearDigits();
  int nHits = hits.size();
  std::vector<int> hitIdx(nHits);
  std::iota(std::begin(hitIdx), std::end(hitIdx), 0);
  // sort hits to improve memory access
  std::sort(hitIdx.begin(), hitIdx.end(), [&hits](auto lhs, auto rhs) {
    return hits[lhs].getDetectorID() < hits[rhs].getDetectorID();
  });
  for (int i : hitIdx | std::views::filter([&](int idx) {
                 if (layer < 0)
                   return true;
                 return detID2Layer(hits[idx].getDetectorID()) == layer;
               })) {
    processHit(hits[i]);
  }
  finalizeDigits();
  writeDigits();
}

void NA6PMuonSpecDigitizer::getHitLocalCoord(const NA6PMuonSpecHit& hit, double xyzLocS[3], double xyzLocE[3])
{
  auto modID = hit.getDetectorID();
  const int geometryIndex = detID2GeometryIndex(modID);
  auto& matrix = mGeoManager.getMatrix(geometryIndex);

  double xyzGloS[3] = {hit.getXIn(), hit.getYIn(), hit.getZIn()};
  double xyzGloE[3] = {hit.getXOut(), hit.getYOut(), hit.getZOut()};

  matrix.MasterToLocal(xyzGloS, xyzLocS);
  matrix.MasterToLocal(xyzGloE, xyzLocE);
}

void NA6PMuonSpecDigitizer::processHit(const NA6PMuonSpecHit& hit)
{
  auto modID = hit.getDetectorID();
  auto& mod = mModules[modID];
  if (mod.isDisabled()) {
    LOGP(info, "Skipping disabled module {}", modID);
    return;
  }
  double xyzLocS[3], xyzLocE[3];
  getHitLocalCoord(hit, xyzLocS, xyzLocE);
  const int geometryIndex = detID2GeometryIndex(modID);
  const auto& mwpc = NA6PMWPCParam::Instance();
  const auto& layout = NA6PLayoutParam::Instance();
  int station = 0;
  int first = 0;
  for (; station < layout.nMSPlanes; ++station) {
    const int n = mwpc.stationGridNX[station] * mwpc.stationGridNY[station];
    if (modID < first + n) break;
    first += n;
  }
  mSegmentation.configure(mwpc, station);
  mSegmentation.setModuleHalfSize(mGeoManager.getModuleHalfX(geometryIndex),
                                  mGeoManager.getModuleHalfY(geometryIndex));
  double deltaX = xyzLocE[0] - xyzLocS[0];
  double deltaY = xyzLocE[1] - xyzLocS[1];
  const float pitch0 = mSegmentation.getStripPitch(0);
  const float pitch1 = mSegmentation.getStripPitch(1);
  const float minPitch = (pitch0 > 0.f && pitch1 > 0.f) ? std::min(pitch0, pitch1) : 1.f;
  int nSteps = std::max(std::abs(deltaX), std::abs(deltaY)) / minPitch;
  if (nSteps < 1)
    nSteps = 1;
  if (nSteps > 100)
    nSteps = 100;
  float chargePerStep = hit.getHitValue() * kGeVToEl / nSteps;
  float stepX = deltaX / nSteps;
  float stepY = deltaY / nSteps;
  float x = static_cast<float>(xyzLocS[0]) + 0.5f * stepX;
  float y = static_cast<float>(xyzLocS[1]) + 0.5f * stepY;
  for (int iStep = 0; iStep < nSteps; ++iStep) {
    NA6PMCComposedLabel lbl(hit.getTrackID(), mEventID, 0);
    // Future detector-response development: instead of assigning the full
    // step charge to one strip, distribute the induced charge over nearby
    // strips using a calibrated response (e.g. a Mathieson function), together
    // with realistic noise and threshold
    for (uint16_t readout = 0; readout < NA6PMuonSpecSegmentation::kNReadoutSystems; ++readout) {
      uint32_t strip = 0;
      if (mSegmentation.localToStrip(x, y, readout, strip)) {
        const auto key = mod.getOrderingKey(readout, strip);
        MSPreDigit* pd = mod.findDigit(key);
        if (!pd) {
          mod.addDigit(key, readout, strip, chargePerStep, lbl);
        } else {
          pd->addContribution(chargePerStep, lbl);
        }
      }
    }
    x += stepX;
    y += stepY;
  }
}

void NA6PMuonSpecDigitizer::finalizeDigits()
{
  for (int jMod = 0; jMod < mNumberOfModules; ++jMod) {
    auto& mod = mModules[jMod];
    if (mod.isDisabled()) {
      LOGP(info, "Skipping disabled module {}", jMod);
      continue;
    }
    auto& buffer = mod.getPreDigits();
    if (buffer.empty()) {
      continue;
    }
    auto itBeg = buffer.begin();
    auto iter = itBeg;
    for (; iter != buffer.end(); ++iter) {
      auto& preDig = iter->second;
      if (preDig.charge >= mThresholds[jMod]) {
        preDig.sortLabelsByEnergy();
        int digID = mDigits.size();
        mDigits.emplace_back(static_cast<uint16_t>(jMod), preDig.stripID,
                             preDig.charge);
        for (auto& plab : preDig.labels) {
          // here we can add selections on the labels to be stored in case of multiple labels per digit
          // e.g. remove delta electrons, remove particles contributing with much smaller energy than the others, etc.
          mMCLabels.addElement(digID, plab.second);
        }
      }
    }
    buffer.erase(itBeg, iter); // erase processed entries; iter == end() for now
  }
}
