// NA6PCCopyright

#ifndef NA6P_MUONSPECMODULAR_H_
#define NA6P_MUONSPECMODULAR_H_

#include "NA6PModule.h"
#include "NA6PMuonSpecModularHit.h"

class TFile;
class TTree;
class TGeoVolume;

class NA6PMuonSpecModular : public NA6PModule
{
 public:
  NA6PMuonSpecModular() : NA6PModule("MuonSpecModular") { setActiveID(1); }
  ~NA6PMuonSpecModular() override = default;
  void createMaterials() override;
  void createGeometry(TGeoVolume* world) override;
  bool stepManager(int volID) override;
  size_t getNHits() const override { return mHits.size(); }
  void createHitsOutput(const std::string& outDir) override;
  void closeHitsOutput() override;
  void writeHits(const std::vector<int>& remapping) override;
  void setAlignableEntries() override;

  void clearHits() override { mHits.clear(); }

  const auto& getHits() const { return mHits; }

  NA6PMuonSpecModularHit* addHit(int trackID, int detID, const TVector3& startPos, const TVector3& endPos, const TVector3& startMom, const TVector3& endMom,
                                 float endTime, float eLoss, unsigned char startStatus, unsigned char endStatus);

 private:
  std::vector<NA6PMuonSpecModularHit> mHits, *hHitsPtr = &mHits;
  TFile* mHitsFile = nullptr;
  TTree* mHitsTree = nullptr;
};

#endif
