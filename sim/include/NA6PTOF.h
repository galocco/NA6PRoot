// NA6PCCopyright

#ifndef NA6P_TOF_H_
#define NA6P_TOF_H_

#include "NA6PModule.h"
#include "NA6PTOFHit.h"

#include <vector>

class TFile;
class TTree;
class TGeoVolume;

class NA6PTOF : public NA6PModule
{
 public:
  NA6PTOF() : NA6PModule("TOF") { setActiveID(2); }
  ~NA6PTOF() override = default;

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

  NA6PTOFHit* addHit(int trackID, int detID, const TVector3& startPos, const TVector3& endPos,
                     const TVector3& startMom, const TVector3& endMom, float endTime,
                     float eLoss, unsigned char startStatus, unsigned char endStatus);

 private:
  static constexpr int NTiles = 4;

  std::vector<NA6PTOFHit> mHits, *hHitsPtr = &mHits;
  TFile* mHitsFile = nullptr;
  TTree* mHitsTree = nullptr;
};

#endif
