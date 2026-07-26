// NA6PCCopyright

#ifndef NA6P_GENHEPMC_H
#define NA6P_GENHEPMC_H

#include "NA6PGenerator.h"
#include "HepMC3/ReaderRootTree.h"

class NA6PGenHepMC : public NA6PGenerator
{
 public:
  using NA6PGenerator::NA6PGenerator;

  NA6PGenHepMC() = default;
  NA6PGenHepMC(const std::string& name, const std::string& filename, bool storeDecayed = true);
  ~NA6PGenHepMC() override = default;
  void generate() override;

  void init() override;

  void SetFileName(const std::string& name) { mFileName = name; }
  long canGenerateMaxEvents() const override;
  void skipEvents(long nEvents) override;
  void TrackSpectators() { mTrackSpectators = true; }

 protected:
  inline static const std::string HEPTreeName{"hepmc3_tree"};
  std::string mFileName = {};
  bool mStoreDecayedPrimaries = true;
  bool mTrackSpectators = false;
  long mNEvInTree = -1;
  long mReadEvents = 0;
  std::unique_ptr<HepMC3::ReaderRootTree> mHEPRootFileReader;

  ClassDefOverride(NA6PGenHepMC, 1); //
};

#endif
