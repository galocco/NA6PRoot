// NA6PCCopyright

#ifndef NA6P_MWPC_TEST_MODULE_H_
#define NA6P_MWPC_TEST_MODULE_H_

#include "NA6PModule.h"

#include <cstddef>

class TGeoVolume;

// Minimal active module used only by the isolated MWPC transport test.
// It deliberately reuses NA6PMWPCChamber instead of duplicating the chamber geometry.
class MWPCTestModule final : public NA6PModule
{
 public:
  struct EventStats {
    int gasSteps = 0;
    int gasEntries = 0;
    int gasExits = 0;
    int primarySteps = 0;
    int primaryEntries = 0;
    double energyDepositGeV = 0.;
  };

  explicit MWPCTestModule(int activeID);
  ~MWPCTestModule() override = default;

  void createMaterials() override;
  void createGeometry(TGeoVolume* world) override;
  bool stepManager(int volID) override;

  size_t getNHits() const override { return static_cast<size_t>(mStats.gasEntries); }
  void clearHits() override { mStats = {}; }

  // The test writes its own compact CSV, so no detector hit tree is needed here.
  void createHitsOutput(const std::string&) override {}
  void closeHitsOutput() override {}
  void writeHits(const std::vector<int>&) override {}

  const EventStats& stats() const { return mStats; }

 private:
  EventStats mStats{};
};

#endif
