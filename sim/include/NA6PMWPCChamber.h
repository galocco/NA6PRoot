// NA6PCCopyright

#ifndef NA6P_MWPC_CHAMBER_H_
#define NA6P_MWPC_CHAMBER_H_

#include <array>
#include <string>
#include <vector>

class NA6PModule;
class TGeoVolume;
class TGeoVolumeAssembly;

class NA6PMWPCChamber
{
 public:
  struct Materials {
    std::string fr4;
    std::string copper;
    std::string honeycomb;
    std::string gas;
  };

  struct Placement {
    double x = 0.;
    double y = 0.;
    double z = 0.;
  };

  // Negative body overrides mean: use the values from NA6PMWPCParam.
  // The override is useful for the present narrow-MS0 study while keeping
  // exactly the same chamber construction code for all six stations.
  NA6PMWPCChamber(const NA6PModule& module, Materials materials,
                  double bodyXOverride = -1., double bodyYOverride = -1.);

  void createMaterials() const;
  TGeoVolumeAssembly* addTo(TGeoVolume* parent, int chamberID, int localCopyID, const Placement& placement) const;

  std::array<double, 3> envelopeFullSize() const;
  std::array<double, 3> gasFullSize() const;
  std::array<double, 3> gasCentre() const;

 private:
  enum class MaterialKind { FR4, Copper, Honeycomb, Gas };

  struct Part {
    std::string name;
    MaterialKind material;
    std::array<double, 3> fullSize;
    std::array<double, 3> centre;
    bool sensitive = false;
  };

  std::vector<Part> buildParts() const;
  void validate(const std::vector<Part>& parts) const;
  static bool overlaps(const Part& a, const Part& b);
  static void requirePositive(double value, const char* name, bool zeroAllowed = false);

  const NA6PModule& mModule;
  Materials mMaterials;
  double mBodyXOverride = -1.;
  double mBodyYOverride = -1.;
};

#endif
