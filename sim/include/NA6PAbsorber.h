// NA6PCCopyright

#ifndef NA6P_ABSORBER_H_
#define NA6P_ABSORBER_H_

#include "NA6PModule.h"

class TGeoVolume;

class NA6PAbsorber : public NA6PModule
{
 public:
  static constexpr int PlugCopyOffset = 200;

  NA6PAbsorber() : NA6PModule("Absorber") {}
  ~NA6PAbsorber() override = default;
  void createMaterials() override;
  void createGeometry(TGeoVolume* world) override;
};

#endif
