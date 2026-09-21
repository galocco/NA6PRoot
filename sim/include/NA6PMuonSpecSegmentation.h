// NA6PCCopyright

#ifndef NA6P_MUONSPEC_SEGMENTATION_H
#define NA6P_MUONSPEC_SEGMENTATION_H

#include <Rtypes.h>
#include <cstdint>

struct NA6PMWPCParam;

// Mapping from chamber-local coordinates (cm) to the two MWPC strip systems.
// It contains no module subdivisions or dead-zone model.
class NA6PMuonSpecSegmentation
{
 public:
  static constexpr uint16_t kNReadoutSystems = 2;
  NA6PMuonSpecSegmentation() = default;
  void configure(const NA6PMWPCParam& param, int station = -1);
  void setModuleHalfSize(float halfX, float halfY) { mHalfX = halfX; mHalfY = halfY; }
  void setStripPitch(uint16_t readout, float pitch);
  void setStripAngleDeg(uint16_t readout, float angle);
  void setStripOrigin(uint16_t readout, float origin);
  void setNumberOfStrips(uint16_t readout, uint32_t count);
  float getStripPitch(uint16_t readout) const;
  float getStripAngleDeg(uint16_t readout) const;
  bool localToStrip(float x, float y, uint16_t readout, uint32_t& strip) const;
  bool stripToLocalProjection(uint16_t readout, uint32_t strip, float& u) const;

 private:
  bool validReadout(uint16_t readout) const { return readout < kNReadoutSystems; }
  float mHalfX = 0.f;
  float mHalfY = 0.f;
  // Temporary engineering default (1 mm). It is not a final detector
  // parameter and should be replaced when the strip design is defined.
  float mPitch[kNReadoutSystems] = {0.1f, 0.1f};
  float mAngleDeg[kNReadoutSystems] = {0.f, 90.f};
  float mOrigin[kNReadoutSystems] = {0.f, 0.f};
  bool mOriginSet[kNReadoutSystems] = {false, false};
  uint32_t mNStrips[kNReadoutSystems] = {0, 0};
  ClassDefNV(NA6PMuonSpecSegmentation, 1);
};

#endif
