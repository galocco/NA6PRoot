// NA6PCCopyright

#include "NA6PMuonSpecSegmentation.h"
#include "NA6PMWPCParam.h"
#include <cmath>

namespace { constexpr float kPi = 3.14159265358979323846f; }

void NA6PMuonSpecSegmentation::configure(const NA6PMWPCParam& param, int station)
{
  for (uint16_t readout = 0; readout < kNReadoutSystems; ++readout) {
    const bool hasStationPitch = station >= 0 && station < NA6PMWPCParam::MaxStations &&
                                 param.stationStripPitch[station][readout] > 0.f;
    setStripPitch(readout, hasStationPitch ? param.stationStripPitch[station][readout]
                                           : param.stripPitch[readout]);
    setStripAngleDeg(readout, param.stripAngleDeg[readout]);
    setNumberOfStrips(readout, param.stripCount[readout] > 0
                                      ? static_cast<uint32_t>(param.stripCount[readout])
                                      : 0u);
  }
}

void NA6PMuonSpecSegmentation::setStripPitch(uint16_t r, float v) { if (validReadout(r)) mPitch[r] = v; }
void NA6PMuonSpecSegmentation::setStripAngleDeg(uint16_t r, float v) { if (validReadout(r)) mAngleDeg[r] = v; }
void NA6PMuonSpecSegmentation::setStripOrigin(uint16_t r, float v) { if (validReadout(r)) { mOrigin[r] = v; mOriginSet[r] = true; } }
void NA6PMuonSpecSegmentation::setNumberOfStrips(uint16_t r, uint32_t v) { if (validReadout(r)) mNStrips[r] = v; }
float NA6PMuonSpecSegmentation::getStripPitch(uint16_t r) const { return validReadout(r) ? mPitch[r] : 0.f; }
float NA6PMuonSpecSegmentation::getStripAngleDeg(uint16_t r) const { return validReadout(r) ? mAngleDeg[r] : 0.f; }

bool NA6PMuonSpecSegmentation::localToStrip(float x, float y, uint16_t r, uint32_t& strip) const
{
  if (!validReadout(r) || mPitch[r] <= 0.f || std::abs(x) > mHalfX || std::abs(y) > mHalfY) return false;
  const float a = mAngleDeg[r] * kPi / 180.f;
  const float u = x * std::cos(a) + y * std::sin(a);
  const float origin = mOriginSet[r] ? mOrigin[r] :
    (-mHalfX * std::abs(std::cos(a)) - mHalfY * std::abs(std::sin(a)));
  const auto i = static_cast<int64_t>(std::floor((u - origin) / mPitch[r]));
  if (i < 0 || (mNStrips[r] && i >= static_cast<int64_t>(mNStrips[r]))) return false;
  strip = static_cast<uint32_t>(i);
  return true;
}

bool NA6PMuonSpecSegmentation::stripToLocalProjection(uint16_t r, uint32_t strip, float& u) const
{
  if (!validReadout(r) || mPitch[r] <= 0.f || (mNStrips[r] && strip >= mNStrips[r])) return false;
  const float a = mAngleDeg[r] * kPi / 180.f;
  const float origin = mOriginSet[r] ? mOrigin[r] :
    (-mHalfX * std::abs(std::cos(a)) - mHalfY * std::abs(std::sin(a)));
  u = origin + (static_cast<float>(strip) + 0.5f) * mPitch[r];
  return true;
}
