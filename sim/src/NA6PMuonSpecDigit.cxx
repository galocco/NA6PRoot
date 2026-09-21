// NA6PCCopyright

#include "NA6PMuonSpecDigit.h"
#include <fmt/format.h>
#include <fairlogger/Logger.h>

NA6PMuonSpecDigit::NA6PMuonSpecDigit(uint16_t detID, const MSStripID& id, float charge)
  : mDetectorID(detID), mStripID(id), mCharge(charge)
{
}

NA6PMuonSpecDigit::NA6PMuonSpecDigit(uint16_t detID, uint16_t readout, uint32_t strip, float charge)
  : mDetectorID(detID), mCharge(charge)
{
  mStripID.setReadout(readout);
  mStripID.setStrip(strip);
}

std::string NA6PMuonSpecDigit::asString() const
{
  return fmt::format("Digit: Det:{} Readout:{} Strip:{} Charge:{}",
                     mDetectorID, getReadout(), getStrip(), mCharge);
}

void NA6PMuonSpecDigit::print() const
{
  LOGP(info, "{}", asString());
}
