// NA6PCCopyright

#include "NA6PTOFMatch.h"

#include <fairlogger/Logger.h>
#include <fmt/format.h>

std::string NA6PTOFMatch::asString() const
{
  return fmt::format("TOFMatch: trackID:{} TOFid:{} time:{:.3e} s length:{:.3f} cm beta:{:.4f} chi2:{:.4f} {}",
                     mIndexTrack, mIndexTOF, mTOF, mPathLength, getBeta(), mMatchChi2,
                     static_cast<const NA6PTrackParCov*>(this)->asString());
}

void NA6PTOFMatch::print() const
{
  LOGP(info, "{}", asString());
}
