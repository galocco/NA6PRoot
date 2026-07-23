// NA6PCCopyright

#include <fmt/format.h>
#include <fairlogger/Logger.h>
#include "NA6PVerTelCluster.h"
#include "NA6PMuonSpecCluster.h"
#include "NA6PTrack.h"
#include "NA6PLayoutParam.h"

//_______________________________________________________________________
NA6PTrack::NA6PTrack()
{
  mClusterIndices.fill(-1);
  invalidate();
  mOuter.invalidate();
  mConstrained.invalidate();
}

//_______________________________________________________________________
NA6PTrack::NA6PTrack(const float* xyz, const float* pxyz, int sign, float errLoose)
{
  mClusterIndices.fill(-1);
  init(xyz, pxyz, sign, errLoose);
  mOuter.invalidate();
  mConstrained.invalidate();
}

//_______________________________________________________________________
void NA6PTrack::reset()
{
  mChi2 = 0;
  mChi2Outer = 0;
  resetCovariance();
  resetClusters();
  mStatusMS = kNotMS;
  invalidate();
  mOuter.invalidate();
  mConstrained.invalidate();
}
//_______________________________________
void NA6PTrack::resetClusters()
{
  mNClusters = 0;
  mClusterMap = 0;
  mClusterIndices.fill(-1);
}

//_______________________________________________________________________
std::string NA6PTrack::asString() const
{
  auto pxyz = getPXYZ();
  return fmt::format("Track: Nclusters:{} chi2:{} chi2:{} pos:{:.4f},{:.4f},{:.4f} mom:{:.3f},{:.3f},{:.3f}",
                     mNClusters, mChi2, mChi2Outer, getX(), getY(), getZ(), pxyz[0], pxyz[1], pxyz[2]);
}

//_______________________________________________________________________
void NA6PTrack::print() const
{
  LOGP(info, "{}", asString());
}

//_______________________________________
template <typename ClusterType>
void NA6PTrack::addCluster(const ClusterType* clu)
{

  mNClusters++;
  int nLay = clu->getLayer();
  mClusterIndices[nLay] = clu->getClusterIndex();
  mClusterMap |= (1 << nLay);
}

template void NA6PTrack::addCluster<NA6PBaseCluster>(const NA6PBaseCluster*);
template void NA6PTrack::addCluster<NA6PVerTelCluster>(const NA6PVerTelCluster*);
template void NA6PTrack::addCluster<NA6PMuonSpecCluster>(const NA6PMuonSpecCluster*);
