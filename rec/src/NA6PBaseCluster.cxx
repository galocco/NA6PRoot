// NA6PCCopyright

#include "NA6PBaseCluster.h"
#include <fmt/format.h>
#include <fairlogger/Logger.h>

NA6PBaseCluster::NA6PBaseCluster(float x, float y, float z, int clusiz, int layer)
  : mXYZ{x, y, z}, mCluSiz{clusiz}, mLayer{int8_t(layer)}
{
}

std::string NA6PBaseCluster::asString() const
{
  return fmt::format("Cluster: Det:{} XYZ:{:+.3e},{:+.3e},{:+.3e} Sig2:{:.3e},{:.3e},{:.3e} cluster size:{}",
                     mDetectorID, getX(), getY(), getZ(), getSigXX(), getSigYX(), getSigYY(), mCluSiz);
}

void NA6PBaseCluster::print() const
{
  LOGP(info, "{}", asString());
}
