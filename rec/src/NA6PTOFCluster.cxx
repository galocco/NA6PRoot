// NA6PCCopyright

#include "NA6PTOFCluster.h"

NA6PTOFCluster::NA6PTOFCluster(float x, float y, float z, int clusiz, int layer, float time)
  : NA6PBaseCluster(x, y, z, clusiz, layer), mTime(time)
{
}
