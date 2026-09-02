// NA6PCCopyright

#ifndef NA6P_TOF_CLUSTER_H
#define NA6P_TOF_CLUSTER_H

#include "NA6PBaseCluster.h"

class NA6PTOFCluster : public NA6PBaseCluster
{
 public:
  NA6PTOFCluster() = default;
  NA6PTOFCluster(float x, float y, float z, int clusiz, int layer, float time);
  NA6PTOFCluster(const NA6PTOFCluster&) = default;
  NA6PTOFCluster& operator=(const NA6PTOFCluster&) = default;

  float getTime() const { return mTime; }
  void setTime(float time) { mTime = time; }

 private:
  float mTime = 0.f; // measured arrival time [s]

  ClassDefNV(NA6PTOFCluster, 1);
};

#endif
