// NA6PCCopyright

// Based on:
// Copyright 2019-2020 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

#ifndef NA6P_BASE_CLUSTER_H
#define NA6P_BASE_CLUSTER_H

#include <string>
#include <array>
#include <Rtypes.h>

// Mother class of all cluster classes

class NA6PBaseCluster
{
 public:
  NA6PBaseCluster() = default;
  NA6PBaseCluster(float x, float y, float z, int clusiz, int layer);
  NA6PBaseCluster(const NA6PBaseCluster&) = default;
  NA6PBaseCluster& operator=(const NA6PBaseCluster&) = default;
  ~NA6PBaseCluster() = default;

  int getLayer() const { return mLayer; }
  int getClusterSize() const { return mCluSiz; }
  auto getX() const { return mXYZ[0]; }
  auto getY() const { return mXYZ[1]; }
  auto getZ() const { return mXYZ[2]; }
  const auto& getXYZ() const { return mXYZ; }
  auto& getXYZ() { return mXYZ; }
  auto getSigXX() const { return mCov[0]; }
  auto getSigYX() const { return mCov[1]; }
  auto getSigYY() const { return mCov[2]; }
  const auto& getCov() const { return mCov; }
  auto& getCov() { return mCov; }

  auto getDetectorID() const { return mDetectorID; }
  auto getHitID() const { return mHitID; }
  int getClusterIndex() const { return mClusterIndex; }

  void setDetectorID(int id) { mDetectorID = id; }
  void setHitID(int id) { mHitID = id; }
  void setPos(float x, float y, float z) { mXYZ = {x, y, z}; }
  void setX(float v) { mXYZ[0] = v; }
  void setY(float v) { mXYZ[1] = v; }
  void setZ(float v) { mXYZ[2] = v; }
  void setErr(float sxx, float syx, float syy) { mCov = {sxx, syx, syy}; }
  void setClusterIndex(int idx) { mClusterIndex = idx; }

  void print() const;
  std::string asString() const;

 protected:
  std::array<float, 3> mXYZ{};
  std::array<float, 3> mCov{};

  int mCluSiz = 0;       // cluster size
  int mHitID = -1;       // hit ID (for test of hitsToRecPoints)
  short mDetectorID = 0; // the detector/sensor id
  int8_t mLayer = -1;
  int mClusterIndex = -1; //! Transient index for track-cluster association

  ClassDefNV(NA6PBaseCluster, 3);
};

#endif
