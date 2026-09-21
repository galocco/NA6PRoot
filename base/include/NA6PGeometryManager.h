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

#ifndef NA6P_GEOMETRY_MANAGER_H
#define NA6P_GEOMETRY_MANAGER_H

#include <Rtypes.h>
#include <TGeoMatrix.h>
#include <vector>

class TGeoVolume;

// Helper class to access cached geometry matrices and sizes of alignable sensors.

class NA6PGeometryManager
{
 public:
  static constexpr int kNVTModulesPerLayer = 4;

  static int getMuonSpecGeometryIndex(int detectorID);

  bool loadGeometry(const char* filename = "geometry.root", const char* geoname = "NA6P");

  // Integer access uses the TGeo alignable-entry index. This is also the
  // VerTel detector ID as long as the VT entries are registered first and in
  // detector-ID order.
  const TGeoHMatrix& getMatrix(int jMod) const { return mMatrices[jMod]; }
  float getModuleHalfX(int jMod) const { return mModuleHalfX[jMod]; }
  float getModuleHalfY(int jMod) const { return mModuleHalfY[jMod]; }

 private:
  bool fillModuleSize(int jMod, TGeoVolume* vol);

  std::vector<TGeoHMatrix> mMatrices;   ///< Local-to-global transforms
  std::vector<float> mModuleHalfX;      ///< Module half length along x
  std::vector<float> mModuleHalfY;      ///< Module half length along y
  bool mGeoLoaded{false};               ///< Flag for successful geometry load
  ClassDefNV(NA6PGeometryManager, 1);
};

#endif
