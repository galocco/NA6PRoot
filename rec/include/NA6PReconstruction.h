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

#ifndef NA6P_RECONSTRUCTION_H
#define NA6P_RECONSTRUCTION_H

#include <string>
#include <Rtypes.h>

// Class to steer the reconstruction
class NA6PVertex;
class NA6PTrack;
class NA6PMCTruthContainer;
class NA6PMCComposedLabel;

class NA6PReconstruction
{
 public:
  NA6PReconstruction(const std::string& name) : mName(name) {}
  virtual ~NA6PReconstruction() = default;

  const std::string& getName() const { return mName; }
  /*
  RSREM RecoParam is a singleton. If one can have multiple NA6PReconstruction-based objects, they should not set each one its own recoparam,
  since they will interfere with eash other. Instead, the recoparam must be configured upstream.
  The same is true for the geometry and the field: these are static objects, one 1 copy of each must be initialized, do this upstream.
  virtual bool init(const std::string& filename, const std::string& geoname = "NA6P");
  void setRecoParamFile(const std::string& recparfile) { mRecoParFilName = recparfile; } // RSTODO is this needed here? Init upstream
  void setGeometryFile(const std::string& geofile, const std::string& geoname = "NA6P") // RSTODO is this needed here? Init upstream
  {
    mGeoFilName = geofile;
    mGeoObjName = geoname;
  }
  */

  // methods to steer cluster reconstruction
  virtual void createClustersOutput();
  virtual void clearClusters();
  virtual void writeClusters();
  virtual void closeClustersOutput();
  // methods to steer vertexing
  virtual void createVerticesOutput();
  virtual void clearVertices();
  virtual void writeVertices();
  virtual void closeVerticesOutput();
  // methods to steer tracking
  virtual void createTracksOutput();
  virtual void clearTracks();
  virtual void writeTracks();
  virtual void closeTracksOutput();

  void setPrimaryVertex(const NA6PVertex* v) { mPrimaryVertex = v; }
  virtual void clearEvent() { mPrimaryVertex = nullptr; }
  void setReadMCTruth(bool val = true) { mReadMCTruth = val; }
  void assignMCLabels(std::vector<NA6PTrack>& trk,
                      std::vector<NA6PMCComposedLabel>& mcTrkLabels,
                      const NA6PMCTruthContainer& mcCluLabels);

 protected:
  std::string mName{"MothClass"}; // detector name
  /*
  // RSREM : see comment above for getters
  std::string mGeoFilName{"geometry.root"};   // name of geometry file
  std::string mGeoObjName{"NA6P"};            // name of geometry object
  std::string mRecoParFilName{""};            // name of reco param file
  bool mIsInitialized = false;                // flag for initialization
  */
  const NA6PVertex* mPrimaryVertex = nullptr; // primary vertex
  bool mReadMCTruth = true;                   // flag to enable/disable the MC info

  ClassDefNV(NA6PReconstruction, 1);
};

#endif
