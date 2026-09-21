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

#ifndef NA6P_MUONSPEC_DIGITIZER_H
#define NA6P_MUONSPEC_DIGITIZER_H

#include "NA6PMuonSpecPreDigitContainer.h"
#include "NA6PMuonSpecHit.h"
#include "NA6PMuonSpecDigit.h"
#include "NA6PMuonSpecSegmentation.h"
#include "NA6PGeometryManager.h"
#include "NA6PMWPCParam.h"
#include "NA6PLayoutParam.h"
#include "NA6PMCTruthContainer.h"
#include <Rtypes.h>
#include <TGeoMatrix.h>

class TFile;
class TTree;

// steers hits -> digits step

class NA6PMuonSpecDigitizer
{
 public:
  static constexpr float kDefaultThresholdkeV = 0.084f;

  static constexpr float kGeVTokeV = 1.e6f;
  static constexpr float kGeVToEl = 3.57e7f; // W = 28 eV, Ar/CO2
  static constexpr float kDefaultThresholdEl = kDefaultThresholdkeV * kGeVTokeV / kGeVToEl;

  NA6PMuonSpecDigitizer() = default;
  NA6PMuonSpecDigitizer(const std::string& name) : mName(name) {}
  ~NA6PMuonSpecDigitizer() = default;

  const std::string& getName() const { return mName; }

  // Hit detector IDs are chamber IDs, sequential across the MS stations.
  // Convert them to the global layer namespace used by reconstruction, where
  // the MS layers follow the Vertex Telescope layers.
  static int detID2Layer(int detID)
  {
    if (detID < 0) {
      return -1;
    }
    const auto& layout = NA6PLayoutParam::Instance();
    const auto& mwpc = NA6PMWPCParam::Instance();
    int first = 0;
    for (int station = 0; station < layout.nMSPlanes; ++station) {
      const int nModules = mwpc.stationGridNX[station] * mwpc.stationGridNY[station];
      if (detID < first + nModules) {
        return layout.nVerTelPlanes + station;
      }
      first += nModules;
    }
    return -1;
  }

  // MuonSpec detector IDs are chamber-local (0..N-1), while the geometry
  // cache stores MWPC entries after all Vertex Telescope modules.
  static int detID2GeometryIndex(int detID)
  {
    return NA6PGeometryManager::getMuonSpecGeometryIndex(detID);
  }

  void init(const char* filename = "geometry.root", const char* geoname = "NA6P");
  void initGeometry(const char* filename = "geometry.root", const char* geoname = "NA6P");
  void process(const std::vector<NA6PMuonSpecHit>& hits, int layer = -1);
  void processHit(const NA6PMuonSpecHit& hit);
  void setEventMetaData(int evID) { mEventID = evID; }
  void finalizeDigits();

  size_t getNDigits() const { return mDigits.size(); }
  void createDigitsOutput(const std::string& outDir = "");
  void closeDigitsOutput();
  void writeDigits();
  void clearDigits()
  {
    mDigits.clear();
    mMCLabels.clear_andfreememory();
  }
  const auto& getDigits() const { return mDigits; }

  void getHitLocalCoord(const NA6PMuonSpecHit& hit, double xyzLocS[3], double xyzLocE[3]);

  void SetThreshold(int modID, float thr)
  {
    if (modID < 0 || modID >= mNumberOfModules) {
      LOGP(error, "Module ID [] out of range", modID);
    } else {
      mThresholds[modID] = thr;
    }
  }

 protected:
  std::string mName{"MuonSpec"};                       ///< detector name
  int mNumberOfModules = 0;                            ///< number of modules
  std::vector<NA6PMuonSpecPreDigitContainer> mModules; ///< Array of module pre-digits containers
  NA6PMuonSpecSegmentation mSegmentation;              ///< segmentation class
  NA6PGeometryManager mGeoManager;                     ///< geometry manager
  std::vector<float> mThresholds;                      ///< Threshold (per module)
  std::vector<NA6PMuonSpecDigit> mDigits, *hDigitsPtr = &mDigits;
  NA6PMCTruthContainer mMCLabels, *hMCLabelsPtr = &mMCLabels;
  TFile* mDigitsFile = nullptr;
  TTree* mDigitsTree = nullptr;
  // info about event being digitized
  int mEventID = 0; // at the moment only event ID, later need at least timestamp for the pile-up

  ClassDefNV(NA6PMuonSpecDigitizer, 1);
};

#endif
