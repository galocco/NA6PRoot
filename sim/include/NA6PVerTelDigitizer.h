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

#ifndef NA6P_VERTEL_DIGITIZER_H
#define NA6P_VERTEL_DIGITIZER_H

#include "NA6PVerTelPreDigitContainer.h"
#include "NA6PVerTelHit.h"
#include "NA6PVerTelDigit.h"
#include "NA6PVerTelSegmentation.h"
#include "NA6PGeometryManager.h"
#include "NA6PMCTruthContainer.h"
#include <Rtypes.h>
#include <TGeoMatrix.h>

class TFile;
class TTree;

// steers hits -> digits step

class NA6PVerTelDigitizer
{
 public:
  static constexpr float kDefaultThresholdkeV = 0.4f;
  static constexpr float kDefaultThresholdEl = 100.f;
  static constexpr float kGeVTokeV = 1.e6;
  static constexpr float kGeVToEl = 2.76e8;

  NA6PVerTelDigitizer() = default;
  NA6PVerTelDigitizer(const std::string& name) : mName(name) {}
  ~NA6PVerTelDigitizer() = default;

  const std::string& getName() const { return mName; }

  static int detID2Layer(int detID) { return detID / NA6PGeometryManager::kNVTModulesPerLayer; }

  void init(const char* filename = "geometry.root", const char* geoname = "NA6P");
  void initGeometry(const char* filename = "geometry.root", const char* geoname = "NA6P");
  void process(const std::vector<NA6PVerTelHit>& hits, int layer = -1);
  void processHit(const NA6PVerTelHit& hit);
  void setEventMetaData(int evID) { mEventID = evID; }
  void finalizeDigits();

  size_t getNDigits() const { return mDigits.size(); }
  void createDigitsOutput();
  void closeDigitsOutput();
  void writeDigits();
  void clearDigits()
  {
    mDigits.clear();
    mMCLabels.clear_andfreememory();
  }
  const auto& getDigits() const { return mDigits; }

  void getHitLocalCoord(const NA6PVerTelHit& hit, double xyzLocS[3], double xyzLocE[3]);

  void SetThreshold(int modID, float thr)
  {
    if (modID < 0 || modID >= mNumberOfModules) {
      LOGP(error, "Module ID [] out of range", modID);
    } else {
      mThresholds[modID] = thr;
    }
  }

 protected:
  std::string mName{"VerTel"};                       ///< detector name
  int mNumberOfModules = 0;                          ///< number of modules
  std::vector<NA6PVerTelPreDigitContainer> mModules; ///< Array of module pre-digits containers
  NA6PVerTelSegmentation mSegmentation;              ///< segmentation class
  NA6PGeometryManager mGeoManager;                   ///< geometry manager
  std::vector<float> mThresholds;                    ///< Threshold (per tile)
  std::vector<NA6PVerTelDigit> mDigits, *hDigitsPtr = &mDigits;
  NA6PMCTruthContainer mMCLabels, *hMCLabelsPtr = &mMCLabels;
  TFile* mDigitsFile = nullptr;
  TTree* mDigitsTree = nullptr;
  // info about event beint digitized
  int mEventID = 0; // at the moment only event ID, later need at least timestamp for the pile-up

  ClassDefNV(NA6PVerTelDigitizer, 1);
};

#endif
