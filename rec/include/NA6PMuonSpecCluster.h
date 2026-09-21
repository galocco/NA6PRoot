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

#ifndef NA6P_MUONSPEC_CLUSTER_H
#define NA6P_MUONSPEC_CLUSTER_H

#include "NA6PBaseCluster.h"
#include <array>
#include <cstdint>

// Muon Spectrometer cluster class

class NA6PMuonSpecCluster : public NA6PBaseCluster
{
 public:

  NA6PMuonSpecCluster() = default;
  NA6PMuonSpecCluster(float x, float y, float z, int clusiz, int layer);
  NA6PMuonSpecCluster(const NA6PMuonSpecCluster&) = default;
  NA6PMuonSpecCluster& operator=(const NA6PMuonSpecCluster&) = default;

  float getReadoutCharge(uint16_t readout) const { return readout < 2 ? mReadoutCharge[readout] : 0.f; }
  uint16_t getReadoutSize(uint16_t readout) const { return readout < 2 ? mReadoutSize[readout] : 0; }
  void setReadoutCharge(uint16_t readout, float charge) { if (readout < 2) mReadoutCharge[readout] = charge; }
  void setReadoutSize(uint16_t readout, uint16_t size) { if (readout < 2) mReadoutSize[readout] = size; }

 private:
  std::array<float, 2> mReadoutCharge{};
  std::array<uint16_t, 2> mReadoutSize{};
  
  ClassDefNV(NA6PMuonSpecCluster, 2);
};

#endif
