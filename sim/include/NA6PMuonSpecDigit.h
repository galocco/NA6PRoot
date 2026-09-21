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

#ifndef NA6P_MUONSPEC_DIGIT_H
#define NA6P_MUONSPEC_DIGIT_H

#include <Rtypes.h>

// One analog strip digit from an MWPC chamber.  There are no tile/RSU/module
// subdivisions; the chamber itself is identified by detectorID.

struct MSStripID {
  uint16_t mReadout = 0; // one of the two strip-coordinate systems
  uint32_t mStrip = 0;   // strip number (count and pitch are TBD)
  uint16_t getReadout() const { return mReadout; }
  uint32_t getStrip() const { return mStrip; }
  void setReadout(uint16_t v) { mReadout = v; }
  void setStrip(uint32_t v) { mStrip = v; }
  ClassDefNV(MSStripID, 1);
};

class NA6PMuonSpecDigit
{
 public:
  NA6PMuonSpecDigit() = default;
  static constexpr uint16_t kNReadoutSystems = 2;
  NA6PMuonSpecDigit(uint16_t detID, const MSStripID& id, float charge);
  NA6PMuonSpecDigit(uint16_t detID, uint16_t readout, uint32_t strip, float charge);

  uint16_t getDetectorID() const { return mDetectorID; }
  const MSStripID& getStripID() const { return mStripID; }
  uint16_t getReadout() const { return mStripID.getReadout(); }
  uint32_t getStrip() const { return mStripID.getStrip(); }
  float getCharge() const { return mCharge; }

  void setDetectorID(uint16_t id) { mDetectorID = id; }
  void setStripID(const MSStripID& id) { mStripID = id; }
  void setReadout(uint16_t id) { mStripID.setReadout(id); }
  void setStrip(uint32_t id) { mStripID.setStrip(id); }
  void setCharge(float charge) { mCharge = charge; }

  void print() const;
  std::string asString() const;

 protected:
  uint16_t mDetectorID = 0; // the detector/sensor id
  MSStripID mStripID;           // strip identifier
  float mCharge = 0.f;          // analog amplitude (calibration units TBD)

  ClassDefNV(NA6PMuonSpecDigit, 1);
};

#endif
