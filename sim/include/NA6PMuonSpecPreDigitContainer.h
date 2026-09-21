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

#ifndef NA6P_MUONSPEC_PREDIGITCONTAINER_H
#define NA6P_MUONSPEC_PREDIGITCONTAINER_H

#include <Rtypes.h>
#include <algorithm>
#include <cstdint>
#include <map>
#include <utility>
#include <vector>
#include "NA6PMuonSpecDigit.h"
#include "NA6PMCComposedLabel.h"

// Accumulated analog signal on one MWPC strip before final digit creation.
struct MSPreDigit {
  MSStripID stripID;
  float charge = 0.f;
  std::vector<std::pair<float, NA6PMCComposedLabel>> labels;

  MSPreDigit(uint16_t readout = 0, uint32_t strip = 0, float value = 0.f,
             NA6PMCComposedLabel label = {})
    : charge(value)
  {
    stripID.setReadout(readout);
    stripID.setStrip(strip);
    labels.emplace_back(value, label);
  }

  void addContribution(float value, const NA6PMCComposedLabel& label)
  {
    charge += value;
    labels.emplace_back(value, label);
  }

  void sortLabelsByEnergy()
  {
    std::sort(labels.begin(), labels.end(),
              [](const auto& a, const auto& b) {
                return a.first > b.first;
              });
  }

  ClassDefNV(MSPreDigit, 1);
};

class NA6PMuonSpecPreDigitContainer
{
 public:
  NA6PMuonSpecPreDigitContainer(UShort_t idx = 0) : mDetectorIndex(idx){};
  ~NA6PMuonSpecPreDigitContainer() = default;

  std::map<ULong64_t, MSPreDigit>& getPreDigits() { return mPreDigits; }
  MSPreDigit* findDigit(ULong64_t key);
  void addDigit(ULong64_t key, uint16_t readout, uint32_t strip,
                float charge, const NA6PMCComposedLabel& label);

  UShort_t getDetectorIndex() const { return mDetectorIndex; }
  void setDetectorIndex(UShort_t ind) { mDetectorIndex = ind; }

  bool isEmpty() const { return mPreDigits.empty(); }
  static ULong64_t getOrderingKey(uint16_t readout, uint32_t strip)
  {
    return (static_cast<ULong64_t>(readout) << 32) | strip;
  }

  bool isDisabled() const { return mDisabled; }
  void disable(bool v) { mDisabled = v; }

 protected:
  UShort_t mDetectorIndex = 0;              ///< Detector index
  bool mDisabled = false;                   ///< flag to disable module
  std::map<ULong64_t, MSPreDigit> mPreDigits;

  ClassDefNV(NA6PMuonSpecPreDigitContainer, 1);
};

//_______________________________________________________________________
inline MSPreDigit* NA6PMuonSpecPreDigitContainer::findDigit(ULong64_t key)
{
  // finds the digit corresponding to global key
  auto digitentry = mPreDigits.find(key);
  return digitentry != mPreDigits.end() ? &(digitentry->second) : nullptr;
}
//_______________________________________________________________________
inline void NA6PMuonSpecPreDigitContainer::addDigit(ULong64_t key, uint16_t readout,
                                                    uint32_t strip, float charge,
                                                    const NA6PMCComposedLabel& label)
{
  mPreDigits.emplace(std::make_pair(key, MSPreDigit(readout, strip, charge, label)));
}

#endif
