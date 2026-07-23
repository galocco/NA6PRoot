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

#ifndef NA6P_MATCH_H
#define NA6P_MATCH_H

#include <string>
#include <Rtypes.h>
#include "NA6PTrackParCov.h"

class NA6PMatch : public NA6PTrackParCov
{
 public:
  using NA6PTrackParCov::NA6PTrackParCov;

  void setChi2Refit(float c) { mChi2Refit = c; }
  float getChi2Refit() const { return mChi2Refit; }
  float getChi2Norm(int ndof = 5) const { return mChi2Refit / std::max(1, 2 * mNClusters - ndof); }

  void setChi2Match(float c) { mChi2Match = c; }
  float getChi2Match() const { return mChi2Match; }

  int getNClusters() const { return mNClusters; }
  void setNClusters(int n) { mNClusters = n; }

  int getIndexVT() const { return mIndexVT; }
  int getIndexMS() const { return mIndexMS; }
  void setIndexVT(int i) { mIndexVT = i; }
  void setIndexMS(int i) { mIndexMS = i; }

  void print() const;
  std::string asString() const;

 protected:
  float mChi2Match = 0.f; // matching chi2
  float mChi2Refit = 0.f; // total chi2 of refit
  int mIndexVT = -1;      // index of the VT Track
  int mIndexMS = -1;      // index of the MS Track
  int8_t mNClusters = 0;  // total number of clusters

  ClassDefNV(NA6PMatch, 2)
};

#endif
