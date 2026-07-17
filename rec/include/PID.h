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

/// \file PID.h
/// \brief particle ids, masses, names class definition
/// \author ruben.shahoyan@cern.ch

#ifndef NA6P_TRACK_PID_H_
#define NA6P_TRACK_PID_H_

#include <cstdint>
#include <limits>

#include "Rtypes.h"
#include "NA6PPhysConst.h"

class PID
{
 public:
  // particle identifiers, continuos starting from 0
  typedef uint8_t ID;
  static constexpr ID Electron = 0;
  static constexpr ID Muon = 1;
  static constexpr ID Pion = 2;
  static constexpr ID Kaon = 3;
  static constexpr ID Proton = 4;
  static constexpr ID Deuteron = 5;
  static constexpr ID Triton = 6;
  static constexpr ID Helium3 = 7;
  static constexpr ID Alpha = 8;
  static constexpr ID Pi0 = 9;
  static constexpr ID Photon = 10;
  static constexpr ID K0 = 11;
  static constexpr ID Lambda = 12;
  static constexpr ID HyperTriton = 13;
  static constexpr ID Hyperhydrog4 = 14;
  static constexpr ID XiMinus = 15;
  static constexpr ID OmegaMinus = 16;
  static constexpr ID HyperHelium4 = 17;
  static constexpr ID HyperHelium5 = 18;

  static constexpr ID First = Electron;
  static constexpr ID Last = HyperHelium5; ///< if extra IDs added, update this !!!
  static constexpr ID NIDs = Last + 1;     ///< number of defined IDs
  static constexpr const char* sNames[NIDs + 1] = {"Electron", "Muon", "Pion", "Kaon", "Proton", "Deuteron", "Triton", "Helium3", "Alpha", "Pi0", "Photon", "K0", "Lambda", "HyperTriton", "Hyperhydrog4", "XiMinus", "OmegaMinus", "HyperHelium4", "HyperHelium5", nullptr};
  static constexpr float sMasses[NIDs] = {phys_const::MassElectron, phys_const::MassMuon, phys_const::MassPionCharged, phys_const::MassKaonCharged, phys_const::MassProton, phys_const::MassDeuteron, phys_const::MassTriton, phys_const::MassHelium3, phys_const::MassAlpha, phys_const::MassPi0, phys_const::MassPhoton, phys_const::MassK0, phys_const::MassLambda, phys_const::MassHyperTriton, phys_const::MassHyperhydrog4, phys_const::MassXiMinus, phys_const::MassOmegaMinus, phys_const::MassHyperHelium4, phys_const::MassHyperHelium5};
  static constexpr float sMasses2[NIDs] = {phys_const::MassElectron * phys_const::MassElectron, phys_const::MassMuon* phys_const::MassMuon, phys_const::MassPionCharged* phys_const::MassPionCharged, phys_const::MassKaonCharged* phys_const::MassKaonCharged, phys_const::MassProton* phys_const::MassProton, phys_const::MassDeuteron* phys_const::MassDeuteron, phys_const::MassTriton* phys_const::MassTriton, phys_const::MassHelium3* phys_const::MassHelium3, phys_const::MassAlpha* phys_const::MassAlpha, phys_const::MassPi0* phys_const::MassPi0, phys_const::MassPhoton* phys_const::MassPhoton, phys_const::MassK0* phys_const::MassK0, phys_const::MassLambda* phys_const::MassLambda, phys_const::MassHyperTriton* phys_const::MassHyperTriton, phys_const::MassHyperhydrog4* phys_const::MassHyperhydrog4, phys_const::MassXiMinus* phys_const::MassXiMinus, phys_const::MassOmegaMinus* phys_const::MassOmegaMinus, phys_const::MassHyperHelium4* phys_const::MassHyperHelium4, phys_const::MassHyperHelium5* phys_const::MassHyperHelium5};
  static constexpr float sMassesI[NIDs] = {1.f / phys_const::MassElectron, 1.f / phys_const::MassMuon, 1.f / phys_const::MassPionCharged, 1.f / phys_const::MassKaonCharged, 1.f / phys_const::MassProton, 1.f / phys_const::MassDeuteron, 1.f / phys_const::MassTriton, 1.f / phys_const::MassHelium3, 1.f / phys_const::MassAlpha, std::numeric_limits<float>::infinity(), std::numeric_limits<float>::infinity(), std::numeric_limits<float>::infinity(), std::numeric_limits<float>::infinity(), 1.f / phys_const::MassHyperTriton, 1.f / phys_const::MassHyperhydrog4, 1.f / phys_const::MassXiMinus, 1.f / phys_const::MassOmegaMinus, 1.f / phys_const::MassHyperHelium4, 1.f / phys_const::MassHyperHelium5};
  static constexpr float sMasses2Z[NIDs] = {
    phys_const::MassElectron,
    phys_const::MassMuon,
    phys_const::MassPionCharged,
    phys_const::MassKaonCharged,
    phys_const::MassProton,
    phys_const::MassDeuteron,
    phys_const::MassTriton,
    phys_const::MassHelium3 / 2.f,
    phys_const::MassAlpha / 2.f,
    std::numeric_limits<float>::infinity(), // Pi0, Z=0
    std::numeric_limits<float>::infinity(), // Photon, Z=0
    std::numeric_limits<float>::infinity(), // K0, Z=0
    std::numeric_limits<float>::infinity(), // Lambda, Z=0
    phys_const::MassHyperTriton,
    phys_const::MassHyperhydrog4,
    phys_const::MassXiMinus,
    phys_const::MassOmegaMinus,
    phys_const::MassHyperHelium4 / 2.f,
    phys_const::MassHyperHelium5 / 2.f};
  static constexpr float sCharges[NIDs] = {1, 1, 1, 1, 1, 1, 1, 2, 2, 0, 0, 0, 0, 1, 1, 1, 1, 2, 2};

  PID() = default;
  PID(ID id) : mID(id) {}
  PID(const char* name);
  PID(const PID& src) = default;
  PID& operator=(const PID& src) = default;

  static PID PDG2PID(int pdg);
  ID getID() const { return mID; }
  operator ID() const { return getID(); }
  void setFromPDG(int pdg) { *this = PDG2PID(pdg); }
  float getMass() const { return getMass(mID); }
  float getMass2() const { return getMass2(mID); }
  float getMassInv() const { return getMassInv(mID); }
  float getMass2Z() const { return getMass2Z(mID); }
  const char* getName() const { return getName(mID); }
  int getCharge() const { return getCharge(mID); }

  static float getMass(ID id) { return sMasses[id]; }
  static float getMass2(ID id) { return sMasses2[id]; }
  static float getMassInv(ID id) { return sMassesI[id]; }
  static float getMass2Z(ID id) { return sMasses2Z[id]; }
  static int getCharge(ID id) { return sCharges[id]; }
  static const char* getName(ID id) { return sNames[id]; }

 private:
  ID mID = Pion;

  // are 2 strings equal ? (trick from Giulio)
  static constexpr bool sameStr(char const* x, char const* y) { return !*x && !*y ? true : /* default */ (*x == *y && sameStr(x + 1, y + 1)); }
  static constexpr ID nameToID(char const* name, ID id) { return id > Last ? id : sameStr(name, sNames[id]) ? id
                                                                                                            : nameToID(name, id + 1); }

  ClassDefNV(PID, 1);
};

#endif
