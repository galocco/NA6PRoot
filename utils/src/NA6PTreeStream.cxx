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

#include "NA6PTreeStream.h"
#include <TBranch.h>

//_________________________________________________
NA6PTreeStream::NA6PTreeStream(const char* treename) : mTree(treename, treename)
{
  //
  // Standard ctor
}

//_________________________________________________
int NA6PTreeStream::CheckIn(Char_t type, const void* pointer)
{
  // Insert object

  if (mCurrentIndex >= static_cast<int>(mElements.size())) {
    auto& element = mElements.emplace_back();
    element.type = type;
    TString name = mNextName;
    if (name.Length()) {
      if (mNextNameCounter > 0) {
        name += mNextNameCounter;
      }
    } else {
      name = TString::Format("B%d.", static_cast<int>(mElements.size()));
    }
    element.name = name.Data();
    element.ptr = pointer;
    element.arsize = mNextArraySize;
    mNextArraySize = 1; // reset
  } else {
    auto& element = mElements[mCurrentIndex];
    if (element.type != type) {
      mStatus++;
      return 1; // mismatched data element
    }
    element.ptr = pointer;
  }
  mCurrentIndex++;
  return 0;
}

//_________________________________________________
void NA6PTreeStream::BuildTree()
{
  // Build the Tree

  int entriesFilled = mTree.GetEntries();
  if (mBranches.size() < mElements.size()) {
    mBranches.resize(mElements.size());
  }

  TString name;
  TBranch* br = nullptr;
  for (int i = 0; i < static_cast<int>(mElements.size()); i++) {
    //
    auto& element = mElements[i];
    if (mBranches[i]) {
      continue;
    }
    name = element.name.data();
    if (name.IsNull()) {
      name = TString::Format("B%d", i);
    }
    if (element.cls) {
      br = mTree.Branch(name.Data(), element.cls->GetName(), const_cast<void**>(&element.ptr));
      mBranches[i] = br;
      if (entriesFilled) {
        br->SetAddress(nullptr);
        for (int ientry = 0; ientry < entriesFilled; ientry++) {
          br->Fill();
        }
        br->SetAddress(const_cast<void**>(&element.ptr));
      }
    }

    if (element.type > 0) {
      TString nameC;
      if (element.arsize > 1) {
        nameC = TString::Format("%s[%d]/%c", name.Data(), element.arsize,
                                element.type);
      } else {
        nameC = TString::Format("%s/%c", name.Data(), element.type);
      }
      br = mTree.Branch(name.Data(), const_cast<void*>(element.ptr), nameC.Data());
      if (entriesFilled) {
        br->SetAddress(nullptr);
        for (int ientry = 0; ientry < entriesFilled; ientry++) {
          br->Fill();
        }
        br->SetAddress(const_cast<void*>(element.ptr));
      }
      mBranches[i] = br;
    }
  }
}

//_________________________________________________
void NA6PTreeStream::Fill()
{
  // Fill the tree

  int entries = mElements.size();
  if (entries > mTree.GetNbranches()) {
    BuildTree();
  }
  for (int i = 0; i < entries; i++) {
    auto& element = mElements[i];
    if (!element.type) {
      continue;
    }
    auto br = mBranches[i];
    if (br) {
      if (element.type) {
        br->SetAddress(const_cast<void*>(element.ptr));
      }
    }
  }
  if (!mStatus) {
    mTree.Fill(); // fill only in case of non conflicts
  }
  mStatus = 0;
}

//_________________________________________________
NA6PTreeStream& NA6PTreeStream::Endl()
{
  // Perform pseudo endl operation

  if (mTree.GetNbranches() == 0) {
    BuildTree();
  }
  Fill();
  mStatus = 0;
  mCurrentIndex = 0;
  return *this;
}

//_________________________________________________
NA6PTreeStream& NA6PTreeStream::operator<<(const Char_t* name)
{
  // Stream the branch name
  if (name[0] == '\n') {
    return Endl();
  }

  // if tree was already defined ignore
  if (mTree.GetEntries() > 0) {
    return *this;
  }

  int arsize = 1;

  // check branch name if tree was not
  Int_t last = 0;
  for (last = 0;; last++) {
    if (name[last] == 0) {
      break;
    }
  }
  if (last > 0 && name[last - 1] == '=') {
    mNextName = name;
    mNextName[last - 1] = 0; // remove '=' from string
    mNextNameCounter = 0;

    TString inName{name};
    auto brkStaPos = inName.Index('[');

    if (brkStaPos != kNPOS) {
      auto brkEndPos = inName.Index(']');
      if (brkEndPos != kNPOS && brkEndPos > brkStaPos + 1) {
        TString size = inName(brkStaPos + 1, brkEndPos - brkStaPos - 1);
        arsize = size.Atoi();
        mNextName = inName(0, brkStaPos); // use parsed name
      }
    }
  }

  mNextArraySize = arsize;

  return *this;
}
