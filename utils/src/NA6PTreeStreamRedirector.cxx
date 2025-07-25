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

#include "NA6PTreeStreamRedirector.h"
#include <TFile.h>
#include <TLeaf.h>
#include <cstring>

//_________________________________________________
NA6PTreeStreamRedirector::NA6PTreeStreamRedirector(const char* fname, const char* option)
{
  // Constructor

  TString name(fname);
  if (!name.IsNull()) {
    mOwnDirectory = std::unique_ptr<TDirectory>(TFile::Open(fname, option));
    mDirectory = mOwnDirectory.get();
  } else {
    mDirectory = gDirectory;
  }
}

//_________________________________________________
NA6PTreeStreamRedirector::~NA6PTreeStreamRedirector()
{
  // Destructor
  Close(); // write the tree to the selected file
}

//_________________________________________________
void NA6PTreeStreamRedirector::SetFile(TFile* sfile)
{
  // set the external file
  SetDirectory(sfile);
}

//_________________________________________________
void NA6PTreeStreamRedirector::SetDirectory(TDirectory* sfile)
{
  // Set the external directory
  // In case other directory already attached old file is closed before
  // Redirector will be the owner of file ?

  if (mOwnDirectory) {
    mDirectory->Close();
    mOwnDirectory.reset();
  }
  mDirectory = sfile;
}

//_____________________________________________________
NA6PTreeStream& NA6PTreeStreamRedirector::operator<<(Int_t id)
{
  // return reference to the data layout with given identifier
  // if not existing - creates new

  for (auto& layout : mDataLayouts) {
    if (layout->getID() == id) {
      return *layout.get();
    }
  }

  TDirectory* backup = gDirectory;
  mDirectory->cd();
  mDataLayouts.emplace_back(std::unique_ptr<NA6PTreeStream>(new NA6PTreeStream(Form("Tree%d", id))));
  auto layout = mDataLayouts.back().get();
  layout->setID(id);
  if (backup) {
    backup->cd();
  }
  return *layout;
}

//_________________________________________________
NA6PTreeStream& NA6PTreeStreamRedirector::operator<<(const char* name)
{
  // return reference to the data layout with given identifier
  // if not existing - creates new

  for (auto& layout : mDataLayouts) {
    if (!std::strcmp(layout->getName(), name)) {
      return *layout.get();
    }
  }

  // create new
  TDirectory* backup = gDirectory;
  mDirectory->cd();
  mDataLayouts.emplace_back(std::unique_ptr<NA6PTreeStream>(new NA6PTreeStream(name)));
  auto layout = mDataLayouts.back().get();
  layout->setID(-1);
  if (backup) {
    backup->cd();
  }
  return *layout;
}

//_________________________________________________
void NA6PTreeStreamRedirector::Close()
{
  // flush and close
  if (!mDirectory) {
    return;
  }
  TDirectory* backup = gDirectory;
  mDirectory->cd();
  for (auto& layout : mDataLayouts) {
    layout->getTree().Write(layout->getName(), TObject::kOverwrite);
  }
  mDataLayouts.clear();
  if (backup) {
    backup->cd();
  }

  if (mOwnDirectory) {
    mDirectory->Close();
  } else {
    mDirectory = nullptr;
  }
}

//_________________________________________________
void NA6PTreeStreamRedirector::FixLeafNameBug(TTree* tree)
{
  // On the fly BUG FIX for name and titles of branches and Leave:
  //     renaming Leaves and changing title of branches to be the same as Branch name
  // Explanation of FIX:
  //    In the  friend tree Join logic it is assumed leaf names are used to find appropraiat primary/secondary keys
  //    For the standard queries however the branch names are used to identify data
  //    Hovewer in the Branch constructor it is not checked
  // As a consequence  - in case the name of the leave and and the name of branch is not the same  + freind trees are
  // sparse
  //    wrong joins ( unrelated pair of information) are used
  // FIX:
  //   To be able to use friend trees with proper indexing (in case of sarse trees) branches and leaves has to be named
  //   consistently
  //   In this routine bnrach name is taken as a reference and branch title and leave name titles are renamed
  //   After the fix unit test code with pairs of sprse friend trees worked properly
  // Side effects of fix:
  //
  if (!tree) {
    return;
  }
  TObjArray* brArray = tree->GetListOfBranches();
  TObjArray* lArray = tree->GetListOfLeaves();
  for (int i = 0; i < brArray->GetLast(); i++) {
    TBranch* br = (TBranch*)brArray->At(i);
    TString brTitle(br->GetTitle());
    if (!brTitle.Contains(br->GetName())) {
      int pos = brTitle.First("/");
      TString leafName;
      if (pos < brTitle.Length()) {
        brTitle[pos] = 0;
        leafName = TString::Format("%s", brTitle.Data()).Data();
        TLeaf* leaf = (TLeaf*)lArray->FindObject(leafName);
        if (leaf) {
          leaf->SetName(br->GetName());
          leaf->SetTitle(br->GetName());
          br->SetTitle(TString::Format("%s/%s", br->GetName(), &(brTitle.Data()[pos + 1])).Data());
        }
      }
    }
  }
}
