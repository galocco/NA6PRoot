#include <TFile.h>
#include <TTree.h>
#include <TClonesArray.h>
#include <TH1D.h>
#include <TCanvas.h>
#include <TLorentzVector.h>
#include <iostream>
#include "NA6PTrack.h"

// ⚠️ includi l'header della classe Matching


void InvariantMassPlot(const char* filename = "TracksMuonSpecMatching_MC.root",
                   const char* outputname = "jpsiMuonSpec.png",
                   const char* treename = "tracksMuonSpec",
                   const char* branchname = "MuonSpec")
{
  // =========================
  // Apri il file
  // =========================
  TFile* f = TFile::Open(filename);
  if (!f || f->IsZombie()) {
    std::cerr << "Errore apertura file ROOT" << std::endl;
    return;
  }

  // =========================
  // Prendi il TTree
  // =========================
  TTree* tree = (TTree*)f->Get(treename);
  if (!tree) {
    std::cerr << "TTree " << treename << " non trovato" << std::endl;
    return;
  }

  // =========================
  // Branch Matching
  // =========================
  std::vector<NA6PTrack>* tracks = nullptr;
  tree->SetBranchAddress(branchname, &tracks);



  const double jpsiMass = 3.0969; // GeV/c^2
  const double muMass = 0.105658; // GeV/c^2
  // =========================
  // Istogramma massa invariante
  // =========================
  TH1D* hMinv = new TH1D(
    "hMinv",
    "Invariant mass of muon pairs;M_{#mu#mu} [GeV/c^{2}];Counts",
    100, jpsiMass - 0.5, jpsiMass + 0.5
  );
  // =========================
  // Loop sugli eventi
  // =========================
  Long64_t nEvents = tree->GetEntries();
  for (Long64_t iEv = 0; iEv < nEvents; ++iEv) {
    tree->GetEntry(iEv);

    int nTracks = tracks->size();
    if (nTracks < 2) continue;

    // =========================
    // Loop sulle coppie
    // =========================
    for (int i = 0; i < nTracks; ++i) {
      auto& t1 = tracks->at(i);
      if(t1.getParticleID() < 0) continue; // solo particelle identificate come muoni
      for (int j = i + 1; j < nTracks; ++j) {
        auto& t2 = tracks->at(j);
        if(t2.getParticleID() < 0) continue; // solo particelle identificate come muoni

        // opzionale: solo coppie opposta carica
        if (t1.getCharge() * t2.getCharge() >= 0) continue;

        TLorentzVector p1, p2;
        double pxyz[3];
        t1.getPXYZ(pxyz);
        p1.SetXYZM(
          pxyz[0],
          pxyz[1],
          pxyz[2],
          muMass
        );

        t2.getPXYZ(pxyz);
        p2.SetXYZM(
          pxyz[0],
          pxyz[1],
          pxyz[2],
          muMass
        );

        double minv = (p1 + p2).M();
        hMinv->Fill(minv);
      }
    }
  }

  // =========================
  // Plot
  // =========================
  TCanvas* c = new TCanvas("c", "Invariant Mass", 800, 600);
  hMinv->Draw();

  c->SaveAs(outputname);

  std::cout << "Fatto. Coppie analizzate e istogramma prodotto." << std::endl;
}

void InvariantMass()
{
  InvariantMassPlot("TracksMuonSpecMatching_MC.root", "jpsiMuonSpec.png", "tracksMatching", "Matching");
  InvariantMassPlot("TracksVerTelMatching_MC.root", "jpsiVerTel.png", "tracksMatching", "Matching");
  InvariantMassPlot("TracksMatching_MC.root", "jpsiMatch.png", "tracksMatching", "Matching");
}