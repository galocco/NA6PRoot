#include <TCanvas.h>
#include <TFile.h>
#include <TTree.h>
#include <TH1D.h>
#include <TStyle.h>
#include <TLatex.h>
#include <TSystem.h>
#include <iostream>
#include <vector>

#include "NA6PTrack.h"

// Plot the fraction of reconstructed muon tracks in each mStatusMS category.
// Usage from ROOT:
//   .L plotMSStatus.C+
//   plotMSStatus("TracksMuonSpec.root", "mStatusMS.png");
void plotMSStatus(const char* fileName = "TracksMuonSpec.root",
                  const char* outputName = "mStatusMS.png")
{
  auto file = TFile::Open(fileName, "READ");
  if (!file || file->IsZombie()) {
    std::cerr << "Cannot open " << fileName << std::endl;
    return;
  }

  auto tree = dynamic_cast<TTree*>(file->Get("tracksMuonSpec"));
  if (!tree) {
    std::cerr << "Tree tracksMuonSpec not found in " << fileName << std::endl;
    return;
  }

  std::vector<NA6PTrack>* tracks = nullptr;
  tree->SetBranchAddress("MuonSpec", &tracks);

  constexpr int nStatus = 6;
  const char* labels[nStatus] = {
    "Not MS", "Full MS+MID", "MS not matched", "Matched, not refitted 1",
    "Matched, not refitted 2", "Matched and refitted"
  };

  TH1D h("hStatusMS", "Muon-spectrometer track status;Status;Fraction (%)",
         nStatus, -0.5, nStatus - 0.5);
  for (int i = 0; i < nStatus; ++i)
    h.GetXaxis()->SetBinLabel(i + 1, labels[i]);

  Long64_t nTracks = 0;
  for (Long64_t i = 0; i < tree->GetEntries(); ++i) {
    tree->GetEntry(i);
    if (!tracks)
      continue;
    for (const auto& track : *tracks) {
      const int status = track.getStatusMS();
      if (status >= 0 && status < nStatus) {
        h.Fill(status);
        ++nTracks;
      }
    }
  }

  if (nTracks > 0)
    h.Scale(100. / nTracks);

  gStyle->SetOptStat(0);
  auto canvas = new TCanvas("cStatusMS", "mStatusMS", 1100, 650);
  canvas->SetBottomMargin(0.25);
  h.SetFillColor(kAzure - 4);
  h.SetBarWidth(0.8);
  h.SetBarOffset(0.1);
  h.SetMinimum(0.);
  h.Draw("BAR");

  TLatex text;
  text.SetTextAlign(22);
  text.SetTextSize(0.03);
  for (int i = 1; i <= nStatus; ++i) {
    if (h.GetBinContent(i) > 0)
      text.DrawLatex(h.GetBinCenter(i), h.GetBinContent(i) + 0.02,
                     Form("%.2f%%", h.GetBinContent(i)));
  }

  canvas->SaveAs(outputName);
  std::cout << "Processed " << nTracks << " tracks" << std::endl;
  std::cout << "Saved " << outputName << std::endl;
}
