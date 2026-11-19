// Run from NA6PRoot after sourcing install/init.sh:
// root -l -b -q 'plotMSWallEnergyLoss.C("HitsMuonSpecModular.root")'
// The default z values are the planes immediately before and after the wall
// in na6pLayout_Jelle.ini. Override them for another geometry.

#include <TCanvas.h>
#include <TFile.h>
#include <TH1D.h>
#include <TLatex.h>
#include <TLine.h>
#include <TStyle.h>
#include <TTree.h>
#include "sim/include/NA6PMuonSpecModularHit.h"

#include <cmath>
#include <iostream>
#include <map>
#include <memory>
#include <vector>

void plotMSWallEnergyLoss(const char* input = "HitsMuonSpecModular.root",
                          const char* output = "ms3_ms4_energy_loss.png",
                          double zMS3 = 593., double zMS4 = 789.5,
                          double zTolerance = 2.)
{
  constexpr double muonMassGeV = 0.1056583755;
  auto file = std::unique_ptr<TFile>(TFile::Open(input, "READ"));
  if (!file || file->IsZombie()) {
    std::cerr << "Cannot open " << input << '\n';
    return;
  }
  auto* tree = file->Get<TTree>("hitsMuonSpecModular");
  if (!tree || !tree->GetBranch("MuonSpecModular")) {
    std::cerr << "Missing hitsMuonSpecModular/MuonSpecModular in " << input << '\n';
    return;
  }

  std::vector<NA6PMuonSpecModularHit>* hits = nullptr;
  tree->SetBranchAddress("MuonSpecModular", &hits);
  auto hist = std::make_unique<TH1D>("hMSWallEnergyLoss",
      "Muon energy loss from MS3 to MS4 (wall included);E_{MS3,out} - E_{MS4,in} (GeV);Tracks",
      150, -1., 14.);
  Long64_t pairs = 0;
  for (Long64_t event = 0; event < tree->GetEntries(); ++event) {
    tree->GetEntry(event);
    if (!hits) continue;
    // One hit on each plane per track; retain the closest hit if a track
    // leaves multiple hits near the same plane.
    struct Pair { const NA6PMuonSpecModularHit* before = nullptr;
                  const NA6PMuonSpecModularHit* after = nullptr; };
    std::map<int, Pair> byTrack;
    for (const auto& hit : *hits) {
      const double z = hit.getZ();
      auto& pair = byTrack[hit.getTrackID()];
      if (std::abs(z - zMS3) <= zTolerance &&
          (!pair.before || std::abs(z - zMS3) < std::abs(pair.before->getZ() - zMS3)))
        pair.before = &hit;
      if (std::abs(z - zMS4) <= zTolerance &&
          (!pair.after || std::abs(z - zMS4) < std::abs(pair.after->getZ() - zMS4)))
        pair.after = &hit;
    }
    for (const auto& [trackID, pair] : byTrack) {
      if (!pair.before || !pair.after) continue;
      const auto* a = pair.before;
      const auto* b = pair.after;
      const double e3 = std::sqrt(a->getPXOut()*a->getPXOut() +
                                  a->getPYOut()*a->getPYOut() +
                                  a->getPZOut()*a->getPZOut() + muonMassGeV*muonMassGeV);
      const double e4 = std::sqrt(b->getPXIn()*b->getPXIn() +
                                  b->getPYIn()*b->getPYIn() +
                                  b->getPZIn()*b->getPZIn() + muonMassGeV*muonMassGeV);
      hist->Fill(e3 - e4);
      ++pairs;
    }
  }
  if (!pairs) {
    std::cerr << "No matched MS3/MS4 hits; check plane z positions and tolerance.\n";
    return;
  }
  // The peak-bin center is a reproducible, binning-dependent MPV estimate.
  // Keep the mean over all matched pairs, including entries outside the view.
  const double mpv = hist->GetBinCenter(hist->GetMaximumBin());
  const double mean = hist->GetMean();
  const double yMax = 1.18 * hist->GetMaximum();

  gStyle->SetOptStat(0);
  TCanvas canvas("cMSWallEnergyLoss", "MS3 to MS4 energy loss", 1100, 750);
  canvas.SetLeftMargin(0.12);
  canvas.SetRightMargin(0.05);
  canvas.SetTopMargin(0.10);
  canvas.SetBottomMargin(0.13);
  hist->SetStats(false);
  hist->SetTitle("Energy loss in the absorber");
  hist->GetXaxis()->SetTitle("E_{MS3,out} - E_{MS4,in} (GeV)");
  hist->GetYaxis()->SetTitle("Counts");
  hist->GetXaxis()->SetRangeUser(0., 4.);
  hist->GetXaxis()->SetTitleSize(0.05);
  hist->GetYaxis()->SetTitleSize(0.05);
  hist->GetXaxis()->SetLabelSize(0.042);
  hist->GetYaxis()->SetLabelSize(0.042);
  hist->GetXaxis()->SetTitleOffset(1.08);
  hist->GetYaxis()->SetTitleOffset(1.12);
  hist->SetLineColor(kAzure + 2);
  hist->SetFillColorAlpha(kAzure - 9, 0.45);
  hist->SetLineWidth(2);
  hist->SetMinimum(0.);
  hist->SetMaximum(yMax);
  hist->Draw("HIST");

  TLine mpvLine(mpv, 0., mpv, 0.92 * yMax);
  mpvLine.SetLineColor(kOrange + 7);
  mpvLine.SetLineStyle(2);
  mpvLine.SetLineWidth(2);
  mpvLine.Draw();
  TLine meanLine(mean, 0., mean, 0.92 * yMax);
  meanLine.SetLineColor(kRed + 1);
  meanLine.SetLineStyle(7);
  meanLine.SetLineWidth(2);
  meanLine.Draw();

  TLatex label;
  label.SetNDC();
  label.SetTextSize(0.043);
  label.SetTextColor(kBlack);
  label.DrawLatex(0.58, 0.80, Form("MPV (peak bin): %.3f GeV", mpv));
  label.DrawLatex(0.58, 0.73, Form("Mean: %.3f GeV", mean));

  canvas.SaveAs(output);
  std::cout << "Plotted " << pairs << " track pairs to " << output << '\n';
}
