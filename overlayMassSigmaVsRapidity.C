#if !defined(__CINT__) || defined(__MAKECINT__)
#include <TFile.h>
#include <TH1D.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TString.h>
#include <iostream>
#include <vector>
#endif

void overlayMassSigmaVsRapidity(
  const std::vector<TString>& files,
  const std::vector<TString>& labels,
  const char* histoName = "hMassSigmaVsRapidity",
  const char* outFile = "overlayD0MassSigmaVsRapidity.png",
  const char* title = "D^{0} mass #sigma from Gaussian fit vs rapidity;Rapidity;#sigma_{M} (GeV/c^{2})")
{
  TH1::AddDirectory(false);
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(0);

  if (files.empty()) {
    std::cerr << "No files provided!" << std::endl;
    return;
  }
  if (files.size() != labels.size()) {
    std::cerr << "Files and labels must have the same size!" << std::endl;
    return;
  }

  int colors[] = {kRed+1, kBlue+1, kGreen+2, kMagenta+1,
                  kOrange+7, kCyan+2, kViolet+1, kPink+2};
  int markers[] = {kFullCircle, kFullSquare, kFullTriangleUp, kFullTriangleDown,
                   kFullDiamond, kFullStar, kFullCross, kOpenCircle};

  // Load all histograms first
  std::vector<TH1D*> histos;
  std::vector<TString> validLabels;
  for (int i = 0; i < (int)files.size(); ++i) {
    TFile* f = TFile::Open(files[i], "READ");
    if (!f || f->IsZombie()) {
      std::cerr << "Cannot open " << files[i] << std::endl;
      continue;
    }
    TH1D* h = dynamic_cast<TH1D*>(f->Get(histoName));
    if (!h) {
      std::cerr << "Histogram " << histoName << " not found in " << files[i] << std::endl;
      f->Close();
      continue;
    }
    h->SetDirectory(nullptr);
    f->Close();

    h->SetName(Form("h_%d", (int)histos.size()));
    h->SetStats(0);
    h->SetLineColor(colors[histos.size() % 8]);
    h->SetMarkerColor(colors[histos.size() % 8]);
    h->SetMarkerStyle(markers[histos.size() % 8]);
    h->SetMarkerSize(1.2);
    h->SetLineWidth(2);
    histos.push_back(h);
    validLabels.push_back(labels[i]);
  }

  if (histos.empty()) {
    std::cerr << "No histograms loaded!" << std::endl;
    return;
  }

  const bool hasRatio = histos.size() > 1;
  const float topFrac   = hasRatio ? 0.65f : 1.0f;
  const float botFrac   = 1.0f - topFrac;

  TCanvas* c = new TCanvas("c", "Mass #sigma vs Rapidity", 900, hasRatio ? 900 : 700);
  c->SetFillColor(0);

  // --- Top pad ---
  TPad* padTop = new TPad("padTop", "", 0., botFrac, 1., 1.);
  padTop->SetLeftMargin(0.14);
  padTop->SetRightMargin(0.05);
  padTop->SetTopMargin(0.08);
  padTop->SetBottomMargin(0.02);
  padTop->SetGrid();
  padTop->Draw();
  padTop->cd();

  TLegend* leg = new TLegend(0.16, 0.55, 0.62, 0.92);
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);
  leg->SetTextSize(0.05);

  for (int i = 0; i < (int)histos.size(); ++i) {
    TH1D* h = histos[i];
    if (i == 0) {
      h->SetTitle(title);
      h->GetXaxis()->SetLabelSize(0.);
      h->GetXaxis()->SetTitleSize(0.);
      h->GetYaxis()->SetTitle("#sigma_{M} (GeV/c^{2})");
      h->GetYaxis()->SetTitleOffset(1.2);
      h->GetYaxis()->SetTitleSize(0.06);
      h->GetYaxis()->SetLabelSize(0.05);
      h->Draw("E P");
    } else {
      h->Draw("E P SAME");
    }
    leg->AddEntry(h, validLabels[i], "lp");
  }
  leg->Draw();

  // --- Bottom pad (ratios) ---
  if (hasRatio) {
    c->cd();
    TPad* padBot = new TPad("padBot", "", 0., 0., 1., botFrac);
    padBot->SetLeftMargin(0.14);
    padBot->SetRightMargin(0.05);
    padBot->SetTopMargin(0.02);
    padBot->SetBottomMargin(0.28);
    padBot->SetGrid();
    padBot->Draw();
    padBot->cd();

    bool firstRatio = true;
    for (int i = 1; i < (int)histos.size(); ++i) {
      TH1D* hRatio = (TH1D*)histos[i]->Clone(Form("hRatio_%d", i));
      hRatio->Reset();      // keep binning, clear contents

      for (int b = 1; b <= histos[i]->GetNbinsX(); ++b) {
          double A  = histos[i]->GetBinContent(b);
          double B  = histos[0]->GetBinContent(b);
          double eA = histos[i]->GetBinError(b);
          double eB = histos[0]->GetBinError(b);

          if (B != 0.0) {
              double R = A / B;
              double eR = 0.0;

              if (A != 0.0) {
                  eR = R * std::sqrt((eA/A)*(eA/A) + (eB/B)*(eB/B));
              } else {
                  // safer when A=0 to avoid division by zero in relative form
                  eR = eA / B;
              }

              hRatio->SetBinContent(b, R);
              hRatio->SetBinError(b, eR);
          } else {
              hRatio->SetBinContent(b, 0.0);
              hRatio->SetBinError(b, 0.0);
          }
      }
      hRatio->SetStats(0);

      if (firstRatio) {
        hRatio->SetTitle("");
        hRatio->GetXaxis()->SetTitle("Rapidity");
        hRatio->GetXaxis()->SetTitleSize(0.12);
        hRatio->GetXaxis()->SetLabelSize(0.10);
        hRatio->GetYaxis()->SetTitle("Ratio to " + validLabels[0]);
        hRatio->GetYaxis()->SetTitleSize(0.10);
        hRatio->GetYaxis()->SetTitleOffset(0.55);
        hRatio->GetYaxis()->SetLabelSize(0.09);
        hRatio->GetYaxis()->SetNdivisions(505);
        hRatio->SetMinimum(0.5);
        hRatio->SetMaximum(2);
        hRatio->Draw("E P");
        firstRatio = false;

      } else {
        hRatio->Draw("E P SAME");
      }
    }

  }

  c->SaveAs(outFile);
  std::cout << "Saved: " << outFile << std::endl;
}

// Convenience entry point: pass file paths and labels directly from the command line or interactively
void overlayMassSigmaVsRapidity()
{

  std::vector<TString> labels = {
    "No plate",
    "Tube only",
    "Aluminium",
    "CarbonFoam #rho = 0.6 g/cm^{3}"
  };

  // Edit these to match your files and desired labels
  std::vector<TString> files = {
    "D0_No_Plate/D0_analysis.root",
    "D0_Air/D0_analysis.root",
    "D0_Aluminium/D0_analysis.root",
    "D0_CarbonFoam/D0_analysis.root"
  };

  overlayMassSigmaVsRapidity(files, labels);

  // Edit these to match your files and desired labels
  std::vector<TString> filesLc = {
    "LAMBDA_C_No_Plate/LambdaC_analysis.root",
    "LAMBDA_C_Air/LambdaC_analysis.root",
    "LAMBDA_C_Aluminium/LambdaC_analysis.root",
    "LAMBDA_C_CarbonFoam/LambdaC_analysis.root",
  };

  overlayMassSigmaVsRapidity(filesLc, labels, "hLcMassSigmaVsY", "overlayLambdaCMassSigmaVsRapidity.png", "#Lambda_{c} mass #sigma from Gaussian fit vs rapidity;Rapidity;#sigma_{M} (GeV/c^{2})");

  // Edit these to match your files and desired labels
  std::vector<TString> filesJpsi = {
    "JPSI_No_Plate/dimuon_analysisVT.root",
    "JPSI_Air/dimuon_analysisVT.root",
    "JPSI_Aluminium/dimuon_analysisVT.root",
    "JPSI_CarbonFoam/dimuon_analysisVT.root",
  };

  overlayMassSigmaVsRapidity(filesJpsi, labels, "hMassSigmaVsRapidity", "overlayJpsiCMassSigmaVsRapidity.png", "J/#psi mass #sigma from Gaussian fit vs rapidity;Rapidity;#sigma_{M} (GeV/c^{2})");

  // Edit these to match your files and desired labels
  std::vector<TString> filesOMEGA = {
    "OMEGA_No_Plate/dimuon_analysisVT.root",
    "OMEGA_Air/dimuon_analysisVT.root",
    "OMEGA_Aluminium/dimuon_analysisVT.root",
    "OMEGA_CarbonFoam/dimuon_analysisVT.root",
  };

  overlayMassSigmaVsRapidity(filesOMEGA, labels, "hMassSigmaVsRapidity", "overlayOmegaMassSigmaVsRapidity.png", "#omega mass #sigma from Gaussian fit vs rapidity;Rapidity;#sigma_{M} (GeV/c^{2})");

}
