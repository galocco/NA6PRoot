#include <TCanvas.h>
#include <TFile.h>
#include <TF1.h>
#include <TH1.h>
#include <TLegend.h>
#include <TLine.h>
#include <TPad.h>
#include <TStyle.h>
#include <TString.h>

#include <algorithm>
#include <cmath>
#include <iostream>

namespace {

TH1* loadHistogram(TFile& input, const char* fileName, const char* cloneName)
{
  auto* source = dynamic_cast<TH1*>(input.Get("hdimumassConstrained"));
  if (!source) {
    std::cerr << "Missing histogram 'hdimumassConstrained' in " << fileName << std::endl;
    return nullptr;
  }

  auto* histogram = dynamic_cast<TH1*>(source->Clone(cloneName));
  histogram->SetDirectory(nullptr);
  histogram->Sumw2();
  return histogram;
}

bool sameBinning(const TH1& first, const TH1& second)
{
  return first.GetNbinsX() == second.GetNbinsX() &&
         first.GetXaxis()->GetXmin() == second.GetXaxis()->GetXmin() &&
         first.GetXaxis()->GetXmax() == second.GetXaxis()->GetXmax();
}

TString sigmaText(const TH1& histogram)
{
  const auto* fit = histogram.GetFunction("gaus");
  if (fit) {
    return Form("#sigma = %.3g #pm %.2g GeV/c^{2}", std::abs(fit->GetParameter(2)), fit->GetParError(2));
  }
  return Form("RMS = %.3g GeV/c^{2}", histogram.GetRMS());
}

} // namespace

void compareRobertaDimuonMass(
    const char* robertaFile = "debug_daiki/dimuon_analysisMS.root",
    const char* robertaFullFile = "debug_roberta_full/dimuon_analysisMS.root",
    const char* outputFile = "comparison_hdimumassConstrained_roberta_vs_roberta_full.png",
    bool normalize = true)
{
  TFile robertaInput(robertaFile, "READ");
  TFile robertaFullInput(robertaFullFile, "READ");
  if (robertaInput.IsZombie() || robertaFullInput.IsZombie()) {
    std::cerr << "Cannot open one or both input files." << std::endl;
    return;
  }

  auto* roberta = loadHistogram(robertaInput, robertaFile, "hMassRoberta");
  auto* robertaFull = loadHistogram(robertaFullInput, robertaFullFile, "hMassRobertaFull");
  if (!roberta || !robertaFull) {
    return;
  }
  if (!sameBinning(*roberta, *robertaFull)) {
    std::cerr << "The two hdimumassConstrained histograms have incompatible binning." << std::endl;
    return;
  }

  const double robertaIntegral = roberta->Integral();
  const double robertaFullIntegral = robertaFull->Integral();
  std::cout << "Roberta integral: " << robertaIntegral
            << ", Roberta full integral: " << robertaFullIntegral << std::endl;
  const bool normalized = normalize && robertaIntegral > 0. && robertaFullIntegral > 0.;
  if (normalize && !normalized) {
    std::cerr << "At least one histogram is empty: drawing unnormalized counts." << std::endl;
  }
  if (normalized) {
    roberta->Scale(1. / robertaIntegral);
    robertaFull->Scale(1. / robertaFullIntegral);
  }

  gStyle->SetOptStat(0);
  auto* canvas = new TCanvas("cRobertaDimuonMass", "Constrained dimuon mass comparison", 1000, 900);
  auto* topPad = new TPad("topPad", "", 0., 0.30, 1., 1.);
  auto* ratioPad = new TPad("ratioPad", "", 0., 0., 1., 0.30);
  topPad->SetBottomMargin(0.02);
  ratioPad->SetTopMargin(0.03);
  ratioPad->SetBottomMargin(0.32);
  topPad->Draw();
  ratioPad->Draw();

  topPad->cd();
  roberta->SetLineColor(kAzure + 2);
  roberta->SetLineWidth(3);
  robertaFull->SetLineColor(kOrange + 7);
  robertaFull->SetLineWidth(3);
  robertaFull->SetLineStyle(2);
  roberta->SetTitle("Constrained dimuon mass");
  roberta->GetYaxis()->SetTitle(normalized ? "Normalized counts" : "Counts");
  roberta->GetYaxis()->SetTitleSize(0.05);
  roberta->GetYaxis()->SetTitleOffset(1.05);
  roberta->GetXaxis()->SetLabelSize(0.);
  const double topMaximum = std::max(roberta->GetMaximum(), robertaFull->GetMaximum());
  roberta->SetMaximum(topMaximum > 0. ? 1.20 * topMaximum : 1.);
  roberta->Draw("hist");
  robertaFull->Draw("hist same");

  auto* legend = new TLegend(0.47, 0.68, 0.88, 0.88);
  legend->SetBorderSize(0);
  legend->SetFillStyle(0);
  legend->AddEntry(roberta, Form("Roberta (N = %.0f, %s)", robertaIntegral, sigmaText(*roberta).Data()), "l");
  legend->AddEntry(robertaFull, Form("Roberta full (N = %.0f, %s)", robertaFullIntegral, sigmaText(*robertaFull).Data()), "l");
  legend->Draw();

  ratioPad->cd();
  auto* ratio = dynamic_cast<TH1*>(robertaFull->Clone("hMassRobertaFullOverRoberta"));
  ratio->SetDirectory(nullptr);
  ratio->Divide(roberta);
  ratio->SetTitle("");
  ratio->SetLineColor(kBlack);
  ratio->SetLineWidth(2);
  ratio->GetYaxis()->SetTitle("Full / Roberta");
  ratio->GetYaxis()->SetNdivisions(505);
  ratio->GetYaxis()->SetTitleSize(0.10);
  ratio->GetYaxis()->SetTitleOffset(0.48);
  ratio->GetYaxis()->SetLabelSize(0.09);
  ratio->GetXaxis()->SetTitle("m_{#mu#mu} (GeV/c^{2})");
  ratio->GetXaxis()->SetTitleSize(0.12);
  ratio->GetXaxis()->SetTitleOffset(1.05);
  ratio->GetXaxis()->SetLabelSize(0.10);
  ratio->SetMinimum(0.);
  ratio->SetMaximum(2.);
  ratio->Draw("hist");

  const double xMin = ratio->GetXaxis()->GetXmin();
  const double xMax = ratio->GetXaxis()->GetXmax();
  auto* unity = new TLine(xMin, 1., xMax, 1.);
  unity->SetLineStyle(2);
  unity->Draw();

  canvas->SaveAs(outputFile);
  std::cout << "Saved comparison to " << outputFile << std::endl;
}
