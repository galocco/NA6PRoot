#if !defined(__CINT__) || defined(__MAKECINT__)
#include <TCanvas.h>
#include <TClass.h>
#include <TDirectory.h>
#include <TFile.h>
#include <TH1.h>
#include <TH2.h>
#include <TKey.h>
#include <TLegend.h>
#include <TPad.h>
#include <TROOT.h>
#include <TStyle.h>
#include <TSystem.h>

#include <map>
#include <memory>
#include <string>
#include <vector>
#endif

namespace
{

std::string sanitizeName(std::string name);

void collect1DHistograms(TDirectory* dir, const std::string& prefix, std::map<std::string, std::string>& histPaths)
{
  if (!dir) {
    return;
  }

  TIter next(dir->GetListOfKeys());
  TKey* key = nullptr;
  while ((key = static_cast<TKey*>(next()))) {
    TObject* obj = key->ReadObj();
    if (!obj) {
      continue;
    }

    std::string name = obj->GetName();
    std::string path = prefix.empty() ? name : prefix + "/" + name;

    if (obj->InheritsFrom(TDirectory::Class())) {
      collect1DHistograms(static_cast<TDirectory*>(obj), path, histPaths);
    } else if (obj->InheritsFrom(TH1::Class()) && !obj->InheritsFrom(TH2::Class())) {
      histPaths[path] = path;
    }

    delete obj;
  }
}

TH1* loadHistogram(TFile* file, const std::string& path, const char* tag)
{
  if (!file) {
    return nullptr;
  }
  TObject* obj = file->Get(path.c_str());
  if (!obj || !obj->InheritsFrom(TH1::Class()) || obj->InheritsFrom(TH2::Class())) {
    return nullptr;
  }
  std::string cloneName = sanitizeName(path) + "_" + tag;
  TH1* hist = static_cast<TH1*>(obj->Clone(cloneName.c_str()));
  hist->SetDirectory(nullptr);
  if (hist->GetListOfFunctions()) {
    hist->GetListOfFunctions()->Clear();
  }
  return hist;
}

std::string sanitizeName(std::string name)
{
  for (char& ch : name) {
    if (!(std::isalnum(static_cast<unsigned char>(ch)) || ch == '_' || ch == '-')) {
      ch = '_';
    }
  }
  return name;
}

void styleHistogram(TH1* hist, int color)
{
  if (!hist) {
    return;
  }
  hist->SetLineColor(color);
  hist->SetMarkerColor(color);
  hist->SetLineWidth(2);
  hist->SetMarkerStyle(20);
  hist->SetMarkerSize(0.9);
  hist->SetStats(0);
}

void computeOverlayYRange(const TH1* h1, const TH1* h2, double& ymin, double& ymax)
{
  constexpr double kHuge = 1.e300;
  ymin = kHuge;
  ymax = -kHuge;

  auto updateRange = [&](const TH1* h) {
    if (!h) {
      return;
    }
    for (int ib = 1; ib <= h->GetNbinsX(); ++ib) {
      double c = h->GetBinContent(ib);
      double e = h->GetBinError(ib);
      ymin = std::min(ymin, c - e);
      ymax = std::max(ymax, c + e);
    }
  };

  updateRange(h1);
  updateRange(h2);

  if (ymin == kHuge || ymax == -kHuge) {
    ymin = 0.;
    ymax = 1.;
    return;
  }

  if (ymax <= ymin) {
    ymax = ymin + 1.;
  }

  double span = ymax - ymin;
  ymin -= 0.08 * span;
  ymax += 0.12 * span;

  if (ymin > 0.) {
    ymin = 0.;
  }
}

} // namespace

void comparePlotTracks(const char* file1 = "../jpsi_bkg_40gev/TrackingPerformance.root",
                       const char* file2 = "../jpsi_bkg_40gev_plate/TrackingPerformance.root",
                       const char* label1 = "No plate",
                       const char* label2 = "Plate",
                       const char* outDir = "compare_plotTracks", int requiredHits = -1)
{
  gStyle->SetOptStat(0);
  if (requiredHits < 0) {
    comparePlotTracks("../jpsi_bkg_40gev/TrackingPerformance_nHits3.root", "../jpsi_bkg_40gev_plate/TrackingPerformance_nHits3.root", label1, label2, "compare_plotTracks_nHits3", 3);
    comparePlotTracks("../jpsi_bkg_40gev/TrackingPerformance_nHits4.root", "../jpsi_bkg_40gev_plate/TrackingPerformance_nHits4.root", label1, label2, "compare_plotTracks_nHits4", 4);
    comparePlotTracks("../jpsi_bkg_40gev/TrackingPerformance_nHits5.root", "../jpsi_bkg_40gev_plate/TrackingPerformance_nHits5.root", label1, label2, "compare_plotTracks_nHits5", 5);
    return;
  }

  std::unique_ptr<TFile> in1(TFile::Open(file1));
  std::unique_ptr<TFile> in2(TFile::Open(file2));
  if (!in1 || in1->IsZombie()) {
    printf("Cannot open %s\n", file1);
    return;
  }
  if (!in2 || in2->IsZombie()) {
    printf("Cannot open %s\n", file2);
    return;
  }

  gSystem->mkdir(outDir, kTRUE);
  std::string rootOutName = std::string(outDir) + "/comparison.root";
  std::unique_ptr<TFile> outFile(TFile::Open(rootOutName.c_str(), "RECREATE"));

  std::map<std::string, std::string> histPaths1;
  std::map<std::string, std::string> histPaths2;
  collect1DHistograms(in1.get(), "", histPaths1);
  collect1DHistograms(in2.get(), "", histPaths2);

  int compared = 0;
  int missingInSecond = 0;

  for (const auto& [path, _] : histPaths1) {
    auto it2 = histPaths2.find(path);
    if (it2 == histPaths2.end()) {
      missingInSecond++;
      continue;
    }

    std::unique_ptr<TH1> h1(loadHistogram(in1.get(), path, "f1"));
    std::unique_ptr<TH1> h2(loadHistogram(in2.get(), path, "f2"));
    if (!h1 || !h2) {
      continue;
    }
    const TAxis* x1 = h1->GetXaxis();
    const TAxis* x2 = h2->GetXaxis();
    bool sameBinning = (h1->GetNbinsX() == h2->GetNbinsX()) &&
                       (x1->GetXmin() == x2->GetXmin()) &&
                       (x1->GetXmax() == x2->GetXmax());
    if (sameBinning) {
      for (int ib = 1; ib <= h1->GetNbinsX(); ++ib) {
        double e1 = x1->GetBinLowEdge(ib);
        double e2 = x2->GetBinLowEdge(ib);
        double w1 = x1->GetBinWidth(ib);
        double w2 = x2->GetBinWidth(ib);
        if (e1 != e2 || w1 != w2) {
          sameBinning = false;
          break;
        }
      }
    }
    if (!sameBinning) {
      printf("Skipping %s: different X axis definition\n", path.c_str());
      continue;
    }

    compared++;
    styleHistogram(h1.get(), kAzure + 2);
    styleHistogram(h2.get(), kOrange + 7);

    double ymin = 0.;
    double ymax = 1.;
    computeOverlayYRange(h1.get(), h2.get(), ymin, ymax);
    h1->SetMinimum(ymin);
    h1->SetMaximum(ymax);

    std::unique_ptr<TCanvas> canvas(new TCanvas(("c_" + sanitizeName(path)).c_str(), path.c_str(), 1000, 900));
    canvas->Divide(1, 2);

    TPad* pad1 = static_cast<TPad*>(canvas->cd(1));
    pad1->SetPad(0., 0.32, 1., 1.);
    pad1->SetBottomMargin(0.02);
    pad1->SetTickx();
    pad1->SetTicky();

    h1->GetXaxis()->SetLabelSize(0.0);
    h1->GetXaxis()->SetTitleSize(0.0);

    // h1->Draw("hist");
    h1->Draw("e1 p");
    // h2->Draw("hist same");
    h2->Draw("e1 p same");

    TLegend legend(0.62, 0.73, 0.88, 0.88);
    legend.SetBorderSize(0);
    legend.SetFillStyle(0);
    legend.AddEntry(h1.get(), label1, "lep");
    legend.AddEntry(h2.get(), label2, "lep");
    legend.Draw();

    TPad* pad2 = static_cast<TPad*>(canvas->cd(2));
    pad2->SetPad(0., 0., 1., 0.32);
    pad2->SetTopMargin(0.03);
    pad2->SetBottomMargin(0.28);
    pad2->SetTickx();
    pad2->SetTicky();

    std::unique_ptr<TH1> ratio(static_cast<TH1*>(h2->Clone((sanitizeName(path) + "_ratio").c_str())));
    ratio->SetDirectory(nullptr);
    // ratio->Divide(h1.get());
    for (int ib = 1; ib <= ratio->GetNbinsX(); ++ib) {
      double n = ratio->GetBinContent(ib);
      double d = h1->GetBinContent(ib);
      double c = (d != 0.) ? n / d : 0.;
      double e = (d != 0. && n != 0.) ? c*std::sqrt(1. / d + 1. / n) : 0.;

      ratio->SetBinContent(ib, c);
      ratio->SetBinError(ib, 0*e);
    }
    ratio->SetTitle("");
    ratio->GetYaxis()->SetTitle("Ratio");
    ratio->GetYaxis()->SetNdivisions(505);
    ratio->GetYaxis()->SetTitleSize(0.10);
    ratio->GetYaxis()->SetTitleOffset(0.45);
    ratio->GetYaxis()->SetLabelSize(0.08);
    ratio->GetXaxis()->SetTitle(h1->GetXaxis()->GetTitle());
    ratio->GetXaxis()->SetTitleSize(0.10);
    ratio->GetXaxis()->SetTitleOffset(1.0);
    ratio->GetXaxis()->SetLabelSize(0.08);
    ratio->SetLineColor(kBlack);
    ratio->SetMarkerColor(kBlack);
    ratio->SetMinimum(0.);
    ratio->SetMaximum(2);
    ratio->Draw("hist e");

    if (outFile) {
      outFile->cd();
      canvas->Write();
      h1->Write((sanitizeName(path) + "_1").c_str());
      h2->Write((sanitizeName(path) + "_2").c_str());
      ratio->Write((sanitizeName(path) + "_ratio").c_str());
    }

    std::string pngName = std::string(outDir) + "/" + sanitizeName(path) + ".png";
    canvas->SaveAs(pngName.c_str());
  }

  printf("Compared %d histograms\n", compared);
  printf("Missing in second file: %d histograms\n", missingInSecond);
  printf("Output written to %s and %s\n", outDir, rootOutName.c_str());
}
