#include <TCanvas.h>
#include <TFile.h>
#include <TH1.h>
#include <TLegend.h>
#include <TPad.h>
#include <TString.h>
#include <TStyle.h>
#include <TSystem.h>

#include <algorithm>
#include <cerrno>
#include <cstdlib>
#include <dirent.h>
#include <iostream>
#include <memory>
#include <vector>

// Run from the repository root:
// root -l -b -q 'compareDaikiSingleMuonEfficiency.C()'
void compareDaikiSingleMuonEfficiency(
    const char* baseDir = "debug_daiki/layout.medAbsorber[33]=Concrete")
{
  const int colors[] = {kBlue + 1, kRed + 1, kGreen + 2, kMagenta + 1,
                        kOrange + 7, kCyan + 2, kViolet + 1, kGray + 2};
  constexpr int nColors = sizeof(colors) / sizeof(colors[0]);
  const char* prefixes[] = {"reco.msMaxDeltaTanLCellsCA[0]=",
                            "reco.vtMaxDeltaTanLCellsCA[0]="};
  const char* labels[] = {"MS", "VT"};
  struct Selection {
    int cut;
    int type;
    TString directory;
  };
  const char* variables[] = {"Rapidity", "RapidityCM", "Momentum"};
  const char* xTitles[] = {"MC rapidity", "MC y-y_{CM}", "MC p (GeV/c)"};
  const char* measures[] = {"hRecoEffVs", "hEffTimesAccVs"};
  const char* measureNames[] = {"reconstruction", "efficiency_times_acceptance"};
  const char* yTitles[] = {"N_{reco}/N_{reconstructable}", "N_{reco}/N_{generated}"};

  std::vector<Selection> selections;
  DIR* directory = opendir(baseDir);
  if (!directory) {
    std::cerr << "Cannot scan " << baseDir << std::endl;
    return;
  }
  while (auto* entry = readdir(directory)) {
    const TString name(entry->d_name);
    int type = -1;
    for (int i = 0; i < 2; ++i) {
      if (name.BeginsWith(prefixes[i])) { type = i; break; }
    }
    if (type < 0) continue;
    const char* value = entry->d_name + TString(prefixes[type]).Length();
    char* end = nullptr;
    errno = 0;
    const long cut = std::strtol(value, &end, 10);
    if (errno || end == value || *end != '\0' || cut < 0) continue;
    const TString path = Form("%s/%s/dimuon_analysisMS.root", baseDir, entry->d_name);
    if (!gSystem->AccessPathName(path))
      selections.push_back({static_cast<int>(cut), type, name});
  }
  closedir(directory);
  std::sort(selections.begin(), selections.end(), [](const Selection& a, const Selection& b) {
    return a.type == b.type ? a.cut < b.cut : a.type < b.type;
  });
  if (selections.empty()) {
    std::cerr << "No cut directories with dimuon_analysisMS.root found in " << baseDir << std::endl;
    return;
  }

  std::vector<std::unique_ptr<TFile>> inputs;
  for (size_t i = 0; i < selections.size(); ++i) {
    const TString path = Form("%s/%s/dimuon_analysisMS.root",
                              baseDir, selections[i].directory.Data());
    inputs.emplace_back(TFile::Open(path, "READ"));
    if (!inputs.back() || inputs.back()->IsZombie()) {
      std::cerr << "Cannot open " << path << std::endl;
      return;
    }
  }

  gStyle->SetOptStat(0);
  for (int measure = 0; measure < 2; ++measure) {
    for (int variable = 0; variable < 3; ++variable) {
      // Keep detached histograms alive until after SaveAs.
      std::vector<std::unique_ptr<TH1>> curves;
      TCanvas canvas(Form("cCompare%d_%d", measure, variable),
                     "Daiki single-muon efficiency comparison", 1000, 700);
      {
        canvas.cd();
        gPad->SetGridy();
        auto legend = new TLegend(0.48, 0.14, 0.88, 0.34);
        legend->SetBorderSize(0);
        legend->SetFillStyle(0);
        for (size_t selection = 0; selection < selections.size(); ++selection) {
          const TString name = Form("%s%s", measures[measure], variables[variable]);
          auto* source = dynamic_cast<TH1*>(inputs[selection]->Get(name));
          if (!source) {
            std::cerr << "Missing " << name << " in " << inputs[selection]->GetName() << std::endl;
            return;
          }
          auto* hist = dynamic_cast<TH1*>(source->Clone(
              Form("%s_type%d_cut%d", name.Data(), selections[selection].type,
                   selections[selection].cut)));
          hist->SetDirectory(nullptr);
          hist->SetStats(false);
          hist->SetLineColor(colors[selection % nColors]);
          hist->SetMarkerColor(colors[selection % nColors]);
          hist->SetMarkerStyle(20 + selection % 10);
          hist->SetMarkerSize(1.1);
          hist->SetTitle(Form("#mu^{+} + #mu^{-};%s;%s",
                              xTitles[variable], yTitles[measure]));
          hist->SetMinimum(0.);
          hist->SetMaximum(1.05);
          hist->Draw(selection == 0 ? "E1" : "E1 SAME");
          legend->AddEntry(hist, Form("%s CA tanL cut = %d",
                                      labels[selections[selection].type],
                                      selections[selection].cut), "lp");
          curves.emplace_back(hist);
        }
        legend->Draw();
      }
      const TString output = Form("%s/compare_single_muon_%s_vs_%s.png", baseDir,
                                  measureNames[measure], variables[variable]);
      canvas.SaveAs(output);
      std::cout << "Saved " << output << std::endl;
    }
  }
}
