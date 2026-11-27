#include "ThermalMuonMSAnalysis.h"
#include "NA6PMCComposedLabel.h"
#include "NA6PMatch.h"
#include "NA6PTrack.h"
#include <TCanvas.h>
#include <TEfficiency.h>
#include <TFile.h>
#include <TGaxis.h>
#include <TGraph.h>
#include <TGraphAsymmErrors.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TLatex.h>
#include <TLegend.h>
#include <TStyle.h>
#include <TSystem.h>
#include <TTree.h>
#include <TTreeReader.h>
#include <TTreeReaderValue.h>
#include <algorithm>
#include <fstream>
#include <iostream>

R__LOAD_LIBRARY(librecLib)

namespace
{
struct Observation {
  double energy;
  bool fake;
};

struct JointObservation {
  double loss;
  double residual;
  bool fake;
};

void drawFakeRate(const std::string& out, std::ostream& csv, const char* name, const char* title,
                  const char* axisTitle, const char* sample, const std::vector<Observation>& observations,
                  std::vector<double> edges, bool showFakeCounts = false)
{
  if (observations.empty()) {
    throw std::runtime_error(std::string("Empty sample for ") + name);
  }
  const auto [minimum, maximum] = std::minmax_element(observations.begin(), observations.end(),
    [](const auto& a, const auto& b) { return a.energy < b.energy; });
  if (minimum->energy < edges.front()) {
    edges.insert(edges.begin(), std::floor(minimum->energy));
  }
  if (maximum->energy >= edges.back()) {
    edges.push_back(std::floor(maximum->energy) + 1.);
  }
  TEfficiency rate(name, "", edges.size() - 1, edges.data());
  rate.SetConfidenceLevel(.68268949);
  rate.SetStatisticOption(TEfficiency::kFCP);
  int fakes = 0;
  for (const auto& point : observations) {
    rate.Fill(point.fake, point.energy);
    fakes += point.fake;
  }

  TGraphAsymmErrors graph;
  graph.SetLineColor(kAzure + 3);
  graph.SetMarkerColor(kAzure + 3);
  graph.SetMarkerStyle(kFullCircle);
  graph.SetMarkerSize(1.15);
  graph.SetLineWidth(2);
  const auto& all = *rate.GetTotalHistogram();
  const auto& fake = *rate.GetPassedHistogram();
  constexpr double rateMax = 1.04;
  const double countMax = std::max(10., std::ceil(1.2 * fake.GetMaximum() / 10.) * 10.);
  constexpr int countColor = kOrange + 8;
  TGraph fakeCounts;
  fakeCounts.SetMarkerColor(countColor);
  fakeCounts.SetMarkerStyle(kOpenSquare);
  fakeCounts.SetMarkerSize(1.25);
  csv << name << ",all,,," << fakes << ',' << observations.size() << ',' << double(fakes) / observations.size() << ",,\n";
  for (int bin = 1; bin <= all.GetNbinsX(); ++bin) {
    const double low = all.GetXaxis()->GetBinLowEdge(bin), high = all.GetXaxis()->GetBinUpEdge(bin);
    const int total = all.GetBinContent(bin), passed = fake.GetBinContent(bin);
    csv << name << ",bin," << low << ',' << high << ',' << passed << ',' << total << ',';
    if (total == 0) {
      csv << ",,\n";
      continue;
    }
    const double value = rate.GetEfficiency(bin);
    const double errLow = rate.GetEfficiencyErrorLow(bin), errHigh = rate.GetEfficiencyErrorUp(bin);
    const int point = graph.GetN();
    graph.SetPoint(point, (low + high) / 2., value);
    graph.SetPointError(point, (high - low) / 2., (high - low) / 2., errLow, errHigh);
    if (showFakeCounts) {
      fakeCounts.SetPoint(point, (low + high) / 2., rateMax * passed / countMax);
    }
    csv << value << ',' << value - errLow << ',' << value + errHigh << '\n';
  }

  TCanvas canvas(name, title, 1050, 725);
  canvas.SetMargin(.12, showFakeCounts ? .12 : .05, .14, .19);
  canvas.SetTicks(1, !showFakeCounts);
  canvas.SetGridy();
  auto* frame = canvas.DrawFrame(edges.front(), 0., edges.back(), rateMax, Form(";%s;Fake rate", axisTitle));
  frame->GetXaxis()->CenterTitle();
  frame->GetXaxis()->SetTitleSize(.043);
  frame->GetYaxis()->SetTitleSize(.043);
  graph.Draw("P SAME");
  TGaxis countAxis(edges.back(), 0., edges.back(), rateMax, 0., countMax, 506, "+L");
  TLegend legend(.60, .30, .86, .415);
  if (showFakeCounts) {
    fakeCounts.Draw("P SAME");
    countAxis.SetTitle("Fake matches / bin");
    countAxis.SetLineColor(countColor);
    countAxis.SetLabelColor(countColor);
    countAxis.SetTitleColor(countColor);
    countAxis.SetLabelFont(frame->GetYaxis()->GetLabelFont());
    countAxis.SetTitleFont(frame->GetYaxis()->GetTitleFont());
    countAxis.SetLabelSize(frame->GetYaxis()->GetLabelSize());
    countAxis.SetTitleSize(frame->GetYaxis()->GetTitleSize());
    countAxis.SetTitleOffset(1.17);
    countAxis.Draw();
    legend.SetBorderSize(0);
    legend.SetFillStyle(0);
    legend.SetTextSize(.029);
    legend.AddEntry(&graph, "Fake rate (left axis)", "pe");
    legend.AddEntry(&fakeCounts, "Fake count (right axis)", "p");
    legend.Draw();
  }
  TLatex label;
  label.SetNDC();
  label.SetTextFont(62);
  label.SetTextSize(.041);
  label.DrawLatex(.12, .945, title);
  label.SetTextFont(42);
  label.SetTextSize(.029);
  label.DrawLatex(.12, .897, Form("VT threshold #geq5 hits  |  %s", sample));
  label.SetTextSize(.025);
  label.DrawLatex(.12, .855, Form("%d / %zu fake overall  |  68%% binomial intervals", fakes, observations.size()));
  canvas.SaveAs((out + "/" + name + ".png").c_str());
  std::cout << name << " fake=" << fakes << '/' << observations.size()
            << " min=" << minimum->energy << " max=" << maximum->energy << '\n';
}

void drawJointFakeRate(const std::string& out, const std::vector<JointObservation>& observations)
{
  const std::vector<double> lossCuts{1.25, 1.35, 1.45, 1.6, 2., 3.};
  const std::vector<double> residualCuts{-1., -.5, -.25, 0., .25, .5, 1.};
  const int nx = lossCuts.size() + 1, ny = residualCuts.size() + 1;
  const char* name = "fakeRate_vs_MS_energyResidual_and_energyLoss";
  TH2D heatmap(name, ";Estimated MS energy loss (GeV);E_{gen} #minus E_{rec}^{MS}(PV) (GeV);Fake rate",
               nx, 0., nx, ny, 0., ny);
  heatmap.SetDirectory(nullptr);
  heatmap.SetStats(false);
  TEfficiency rate("jointRate", "", nx, 0., nx, ny, 0., ny);
  rate.SetConfidenceLevel(.68268949);
  rate.SetStatisticOption(TEfficiency::kFCP);
  int fakes = 0;
  for (const auto& point : observations) {
    const auto x = std::upper_bound(lossCuts.begin(), lossCuts.end(), point.loss) - lossCuts.begin();
    const auto y = std::upper_bound(residualCuts.begin(), residualCuts.end(), point.residual) - residualCuts.begin();
    rate.Fill(point.fake, x + .5, y + .5);
    fakes += point.fake;
  }
  const auto labelBins = [](TAxis& axis, const std::vector<double>& cuts) {
    axis.SetBinLabel(1, Form("<%g", cuts.front()));
    for (size_t i = 1; i < cuts.size(); ++i) {
      axis.SetBinLabel(i + 1, Form("%g to %g", cuts[i - 1], cuts[i]));
    }
    axis.SetBinLabel(cuts.size() + 1, Form("#geq%g", cuts.back()));
    axis.CenterTitle();
    axis.SetTitleSize(.036);
    axis.SetLabelSize(.029);
  };
  labelBins(*heatmap.GetXaxis(), lossCuts);
  labelBins(*heatmap.GetYaxis(), residualCuts);
  heatmap.GetYaxis()->SetTitleOffset(2.05);
  heatmap.GetZaxis()->SetTitleOffset(1.3);

  std::ofstream csv(out + "/fakeRateMSEnergy2D.csv");
  if (!csv) {
    throw std::runtime_error("Cannot create 2D fake-rate table");
  }
  csv << "loss_low_GeV,loss_high_GeV,residual_low_GeV,residual_high_GeV,fake_matches,all_matches,fake_rate,ci68_low,ci68_high\n";
  const auto writeEdges = [&csv](const std::vector<double>& cuts, int bin) {
    if (bin > 1) csv << cuts[bin - 2];
    else csv << "-inf";
    csv << ',';
    if (bin <= static_cast<int>(cuts.size())) csv << cuts[bin - 1];
    else csv << "inf";
    csv << ',';
  };
  const auto& all = *rate.GetTotalHistogram();
  const auto& fake = *rate.GetPassedHistogram();
  for (int iy = 1; iy <= ny; ++iy) {
    for (int ix = 1; ix <= nx; ++ix) {
      const int bin = rate.GetGlobalBin(ix, iy);
      const int total = all.GetBinContent(bin), passed = fake.GetBinContent(bin);
      writeEdges(lossCuts, ix);
      writeEdges(residualCuts, iy);
      csv << passed << ',' << total << ',';
      if (total == 0) {
        csv << ",,\n";
        continue;
      }
      const double value = rate.GetEfficiency(bin);
      heatmap.SetBinContent(ix, iy, std::max(value, 1.e-9)); // Paint occupied zero-rate cells.
      csv << value << ',' << value - rate.GetEfficiencyErrorLow(bin) << ','
          << value + rate.GetEfficiencyErrorUp(bin) << '\n';
    }
  }
  TCanvas canvas(name, "Fake rate vs. MS energy residual and loss", 1150, 800);
  canvas.SetMargin(.155, .14, .13, .16);
  gStyle->SetPalette(kBird);
  heatmap.SetMinimum(0.);
  heatmap.SetMaximum(1.);
  heatmap.Draw("COLZ");
  TLatex label;
  label.SetTextAlign(22);
  label.SetTextSize(.028);
  for (int iy = 1; iy <= ny; ++iy) {
    for (int ix = 1; ix <= nx; ++ix) {
      const int bin = rate.GetGlobalBin(ix, iy);
      const int total = all.GetBinContent(bin);
      if (total == 0) continue;
      label.SetTextColor(rate.GetEfficiency(bin) < .48 ? kWhite : kBlack);
      label.DrawLatex(ix - .5, iy - .5, Form("%d / %d", static_cast<int>(fake.GetBinContent(bin)), total));
    }
  }
  label.SetNDC();
  label.SetTextAlign(11);
  label.SetTextColor(kBlack);
  label.SetTextFont(62);
  label.SetTextSize(.039);
  label.DrawLatex(.155, .949, "Fake rate vs. MS energy residual and loss");
  label.SetTextFont(42);
  label.SetTextSize(.027);
  label.DrawLatex(.155, .904, Form("VT threshold #geq5 hits  |  %d / %zu fake overall", fakes, observations.size()));
  label.DrawLatex(.155, .868, "Cell labels show fake / total matches");
  canvas.SaveAs((out + "/" + name + ".png").c_str());
}
}

void plotFakeRateMSEnergy(const char* runDir = ".")
{
  const std::string run = runDir;
  const auto options = thermalMuonMS::loadPropagation(run);
  const auto* propagator = Propagator::Instance();
  const double mass = PID::getMass(PID::Muon);
  TFile msFile((run + "/reco/TracksMuonSpec.root").c_str());
  TFile vtFile((run + "/reco/TracksVerTel.root").c_str());
  TFile mcFile((run + "/sim/MCKine.root").c_str());
  TFile matchFile((run + "/reco_5hits/TracksMatching.root").c_str());
  TTree* msTree = msFile.Get<TTree>("tracksMuonSpec");
  TTree* vtTree = vtFile.Get<TTree>("tracksVerTel");
  TTree* mcTree = mcFile.Get<TTree>("mckine");
  TTree* matchTree = matchFile.Get<TTree>("tracksMatching");
  if (!msTree || !vtTree || !mcTree || !matchTree || msTree->GetEntries() != mcTree->GetEntries() ||
      vtTree->GetEntries() != mcTree->GetEntries() || matchTree->GetEntries() != mcTree->GetEntries()) {
    throw std::runtime_error("Missing or inconsistent MS, VT, MC or matching trees");
  }
  TTreeReader ms(msTree), vt(vtTree), mc(mcTree), matching(matchTree);
  TTreeReaderValue<std::vector<NA6PTrack>> msTracks(ms, "MuonSpec"), vtTracks(vt, "VerTel");
  TTreeReaderValue<std::vector<NA6PMCComposedLabel>> msLabels(ms, "MuonSpecMCTruth"), vtLabels(vt, "VerTelMCTruth");
  TTreeReaderValue<std::vector<NA6PMatch>> matches(matching, "Matching");
  TTreeReaderValue<NA6PMCEventHeader> header(mc, "header");
  TTreeReaderValue<std::vector<TParticle>> particles(mc, "tracks");

  std::vector<Observation> energyLoss, residual, thermalResidual;
  std::vector<JointObservation> joint;
  Long64_t events = 0;
  int allMatches = 0, allFake = 0, unavailableLoss = 0, unavailableResidual = 0;
  while (mc.Next()) {
    if (!ms.Next() || !vt.Next() || !matching.Next() ||
        msTracks->size() != msLabels->size() || vtTracks->size() != vtLabels->size()) {
      throw std::runtime_error("Inconsistent event data at entry " + std::to_string(events));
    }
    ++events;
    const float zPV = header->getVZ();
    if (!std::isfinite(zPV)) {
      throw std::runtime_error("Invalid primary vertex z");
    }
    const auto thermal = thermalMuonMS::findThermalMuons(*header, *particles);
    const auto genuineID = [&](const NA6PMCComposedLabel& label) {
      if (!label.isValid() || label.isFake()) {
        return -1;
      }
      const int id = label.getTrackID();
      if (label.getSourceID() != 0 || label.getEventID() != header->getEventID() ||
          id < 0 || static_cast<size_t>(id) >= particles->size()) {
        throw std::runtime_error("Genuine track label outside the corresponding MC event");
      }
      return id;
    };
    for (const auto& match : *matches) {
      const int mi = match.getIndexMS(), vi = match.getIndexVT();
      if (mi < 0 || vi < 0 || static_cast<size_t>(mi) >= msTracks->size() ||
          static_cast<size_t>(vi) >= vtTracks->size() || vtTracks->at(vi).getNHits() < 5) {
        throw std::runtime_error("Invalid index or VT hit count in threshold-5 matches");
      }
      const int id = genuineID(msLabels->at(mi));
      const bool isThermal = id >= 0 && thermal[id];
      const bool fake = !(isThermal && genuineID(vtLabels->at(vi)) == id);
      ++allMatches;
      allFake += fake;

      const auto& track = msTracks->at(mi);
      const double pMS = track.getP();
      NA6PTrackParCov atPV = track;
      atPV.setPID(PID::Muon);
      if (!track.isValid() || !std::isfinite(pMS) || pMS <= 0. ||
          !propagator->propagateToZ(atPV, zPV, options) || atPV.getZ() != zPV ||
          !std::isfinite(atPV.getP()) || atPV.getP() <= 0.) {
        ++unavailableLoss;
        ++unavailableResidual;
        continue;
      }
      const double ePV = std::hypot(atPV.getP(), mass);
      const double loss = ePV - std::hypot(pMS, mass);
      energyLoss.push_back({loss, fake});
      if (id < 0 || !std::isfinite(particles->at(id).Energy())) {
        ++unavailableResidual;
        continue;
      }
      const Observation entry{particles->at(id).Energy() - ePV, fake};
      residual.push_back(entry);
      joint.push_back({loss, entry.energy, fake});
      if (isThermal) {
        thermalResidual.push_back(entry);
      }
    }
  }
  if (allMatches == 0 || events != mcTree->GetEntries() ||
      energyLoss.size() + unavailableLoss != static_cast<size_t>(allMatches) ||
      residual.size() + unavailableResidual != static_cast<size_t>(allMatches)) {
    throw std::runtime_error("Incomplete event read or inconsistent match counts");
  }

  const std::string out = run + "/fake";
  if (gSystem->mkdir(out.c_str(), true) != 0 && gSystem->AccessPathName(out.c_str(), kWritePermission)) {
    throw std::runtime_error("Cannot create fake-rate directory");
  }
  std::ofstream csv(out + "/fakeRateMSEnergy.csv");
  if (!csv) {
    throw std::runtime_error("Cannot create fake-rate table");
  }
  csv << "observable,scope,low_GeV,high_GeV,fake_matches,all_matches,fake_rate,ci68_low,ci68_high\n"
      << "threshold5,all,,," << allFake << ',' << allMatches << ',' << double(allFake) / allMatches << ",,\n";
  const std::vector<double> residualEdges{-4., -2., -1., -.5, -.25, 0., .25, .5, 1., 2., 4.};
  const char* residualAxis = "E_{gen} #minus E_{rec}^{MS}(PV) (GeV)";
  drawFakeRate(out, csv, "fakeRate_vs_MS_energyLoss", "Fake rate vs. estimated MS energy loss",
               "E_{rec}^{MS}(PV) #minus E_{rec}^{MS}(MS plane) (GeV)", "All matches with an energy-loss estimate",
               energyLoss, {1., 1.25, 1.30, 1.35, 1.40, 1.45, 1.50, 1.75, 2., 2.25, 2.5, 3.}, true);
  drawFakeRate(out, csv, "fakeRate_vs_MS_energyResidual", "Fake rate vs. MS energy residual",
               residualAxis, "Matches with an unambiguous MS truth label", residual, residualEdges);
  drawFakeRate(out, csv, "fakeRate_vs_MS_energyResidual_thermal", "Fake rate vs. MS energy residual",
               residualAxis, "Thermal MS muons with an unambiguous truth label", thermalResidual, residualEdges);
  drawJointFakeRate(out, joint);
  std::cout << "events=" << events << " fake=" << allFake << '/' << allMatches
            << " unavailableLoss=" << unavailableLoss << " unavailableResidual=" << unavailableResidual << '\n';
}
