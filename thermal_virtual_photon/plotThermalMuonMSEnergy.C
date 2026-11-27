#include "ThermalMuonMSAnalysis.h"
#include "NA6PMCComposedLabel.h"
#include "NA6PTrack.h"
#include <TCanvas.h>
#include <TFile.h>
#include <TGraphAsymmErrors.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TLegend.h>
#include <TLatex.h>
#include <TLine.h>
#include <TParticle.h>
#include <TSystem.h>
#include <TTree.h>
#include <TTreeReader.h>
#include <TTreeReaderValue.h>
#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>
#include <map>
#include <numeric>
#include <stdexcept>
#include <string>
#include <vector>

R__LOAD_LIBRARY(librecLib)

void plotThermalMuonMSEnergy(const char* runDir = ".")
{
  const std::string run = runDir;
  const auto propagation = thermalMuonMS::loadPropagation(run);
  const auto* propagator = Propagator::Instance();
  TFile reco((run + "/reco/TracksMuonSpec.root").c_str());
  TFile kine((run + "/sim/MCKine.root").c_str());
  TTree* msTree = reco.Get<TTree>("tracksMuonSpec");
  TTree* mcTree = kine.Get<TTree>("mckine");
  if (!msTree || !mcTree || msTree->GetEntries() != mcTree->GetEntries()) {
    throw std::runtime_error("Missing or inconsistent MS and MC trees");
  }
  TTreeReader ms(msTree), mc(mcTree);
  TTreeReaderValue<std::vector<NA6PTrack>> tracks(ms, "MuonSpec");
  TTreeReaderValue<std::vector<NA6PMCComposedLabel>> labels(ms, "MuonSpecMCTruth");
  TTreeReaderValue<NA6PMCEventHeader> header(mc, "header");
  TTreeReaderValue<std::vector<TParticle>> particles(mc, "tracks");
  const double muonMass = PID::getMass(PID::Muon);
  std::vector<double> residuals;
  std::vector<std::pair<double, double>> residualsByY;
  std::map<int, int> tracksByHits;
  double minZ = std::numeric_limits<double>::infinity();
  double maxZ = -minZ;
  Long64_t events = 0;
  int invalidEnergy = 0;
  int failedPropagation = 0;
  while (mc.Next()) {
    if (!ms.Next() || tracks->size() != labels->size()) {
      throw std::runtime_error("Inconsistent MS tracks and labels at entry " + std::to_string(events));
    }
    ++events;
    const float zPV = header->getVZ();
    if (!std::isfinite(zPV)) {
      throw std::runtime_error("Invalid primary vertex z");
    }
    const auto thermal = thermalMuonMS::findThermalMuons(*header, *particles);

    for (size_t i = 0; i < tracks->size(); ++i) {
      const auto& label = labels->at(i);
      if (!label.isValid() || label.isFake()) {
        continue;
      }
      const int id = label.getTrackID();
      if (label.getSourceID() != 0 || label.getEventID() != header->getEventID() ||
          id < 0 || static_cast<size_t>(id) >= thermal.size()) {
        throw std::runtime_error("Genuine MS label does not identify a particle in the corresponding MC event");
      }
      if (!thermal[id]) {
        continue;
      }
      const auto& track = tracks->at(i);
      const auto& particle = particles->at(id);
      const double eGen = particle.Energy();
      if (!track.isValid() || !std::isfinite(track.getP()) || track.getP() <= 0. ||
          !std::isfinite(eGen) || eGen < muonMass) {
        ++invalidEnergy;
        continue;
      }
      NA6PTrackParCov atPV = track;
      atPV.setPID(PID::Muon);
      if (!propagator->propagateToZ(atPV, zPV, propagation) || atPV.getZ() != zPV) {
        ++failedPropagation;
        continue;
      }
      const double pRec = atPV.getP();
      if (!std::isfinite(pRec) || pRec <= 0.) {
        ++invalidEnergy;
        continue;
      }
      const double residual = eGen - std::hypot(pRec, muonMass);
      const double yGen = std::asinh(particle.Pz() / std::hypot(particle.Pt(), muonMass));
      if (!std::isfinite(yGen)) {
        throw std::runtime_error("Invalid generated muon rapidity");
      }
      residuals.push_back(residual);
      residualsByY.emplace_back(yGen, residual);
      ++tracksByHits[track.getNHits()];
      minZ = std::min(minZ, static_cast<double>(atPV.getZ()));
      maxZ = std::max(maxZ, static_cast<double>(atPV.getZ()));
    }
  }
  if (events != mcTree->GetEntries() || residuals.empty()) {
    throw std::runtime_error("Incomplete MC read or no thermal MS tracks");
  }

  std::sort(residuals.begin(), residuals.end());
  const auto quantile = [](const std::vector<double>& values, double fraction) {
    const double index = fraction * (values.size() - 1);
    const size_t low = index;
    const size_t high = std::min(low + 1, values.size() - 1);
    return values[low] + (index - low) * (values[high] - values[low]);
  };
  const double mean = std::accumulate(residuals.begin(), residuals.end(), 0.) / residuals.size();
  const double median = quantile(residuals, .5);
  const std::string out = run + "/thermal_muon_ms_energy";
  if (gSystem->mkdir(out.c_str(), true) != 0 && gSystem->AccessPathName(out.c_str(), kWritePermission)) {
    throw std::runtime_error("Cannot create plot directory");
  }

  const char* residualTitle = "E_{gen} #minus E_{rec}^{MS} at the primary vertex (GeV)";
  const auto draw = [&](const char* name, const char* view, int bins, double low, double high, bool logY) {
    TH1D hist(name, Form(";%s;Tracks / bin", residualTitle), bins, low, high);
    hist.SetDirectory(nullptr);
    hist.SetStats(false);
    for (double value : residuals) {
      hist.Fill(value);
    }
    const int outside = hist.GetBinContent(0) + hist.GetBinContent(bins + 1);
    TCanvas canvas("msEnergyCanvas", "Thermal muon energy residual", 1000, 730);
    canvas.SetMargin(.12, .05, .13, .18);
    canvas.SetTicks(1, 1);
    canvas.SetLogy(logY);
    const double minCount = logY ? .55 : 0.;
    const double maxCount = hist.GetMaximum() * (logY ? 2.6 : 1.15);
    hist.SetMinimum(minCount);
    hist.SetMaximum(maxCount);
    hist.SetLineColor(kAzure + 3);
    hist.SetLineWidth(2);
    hist.SetFillColorAlpha(kAzure + 1, .32);
    hist.GetXaxis()->SetTitleSize(.044);
    hist.GetYaxis()->SetTitleSize(.044);
    hist.GetXaxis()->CenterTitle();
    hist.Draw("HIST");
    TLine zero(0., minCount, 0., maxCount);
    zero.SetLineColor(kGray + 2);
    zero.SetLineStyle(2);
    zero.SetLineWidth(2);
    zero.Draw();
    TLatex text;
    text.SetNDC();
    text.SetTextFont(62);
    text.SetTextSize(.042);
    text.DrawLatex(.12, .948, "Thermal muon energy residual");
    text.SetTextFont(42);
    text.SetTextSize(.030);
    text.DrawLatex(.12, .896, Form("MS-only tracks at generated primary-vertex z  |  %zu tracks  |  %s", residuals.size(), view));
    text.SetTextSize(.025);
    text.DrawLatex(.12, .852, Form("Median %.2f GeV  |  Mean %.2f GeV  |  Bin width %.3g GeV  |  %d tracks outside range",
                                  median, mean, hist.GetBinWidth(1), outside));
    canvas.SaveAs((out + "/" + name + ".png").c_str());
    std::cout << name << " range=[" << low << ", " << high << ") bins=" << bins
              << " outside=" << outside << '\n';
  };

  const double fullLow = std::floor(std::min(0., residuals.front()));
  const double fullHigh = std::floor(std::max(0., residuals.back())) + 1.;
  const double coreLow = std::floor(2. * std::min(0., quantile(residuals, .05))) / 2.;
  const double coreHigh = (std::floor(2. * std::max(0., quantile(residuals, .95))) + 1.) / 2.;
  draw("energyResidual_full", "Full range", 120, fullLow, fullHigh, true);
  draw("energyResidual_core", "Peak region", 64, coreLow, coreHigh, false);

  constexpr double yWidth = .2;
  const auto [minY, maxY] = std::minmax_element(residualsByY.begin(), residualsByY.end());
  const int firstYBin = std::min(5, static_cast<int>(std::floor(minY->first / yWidth)));
  const int endYBin = std::max(20, 1 + static_cast<int>(std::floor(maxY->first / yWidth)));
  const int yBins = endYBin - firstYBin;
  const double yLow = firstYBin * yWidth, yHigh = endYBin * yWidth;
  constexpr double residualLimit = 2.;
  constexpr int residualBins = 30;
  TH2D density("energyResidualVsY", Form(";Generated muon rapidity y_{gen};%s;Tracks / bin", residualTitle),
               yBins, yLow, yHigh, residualBins, -residualLimit, residualLimit);
  density.SetDirectory(nullptr);
  density.SetStats(false);
  std::vector<std::vector<double>> slices(yBins);
  for (const auto& [y, residual] : residualsByY) {
    const int bin = density.GetXaxis()->FindFixBin(y);
    if (bin < 1 || bin > yBins) {
      throw std::runtime_error("Generated rapidity outside histogram");
    }
    density.Fill(y, residual);
    slices[bin - 1].push_back(residual);
  }

  if (density.GetEntries() > 0) {
    const int outside = density.GetEntries() - density.Integral(1, yBins, 1, residualBins);
    TCanvas canvas("msEnergyVsYCanvas", "Thermal muon energy residual versus rapidity", 1050, 750);
    canvas.SetMargin(.12, .17, .13, .16);
    canvas.SetTicks(1, 1);
    canvas.SetLogz();
    density.SetMinimum(1.);
    density.GetXaxis()->CenterTitle();
    density.GetYaxis()->CenterTitle();
    density.GetXaxis()->SetTitleSize(.043);
    density.GetYaxis()->SetTitleSize(.043);
    density.GetYaxis()->SetTitleOffset(1.28);
    density.GetZaxis()->SetTitleOffset(1.25);
    density.Draw("COLZ");
    TLine zero(yLow, 0., yHigh, 0.);
    zero.SetLineStyle(2);
    zero.SetLineColor(kGray + 2);
    zero.Draw();
    TGraphAsymmErrors profile;
    profile.SetMarkerStyle(kFullCircle);
    profile.SetMarkerSize(.9);
    profile.SetLineWidth(2);
    for (int i = 0; i < yBins; ++i) {
      auto& values = slices[i];
      if (values.empty()) {
        continue;
      }
      std::sort(values.begin(), values.end());
      const double q16 = quantile(values, .16);
      const double q50 = quantile(values, .50);
      const double q84 = quantile(values, .84);
      const int point = profile.GetN();
      profile.SetPoint(point, density.GetXaxis()->GetBinCenter(i + 1), q50);
      profile.SetPointError(point, 0., 0., q50 - q16, q84 - q50);
      std::cout << "yGen=[" << density.GetXaxis()->GetBinLowEdge(i + 1) << ", "
                << density.GetXaxis()->GetBinUpEdge(i + 1) << ") tracks=" << values.size()
                << " median=" << q50 << " central68=[" << q16 << ", " << q84 << "]\n";
    }
    profile.Draw("P SAME");
    TLegend legend(.14, .76, .50, .83);
    legend.SetBorderSize(0);
    legend.SetFillStyle(0);
    legend.SetTextSize(.026);
    legend.AddEntry(&profile, "Median and 16th#minus84th percentiles", "pe");
    legend.Draw();
    TLatex text;
    text.SetNDC();
    text.SetTextFont(62);
    text.SetTextSize(.040);
    text.DrawLatex(.12, .945, "Thermal muon energy residual vs. rapidity");
    text.SetTextFont(42);
    text.SetTextSize(.029);
    text.DrawLatex(.12, .895, Form("MS-only at generated primary-vertex z  |  %zu tracks  |  #Deltay = %.1f  |  %d outside view",
                                  residuals.size(), yWidth, outside));
    canvas.SaveAs((out + "/energyResidual_vs_yGen.png").c_str());
    std::cout << "energyResidual_vs_yGen outside displayed residual range=" << outside << '\n';
  }

  std::cout << "events=" << events << " tracks=" << residuals.size() << " invalidEnergy=" << invalidEnergy
            << " failedPropagation=" << failedPropagation << " yRange=[" << yLow << ", " << yHigh << ")"
            << " min=" << residuals.front() << " median=" << median << " mean=" << mean
            << " max=" << residuals.back() << " central68=[" << quantile(residuals, .16) << ", " << quantile(residuals, .84)
            << "] zPV=[" << minZ << ", " << maxZ << "] cm\n";
  for (const auto& [hits, count] : tracksByHits) {
    std::cout << "MS hits=" << hits << " tracks=" << count << '\n';
  }
}
