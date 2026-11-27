// Compare paired MS-only energy estimates at the generated primary vertex.
#include "NA6PLayoutParam.h"
#include "NA6PMCComposedLabel.h"
#include "NA6PMCEventHeader.h"
#include "NA6PRecoParam.h"
#include "NA6PTrack.h"
#include "Propagator.h"
#include <TCanvas.h>
#include <TFile.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TLegend.h>
#include <TLine.h>
#include <TMultiGraph.h>
#include <TParticle.h>
#include <TSystem.h>
#include <TTreeReader.h>
#include <TTreeReaderValue.h>
#include <algorithm>
#include <array>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <numeric>
#include <stdexcept>
#include <string>
#include <vector>

R__LOAD_LIBRARY(librecLib)

namespace
{
struct EnergyResiduals {
  std::vector<double> absolute;
  std::vector<double> relative;
};

struct RapidityBin {
  std::array<EnergyResiduals, 2> residuals;
  double sumRapidity = 0.;
  double sumEnergy = 0.;
};

struct Resolution {
  double mean, median, rms, sigma68;
};

Resolution summarize(std::vector<double>& values)
{
  std::sort(values.begin(), values.end());
  const auto quantile = [&](double fraction) {
    const double index = fraction * (values.size() - 1);
    const size_t low = index;
    return std::lerp(values[low], values[std::min(low + 1, values.size() - 1)], index - low);
  };
  const double mean = std::accumulate(values.begin(), values.end(), 0.) / values.size();
  double variance = 0.;
  for (const double value : values) {
    variance += (value - mean) * (value - mean);
  }
  return {mean, quantile(.50), std::sqrt(variance / values.size()), .5 * (quantile(.84) - quantile(.16))};
}
}

void compareMSMuonEnergy(const char* runDir, Long64_t maxEvents = -1, const char* outputDir = "", double rapidityWidth = .2, int thermalSelection = -1)
{
  if (!std::isfinite(rapidityWidth) || rapidityWidth <= 0.) {
    throw std::invalid_argument("Rapidity bin width must be positive and finite");
  }
  const std::string run = runDir;
  const std::string out = *outputDir ? outputDir : run + "/ms_muon_energy_resolution";
  const std::string rapidityOut = out + "/bias_vs_rapidity" + (rapidityWidth == .2 ? "" : Form("_dy%g", rapidityWidth));
  TFile reco((run + "/reco/TracksMuonSpec.root").c_str());
  TFile kine((run + "/sim/MCKine.root").c_str());
  auto* msTree = reco.Get<TTree>("tracksMuonSpec");
  auto* mcTree = kine.Get<TTree>("mckine");
  if (!msTree || !mcTree || msTree->GetEntries() != mcTree->GetEntries()) {
    throw std::runtime_error("Missing or inconsistent MS and MC trees");
  }

  na6p::conf::ConfigurableParam::updateFromFile(run + "/reco/na6pLayout.ini", "layout");
  na6p::conf::ConfigurableParam::updateFromFile(run + "/reco/na6pRecoParam.ini", "reco");
  if (!Propagator::loadField() || !Propagator::loadGeometry(run + "/sim/geometry.root")) {
    throw std::runtime_error("Cannot initialize field or geometry");
  }
  std::array<Propagator::PropOpt, 2> options;
  for (auto& opt : options) {
    opt.maxStep = NA6PRecoParam::Instance().maxPropagationStep;
    if (!std::isfinite(opt.maxStep) || opt.maxStep <= 0.f) {
      throw std::runtime_error("Invalid maximum propagation step");
    }
  }
  options[1].muonRadiativeLoss = true;
  const auto* propagator = Propagator::Instance();
  TTreeReader ms(msTree), mc(mcTree);
  TTreeReaderValue<std::vector<NA6PTrack>> tracks(ms, "MuonSpec");
  TTreeReaderValue<std::vector<NA6PMCComposedLabel>> labels(ms, "MuonSpecMCTruth");
  TTreeReaderValue<NA6PMCEventHeader> header(mc, "header");
  TTreeReaderValue<std::vector<TParticle>> particles(mc, "tracks");
  constexpr double energyEdges[] = {0., 2., 4., 6., 10., 15., 20., 30., 40., 60., std::numeric_limits<double>::infinity()};
  constexpr size_t nBins = std::size(energyEdges) - 1;
  constexpr size_t overall = nBins;
  constexpr double vertexTolerance = 1.e-3; // cm, compare production point with the MC primary vertex
  std::array<std::array<EnergyResiduals, 2>, nBins + 1> residuals;
  std::map<int, RapidityBin> residualsByRapidity;
  std::array<double, nBins> sumEnergy{};
  std::array<Long64_t, 2> failed{};
  const double mass = PID::getMass(PID::Muon);
  const Long64_t nEvents = maxEvents < 0 ? mcTree->GetEntries() : std::min(maxEvents, mcTree->GetEntries());
  Long64_t selected = 0;
  double sumEnergyShift = 0.;
  for (Long64_t event = 0; event < nEvents; ++event) {
    if (!ms.Next() || !mc.Next() || tracks->size() != labels->size()) {
      throw std::runtime_error("Inconsistent MS and MC data at entry " + std::to_string(event));
    }
    const double xPV = header->getVX(), yPV = header->getVY(), zPV = header->getVZ();
    if (!std::isfinite(xPV) || !std::isfinite(yPV) || !std::isfinite(zPV)) {
      throw std::runtime_error("Invalid primary vertex");
    }
    for (size_t i = 0; i < tracks->size(); ++i) {
      const auto& label = labels->at(i);
      if (!label.isValid() || label.isFake() || label.getSourceID() != 0) {
        continue;
      }
      const int id = label.getTrackID();
      // Event IDs in composed labels are stored modulo 2^nbitsEvID.
      if (label.getEventID() != (header->getEventID() & NA6PMCComposedLabel::maskEvID) ||
          id < 0 || static_cast<size_t>(id) >= particles->size()) {
        throw std::runtime_error("MS label does not identify a particle in the corresponding MC event");
      }
      const auto& particle = particles->at(id);
      const auto& track = tracks->at(i);
      const double eGen = particle.Energy();
      if (std::abs(particle.GetPdgCode()) != 13 ||
          (thermalSelection != 0 && !(std::hypot(particle.Vx() - xPV, particle.Vy() - yPV, particle.Vz() - zPV) <= vertexTolerance)) ||
          !std::isfinite(eGen) || eGen < mass || !track.isValid() || !(track.getZ() > zPV) ||
          !std::isfinite(track.getP()) || track.getP() <= 0.) {
        continue;
      }
      if (thermalSelection >= 0) {
        const auto& generators = header->getGenHeaders();
        const bool thermal = std::any_of(generators.begin(), generators.end(), [id](const auto& generator) {
          const auto& name = generator.getInfo();
          const auto begin = generator.getPrimariesOffset();
          return (name.starts_with("VirtualPhoton_pdg23") || name.starts_with("thermalMuPlus_pdg-13") ||
                  name.starts_with("thermalMuMinus_pdg13")) &&
                 static_cast<unsigned>(id) >= begin && static_cast<unsigned>(id) - begin < generator.getNPrimaries();
        });
        if (thermal != (thermalSelection == 1)) {
          continue;
        }
      }
      ++selected;
      std::array<double, 2> eRec{};
      bool paired = true;
      for (size_t mode = 0; mode < options.size(); ++mode) {
        NA6PTrackParCov atPV = track;
        atPV.setPID(PID::Muon);
        if (!propagator->propagateToZ(atPV, zPV, options[mode]) || !std::isfinite(atPV.getP()) || atPV.getP() <= 0.f) {
          ++failed[mode];
          paired = false;
          continue;
        }
        eRec[mode] = std::hypot(atPV.getP(), mass);
      }
      if (!paired) {
        continue;
      }
      const size_t bin = std::upper_bound(std::begin(energyEdges), std::end(energyEdges), eGen) - std::begin(energyEdges) - 1;
      sumEnergy[bin] += eGen;
      sumEnergyShift += eRec[1] - eRec[0];
      const double yGen = std::asinh(particle.Pz() / std::hypot(particle.Pt(), mass));
      if (!std::isfinite(yGen)) {
        throw std::runtime_error("Invalid generated muon rapidity");
      }
      const double yBin = std::floor(yGen / rapidityWidth);
      if (yBin < std::numeric_limits<int>::min() || yBin >= std::numeric_limits<int>::max()) {
        throw std::runtime_error("Rapidity bin index out of range");
      }
      auto& rapidityBin = residualsByRapidity[static_cast<int>(yBin)];
      rapidityBin.sumRapidity += yGen;
      rapidityBin.sumEnergy += eGen;
      for (size_t mode = 0; mode < options.size(); ++mode) {
        const double delta = eGen - eRec[mode];
        const double relative = 100. * delta / eGen;
        for (const size_t slot : {bin, overall}) {
          residuals[slot][mode].absolute.push_back(delta);
          residuals[slot][mode].relative.push_back(relative);
        }
        rapidityBin.residuals[mode].absolute.push_back(delta);
        rapidityBin.residuals[mode].relative.push_back(relative);
      }
    }
    if ((event + 1) % 10000 == 0) {
      std::cout << "Events processed: " << event + 1 << '/' << nEvents << std::endl;
    }
  }

  const auto paired = residuals[overall][0].absolute.size();
  std::cout << "Events=" << nEvents << " selected=" << selected << " paired=" << paired
            << " failed_ionization=" << failed[0] << " failed_radiative=" << failed[1] << '\n';
  if (paired == 0) {
    throw std::runtime_error("No paired MS muon energies at the primary vertex");
  }
  std::cout << "Mean Erec shift (radiative - ionization) = " << sumEnergyShift / paired << " GeV\n";
  if (gSystem->mkdir(out.c_str(), true) != 0 && gSystem->AccessPathName(out.c_str(), kWritePermission)) {
    throw std::runtime_error("Cannot create resolution output directory");
  }

  std::ofstream csv(out + "/resolution.csv");
  csv << std::setprecision(9) << "sample,e_gen_low_GeV,e_gen_high_GeV,model,count,mean_GeV,median_GeV,rms_GeV,sigma68_GeV,"
         "mean_pct,median_pct,rms_pct,sigma68_pct\n";
  TCanvas canvas("msResolution", "MS energy resolution", 1100, 520);
  canvas.Divide(2, 1);
  TMultiGraph bias, width;
  bias.SetTitle(";Mean generated energy in bin (GeV);Median of (E_{gen} - E_{rec}) / E_{gen} (%)");
  width.SetTitle(";Mean generated energy in bin (GeV);Central 68% half-width (%)");
  TLegend legend(.55, .7, .88, .88);
  constexpr const char* modes[] = {"Ionization", "Ionization + radiative"};
  for (size_t mode = 0; mode < options.size(); ++mode) {
    auto* biasGraph = new TGraph;
    auto* widthGraph = new TGraph;
    for (auto* graph : {biasGraph, widthGraph}) {
      graph->SetMarkerStyle(mode == 0 ? kOpenCircle : kFullSquare);
      graph->SetMarkerColor(mode == 0 ? kAzure + 2 : kRed + 1);
    }
    for (size_t bin = 0; bin <= overall; ++bin) {
      auto& sample = residuals[bin][mode];
      const auto count = sample.absolute.size();
      if (count == 0) {
        continue;
      }
      const auto abs = summarize(sample.absolute);
      const auto rel = summarize(sample.relative);
      if (bin == overall) {
        csv << "all,,,";
        std::cout << modes[mode] << "  mean=" << abs.mean << " GeV  median=" << abs.median
                  << " GeV  RMS=" << abs.rms << " GeV  sigma68=" << abs.sigma68 << " GeV (" << rel.sigma68 << "%)\n";
      } else {
        csv << "bin," << energyEdges[bin] << ',' << energyEdges[bin + 1] << ',';
        if (count >= 10) {
          const double x = sumEnergy[bin] / count;
          biasGraph->AddPoint(x, rel.median);
          widthGraph->AddPoint(x, rel.sigma68);
        }
      }
      csv << modes[mode] << ',' << count << ',' << abs.mean << ',' << abs.median << ',' << abs.rms << ',' << abs.sigma68
          << ',' << rel.mean << ',' << rel.median << ',' << rel.rms << ',' << rel.sigma68 << '\n';
    }
    bias.Add(biasGraph);
    width.Add(widthGraph);
    legend.AddEntry(widthGraph, modes[mode], "p");
  }
  csv.close();
  if (!csv) {
    throw std::runtime_error("Cannot write resolution CSV");
  }
  if (bias.GetListOfGraphs()->First() && static_cast<TGraph*>(bias.GetListOfGraphs()->First())->GetN() > 0) {
    canvas.cd(1);
    gPad->SetGrid();
    bias.Draw("AP");
    canvas.cd(2);
    gPad->SetGrid();
    width.Draw("AP");
    legend.Draw();
    canvas.SaveAs((out + "/resolution_vs_energy.png").c_str());
  }

  std::ofstream rapidityCSV(rapidityOut + ".csv");
  rapidityCSV << std::setprecision(9)
              << "y_low,y_high,y_mean,e_gen_mean_GeV,model,count,mean_MeV,median_MeV,mean_sem_MeV,mean_pct,median_pct,mean_sem_pct\n";
  TCanvas rapidityCanvas("msRapidityBias", "MS energy bias versus generated rapidity", 1120, 520);
  rapidityCanvas.Divide(2, 1);
  TMultiGraph absoluteBias, relativeBias;
  absoluteBias.SetTitle(";Generated muon rapidity y_{gen};Mean E_{gen} - E_{rec} (MeV)");
  relativeBias.SetTitle(";Generated muon rapidity y_{gen};Mean (E_{gen} - E_{rec}) / E_{gen} (%)");
  TLegend rapidityLegend(.20, .15, .57, .28);
  bool hasPlottedBins = false;
  for (size_t mode = 0; mode < options.size(); ++mode) {
    auto* absGraph = new TGraphErrors;
    auto* relGraph = new TGraphErrors;
    for (auto* graph : {absGraph, relGraph}) {
      graph->SetMarkerStyle(mode == 0 ? kOpenCircle : kFullSquare);
      graph->SetMarkerColor(mode == 0 ? kAzure + 2 : kRed + 1);
      graph->SetLineColor(graph->GetMarkerColor());
    }
    for (auto& [bin, values] : residualsByRapidity) {
      auto& sample = values.residuals[mode];
      const auto count = sample.absolute.size();
      const double x = values.sumRapidity / count;
      const auto abs = summarize(sample.absolute);
      const auto rel = summarize(sample.relative);
      // The RMS uses 1/N, hence the standard error of the mean is RMS/sqrt(N-1).
      const double semScale = count > 1 ? 1. / std::sqrt(count - 1.) : std::numeric_limits<double>::quiet_NaN();
      const double absSem = abs.rms * semScale * 1000.;
      const double relSem = rel.rms * semScale;
      rapidityCSV << bin * rapidityWidth << ',' << (bin + 1) * rapidityWidth << ',' << x << ','
                  << values.sumEnergy / count << ',' << modes[mode] << ',' << count << ','
                  << abs.mean * 1000. << ',' << abs.median * 1000. << ',' << absSem << ','
                  << rel.mean << ',' << rel.median << ',' << relSem << '\n';
      if (count >= 10) {
        const int point = absGraph->GetN();
        absGraph->SetPoint(point, x, abs.mean * 1000.);
        absGraph->SetPointError(point, 0., absSem);
        relGraph->SetPoint(point, x, rel.mean);
        relGraph->SetPointError(point, 0., relSem);
        hasPlottedBins = true;
      }
    }
    absoluteBias.Add(absGraph);
    relativeBias.Add(relGraph);
    rapidityLegend.AddEntry(relGraph, modes[mode], "pe");
  }
  rapidityCSV.close();
  if (!rapidityCSV) {
    throw std::runtime_error("Cannot write rapidity bias CSV");
  }
  if (hasPlottedBins) {
    std::array<TMultiGraph*, 2> graphs{&absoluteBias, &relativeBias};
    for (size_t i = 0; i < graphs.size(); ++i) {
      rapidityCanvas.cd(i + 1);
      gPad->SetGrid();
      gPad->SetLeftMargin(.16);
      graphs[i]->Draw("AP");
      auto* axis = graphs[i]->GetXaxis();
      TLine zero;
      zero.SetLineStyle(kDashed);
      zero.DrawLine(axis->GetXmin(), 0., axis->GetXmax(), 0.);
    }
    rapidityLegend.Draw();
    rapidityCanvas.SaveAs((rapidityOut + ".png").c_str());
  }
}
