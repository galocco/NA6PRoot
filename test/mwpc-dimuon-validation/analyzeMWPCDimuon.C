#include <algorithm>
#include <array>
#include <bitset>
#include <cmath>
#include <iostream>
#include <memory>
#include <string>
#include <vector>

#include <TCanvas.h>
#include <TFile.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TLine.h>
#include <TParticle.h>
#include <TString.h>
#include <TStyle.h>
#include <TSystem.h>
#include <TTree.h>

#include "NA6PLayoutParam.h"
#include "NA6PMWPCParam.h"
#include "NA6PMuonSpecModularHit.h"

namespace
{
constexpr int kNStations = 6;

int parentPDG(const std::string& channel)
{
  if (channel == "Jpsi") {
    return 443;
  }
  if (channel == "Omega") {
    return 223;
  }
  if (channel == "Phi") {
    return 333;
  }
  return 0;
}

std::string displayName(const std::string& channel)
{
  if (channel == "Jpsi") {
    return "J/#psi";
  }
  if (channel == "Omega") {
    return "#omega";
  }
  if (channel == "Phi") {
    return "#phi";
  }
  return channel;
}

double momentum(const TParticle& p)
{
  return std::sqrt(p.Px() * p.Px() + p.Py() * p.Py() + p.Pz() * p.Pz());
}

double momentumT(const TParticle& p)
{
  return std::sqrt(p.Px() * p.Px() + p.Py() * p.Py());
}

double rapidity(const TParticle& p)
{
  const double plus = p.Energy() + p.Pz();
  const double minus = p.Energy() - p.Pz();
  if (plus <= 0. || minus <= 0.) {
    return 0.;
  }
  return 0.5 * std::log(plus / minus);
}

int stationFromDetectorID(int detectorID, const NA6PMWPCParam& p, int& localID)
{
  int first = 0;
  for (int station = 0; station < kNStations; ++station) {
    const int n = p.stationGridNX[station] * p.stationGridNY[station];
    if (detectorID >= first && detectorID < first + n) {
      localID = detectorID - first;
      return station;
    }
    first += n;
  }
  localID = -1;
  return -1;
}

int layerFromLocalID(int localID, int nx)
{
  const int row = localID / nx;
  const int col = localID % nx;
  return 2 * (row & 1) + (col & 1); // A,B,C,D -> 0,1,2,3
}

std::unique_ptr<TH1D> makeMomentumHistogram(const char* name, const char* title,
                                            const std::vector<double>& values)
{
  double xmax = 1.;
  if (!values.empty()) {
    xmax = *std::max_element(values.begin(), values.end());
    xmax = std::max(1.e-3, 1.05 * xmax);
  }
  auto h = std::make_unique<TH1D>(name, title, 100, 0., xmax);
  h->SetDirectory(nullptr);
  for (double v : values) {
    h->Fill(v);
  }
  return h;
}

void drawActiveGrid(int station, double xmin, double xmax, double ymin, double ymax,
                    std::vector<std::unique_ptr<TLine>>& lines)
{
  const auto& p = NA6PMWPCParam::Instance();
  const int nx = p.stationGridNX[station];
  const int ny = p.stationGridNY[station];
  const double gasX = (station == 0 && p.useNarrowMS0)
                        ? p.ms0GasX
                        : p.bodyX - 2. * p.innerFrameWidth;
  const double gasY = p.bodyY - 2. * p.innerFrameWidth;
  const double pitchX = gasX - p.activeOverlapX;
  const double pitchY = gasY - p.activeOverlapY;

  auto addVertical = [&](double x) {
    auto line = std::make_unique<TLine>(x, ymin, x, ymax);
    line->SetLineColorAlpha(kGray + 2, 0.35);
    line->SetLineWidth(1);
    line->Draw("same");
    lines.emplace_back(std::move(line));
  };
  auto addHorizontal = [&](double y) {
    auto line = std::make_unique<TLine>(xmin, y, xmax, y);
    line->SetLineColorAlpha(kGray + 2, 0.35);
    line->SetLineWidth(1);
    line->Draw("same");
    lines.emplace_back(std::move(line));
  };

  for (int col = 0; col < nx; ++col) {
    const double xc = (static_cast<double>(col) - 0.5 * (nx - 1)) * pitchX;
    addVertical(xc - 0.5 * gasX);
    addVertical(xc + 0.5 * gasX);
  }
  for (int row = 0; row < ny; ++row) {
    const double yc = (static_cast<double>(row) - 0.5 * (ny - 1)) * pitchY;
    addHorizontal(yc - 0.5 * gasY);
    addHorizontal(yc + 0.5 * gasY);
  }
}
} // namespace

void analyzeMWPCDimuon(const char* runDir = "test_runs/mwpc_dimuon/Jpsi",
                       const char* channelName = "Jpsi",
                       int minBinEntries = 10,
                       double binSizeCm = 10.)
{
  const std::string channel = channelName;
  const int wantedParentPDG = parentPDG(channel);
  if (!wantedParentPDG) {
    std::cerr << "Unknown channel '" << channel
              << "'. Use Jpsi, Omega or Phi." << std::endl;
    return;
  }

  const std::string kineName = std::string(runDir) + "/MCKine.root";
  const std::string hitsName = std::string(runDir) + "/HitsMuonSpecModular.root";

  std::unique_ptr<TFile> kineFile(TFile::Open(kineName.c_str(), "READ"));
  std::unique_ptr<TFile> hitsFile(TFile::Open(hitsName.c_str(), "READ"));
  if (!kineFile || kineFile->IsZombie()) {
    std::cerr << "Cannot open " << kineName << std::endl;
    return;
  }
  if (!hitsFile || hitsFile->IsZombie()) {
    std::cerr << "Cannot open " << hitsName << std::endl;
    return;
  }

  auto* kineTree = dynamic_cast<TTree*>(kineFile->Get("mckine"));
  auto* hitsTree = dynamic_cast<TTree*>(hitsFile->Get("hitsMuonSpecModular"));
  if (!kineTree || !hitsTree) {
    std::cerr << "Required trees mckine / hitsMuonSpecModular were not found" << std::endl;
    return;
  }
  if (kineTree->GetEntries() != hitsTree->GetEntries()) {
    std::cerr << "Entry mismatch: MCKine=" << kineTree->GetEntries()
              << " Hits=" << hitsTree->GetEntries() << std::endl;
    return;
  }

  std::vector<TParticle>* tracks = nullptr;
  std::vector<NA6PMuonSpecModularHit>* hits = nullptr;
  kineTree->SetBranchAddress("tracks", &tracks);
  hitsTree->SetBranchAddress("MuonSpecModular", &hits);

  const auto& p = NA6PMWPCParam::Instance();
  const auto& layout = NA6PLayoutParam::Instance();

  std::array<std::unique_ptr<TH2D>, kNStations> hCount;
  std::array<std::unique_ptr<TH2D>, kNStations> hLayerSum;
  std::array<double, kNStations> xMin{}, xMax{}, yMin{}, yMax{};

  for (int station = 0; station < kNStations; ++station) {
    const int nx = p.stationGridNX[station];
    const int ny = p.stationGridNY[station];
    const double gasX = (station == 0 && p.useNarrowMS0)
                          ? p.ms0GasX
                          : p.bodyX - 2. * p.innerFrameWidth;
    const double gasY = p.bodyY - 2. * p.innerFrameWidth;
    const double pitchX = gasX - p.activeOverlapX;
    const double pitchY = gasY - p.activeOverlapY;
    const double activeX = gasX + (nx - 1) * pitchX;
    const double activeY = gasY + (ny - 1) * pitchY;
    const double margin = 20.;
    xMin[station] = -0.5 * activeX - margin;
    xMax[station] = +0.5 * activeX + margin;
    yMin[station] = -0.5 * activeY - margin;
    yMax[station] = +0.5 * activeY + margin;
    const int nbinX = std::max(10, static_cast<int>(std::ceil((xMax[station] - xMin[station]) / binSizeCm)));
    const int nbinY = std::max(10, static_cast<int>(std::ceil((yMax[station] - yMin[station]) / binSizeCm)));

    hCount[station] = std::make_unique<TH2D>(
      Form("hCount_MS%d", station), "accepted signal muons",
      nbinX, xMin[station], xMax[station], nbinY, yMin[station], yMax[station]);
    hLayerSum[station] = std::make_unique<TH2D>(
      Form("hLayerSum_MS%d", station), "sum of distinct stagger layers",
      nbinX, xMin[station], xMax[station], nbinY, yMin[station], yMax[station]);
    hCount[station]->SetDirectory(nullptr);
    hLayerSum[station]->SetDirectory(nullptr);
  }

  std::vector<double> parentP, parentPt, parentY, daughterP, daughterPt;
  long nGoodPairs = 0;
  long nBadPairs = 0;
  double maxMomentumClosure = 0.;
  double maxEnergyClosure = 0.;

  const Long64_t nEvents = kineTree->GetEntries();
  for (Long64_t iev = 0; iev < nEvents; ++iev) {
    kineTree->GetEntry(iev);
    hitsTree->GetEntry(iev);
    if (!tracks || !hits) {
      continue;
    }

    int parentIndex = -1;
    for (std::size_t i = 0; i < tracks->size(); ++i) {
      if ((*tracks)[i].GetPdgCode() == wantedParentPDG && (*tracks)[i].GetFirstMother() < 0) {
        parentIndex = static_cast<int>(i);
        break;
      }
    }
    if (parentIndex < 0) {
      ++nBadPairs;
      continue;
    }

    const auto& parent = (*tracks)[parentIndex];
    parentP.push_back(momentum(parent));
    parentPt.push_back(momentumT(parent));
    parentY.push_back(rapidity(parent));

    std::vector<int> muonIndices;
    for (std::size_t i = 0; i < tracks->size(); ++i) {
      const auto& tr = (*tracks)[i];
      if (std::abs(tr.GetPdgCode()) == 13 && tr.GetFirstMother() == parentIndex) {
        muonIndices.push_back(static_cast<int>(i));
      }
    }

    if (muonIndices.size() != 2) {
      ++nBadPairs;
      continue;
    }
    ++nGoodPairs;

    double sumPx = 0., sumPy = 0., sumPz = 0., sumE = 0.;
    for (int muonIndex : muonIndices) {
      const auto& mu = (*tracks)[muonIndex];
      daughterP.push_back(momentum(mu));
      daughterPt.push_back(momentumT(mu));
      sumPx += mu.Px();
      sumPy += mu.Py();
      sumPz += mu.Pz();
      sumE += mu.Energy();

      std::array<std::bitset<4>, kNStations> layers;
      std::array<const NA6PMuonSpecModularHit*, kNStations> referenceHit{};
      std::array<double, kNStations> bestDz{};
      bestDz.fill(1.e30);

      for (const auto& hit : *hits) {
        if (hit.getTrackID() != muonIndex) {
          continue;
        }
        int localID = -1;
        const int station = stationFromDetectorID(hit.getDetectorID(), p, localID);
        if (station < 0 || station >= kNStations) {
          continue;
        }
        const int layer = layerFromLocalID(localID, p.stationGridNX[station]);
        layers[station].set(layer);

        const double zNominal = layout.shiftMS[2] + layout.posMSPlaneZ[station];
        const double dz = std::abs(hit.getZIn() - zNominal);
        if (dz < bestDz[station]) {
          bestDz[station] = dz;
          referenceHit[station] = &hit;
        }
      }

      for (int station = 0; station < kNStations; ++station) {
        const auto* hit = referenceHit[station];
        if (!hit || layers[station].none()) {
          continue;
        }
        const double pz = hit->getPZIn();
        if (std::abs(pz) < 1.e-12) {
          continue;
        }
        const double zNominal = layout.shiftMS[2] + layout.posMSPlaneZ[station];
        const double xCentre = layout.shiftMS[0] + layout.posMSPlaneX[station];
        const double yCentre = layout.shiftMS[1] + layout.posMSPlaneY[station];
        const double dz = zNominal - hit->getZIn();
        const double xNominal = hit->getXIn() + dz * hit->getPXIn() / pz - xCentre;
        const double yNominal = hit->getYIn() + dz * hit->getPYIn() / pz - yCentre;
        const double nLayers = static_cast<double>(layers[station].count());
        hCount[station]->Fill(xNominal, yNominal);
        hLayerSum[station]->Fill(xNominal, yNominal, nLayers);
      }
    }

    const double dpx = sumPx - parent.Px();
    const double dpy = sumPy - parent.Py();
    const double dpz = sumPz - parent.Pz();
    maxMomentumClosure = std::max(maxMomentumClosure, std::sqrt(dpx * dpx + dpy * dpy + dpz * dpz));
    maxEnergyClosure = std::max(maxEnergyClosure, std::abs(sumE - parent.Energy()));
  }

  const std::string plotDir = std::string(runDir) + "/plots";
  gSystem->mkdir(plotDir.c_str(), true);
  gStyle->SetOptStat(0);
  gStyle->SetNumberContours(100);

  std::array<std::unique_ptr<TH2D>, kNStations> hMean;
  for (int station = 0; station < kNStations; ++station) {
    hMean[station].reset(static_cast<TH2D*>(hLayerSum[station]->Clone(Form("hMeanLayers_MS%d", station))));
    hMean[station]->Reset("ICES");
    hMean[station]->SetDirectory(nullptr);

    for (int ix = 1; ix <= hMean[station]->GetNbinsX(); ++ix) {
      for (int iy = 1; iy <= hMean[station]->GetNbinsY(); ++iy) {
        const double n = hCount[station]->GetBinContent(ix, iy);
        if (n >= minBinEntries) {
          hMean[station]->SetBinContent(ix, iy, hLayerSum[station]->GetBinContent(ix, iy) / n);
        } else {
          // A value below the displayed range leaves low-statistics bins visually empty.
          hMean[station]->SetBinContent(ix, iy, -1.);
        }
      }
    }

    hMean[station]->SetTitle(Form(
      "MS%d - %s - layer-counted coverage - Ox=Oy=%.1f cm, dz=%.1f cm;"
      "x at nominal station plane [cm];y at nominal station plane [cm];"
      "mean number of A/B/C/D layers hit",
      station, displayName(channel).c_str(), p.activeOverlapX, p.staggerZStep));
    hMean[station]->SetMinimum(0.);
    hMean[station]->SetMaximum(4.);

    TCanvas canvas(Form("cCoverage_MS%d", station), "", 1050, 900);
    canvas.SetRightMargin(0.15);
    hMean[station]->Draw("COLZ");
    std::vector<std::unique_ptr<TLine>> gridLines;
    drawActiveGrid(station, xMin[station], xMax[station], yMin[station], yMax[station], gridLines);
    canvas.SaveAs(Form("%s/coverage_MS%d_%s.png", plotDir.c_str(), station, channel.c_str()));
  }

  auto hParentP = makeMomentumHistogram("hParentP", Form("%s parent |p|;|p| [GeV/c];parents", displayName(channel).c_str()), parentP);
  auto hParentPt = makeMomentumHistogram("hParentPt", Form("%s parent p_{T};p_{T} [GeV/c];parents", displayName(channel).c_str()), parentPt);
  auto hDaughterP = makeMomentumHistogram("hDaughterP", Form("%s daughter muons |p|;|p| [GeV/c];muons", displayName(channel).c_str()), daughterP);
  auto hDaughterPt = makeMomentumHistogram("hDaughterPt", Form("%s daughter muons p_{T};p_{T} [GeV/c];muons", displayName(channel).c_str()), daughterPt);
  auto hParentY = std::make_unique<TH1D>("hParentY", Form("%s parent rapidity;rapidity y;parents", displayName(channel).c_str()), 120, 0., 6.);
  hParentY->SetDirectory(nullptr);
  for (double y : parentY) {
    hParentY->Fill(y);
  }

  TCanvas cKin("cKinematics", "", 1500, 900);
  cKin.Divide(3, 2);
  cKin.cd(1);
  gPad->SetLogy();
  hParentP->Draw("HIST");
  cKin.cd(2);
  hParentPt->Draw("HIST");
  cKin.cd(3);
  hParentY->Draw("HIST");
  cKin.cd(4);
  gPad->SetLogy();
  hDaughterP->Draw("HIST");
  cKin.cd(5);
  hDaughterPt->Draw("HIST");
  cKin.SaveAs(Form("%s/kinematics_%s.png", plotDir.c_str(), channel.c_str()));

  TCanvas cRap("cParentRapidity", "", 900, 700);
  hParentY->Draw("HIST");
  cRap.SaveAs(Form("%s/parent_rapidity_%s.png", plotDir.c_str(), channel.c_str()));

  TFile out(Form("%s/validation_%s.root", plotDir.c_str(), channel.c_str()), "RECREATE");
  hParentP->Write();
  hParentPt->Write();
  hParentY->Write();
  hDaughterP->Write();
  hDaughterPt->Write();
  for (int station = 0; station < kNStations; ++station) {
    hCount[station]->Write();
    hLayerSum[station]->Write();
    hMean[station]->Write();
  }
  out.Close();

  std::cout << "\n========== MWPC DIMUON VALIDATION: " << channel << " ==========\n";
  std::cout << "events                 = " << nEvents << '\n';
  std::cout << "good dimuon pairs      = " << nGoodPairs << '\n';
  std::cout << "events without 2 muons = " << nBadPairs << '\n';
  std::cout << "saved daughter muons   = " << daughterP.size() << '\n';
  std::cout << "max |sum p_mu - p_M|   = " << maxMomentumClosure << " GeV/c\n";
  std::cout << "max |sum E_mu - E_M|   = " << maxEnergyClosure << " GeV\n";
  for (int station = 0; station < kNStations; ++station) {
    std::cout << "MS" << station << " accepted daughter muons = "
              << static_cast<long>(hCount[station]->GetEntries()) << '\n';
  }
  std::cout << "plots: " << plotDir << "\n";
  std::cout << "========================================================\n";
}
