#include <algorithm>
#include <cmath>
#include <iostream>
#include <memory>
#include <string>
#include <vector>

#include <TCanvas.h>
#include <TFile.h>
#include <TH1D.h>
#include <TParticle.h>
#include <TString.h>
#include <TStyle.h>
#include <TSystem.h>
#include <TTree.h>

namespace
{
int parentPDG(const std::string& channel)
{
  if (channel == "Jpsi") return 443;
  if (channel == "Omega") return 223;
  if (channel == "Phi") return 333;
  return 0;
}

std::string displayName(const std::string& channel)
{
  if (channel == "Jpsi") return "J/#psi";
  if (channel == "Omega") return "#omega";
  if (channel == "Phi") return "#phi";
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
  if (plus <= 0. || minus <= 0.) return 0.;
  return 0.5 * std::log(plus / minus);
}

std::unique_ptr<TH1D> makeHistogram(const char* name, const char* title,
                                    const std::vector<double>& values,
                                    int bins = 100)
{
  double xmax = 1.;
  if (!values.empty()) {
    xmax = std::max(1.e-3, 1.05 * *std::max_element(values.begin(), values.end()));
  }
  auto h = std::make_unique<TH1D>(name, title, bins, 0., xmax);
  h->SetDirectory(nullptr);
  for (double v : values) h->Fill(v);
  return h;
}
} // namespace

void analyzeMWPCDimuon(const char* runDir = "test_runs/mwpc_dimuon/Jpsi",
                       const char* channelName = "Jpsi")
{
  const std::string channel = channelName;
  const int wantedParentPDG = parentPDG(channel);
  if (!wantedParentPDG) {
    std::cerr << "Unknown channel '" << channel << "'. Use Jpsi, Omega or Phi.\n";
    return;
  }

  const std::string kineName = std::string(runDir) + "/MCKine.root";
  std::unique_ptr<TFile> kineFile(TFile::Open(kineName.c_str(), "READ"));
  if (!kineFile || kineFile->IsZombie()) {
    std::cerr << "Cannot open " << kineName << '\n';
    return;
  }

  auto* kineTree = dynamic_cast<TTree*>(kineFile->Get("mckine"));
  if (!kineTree) {
    std::cerr << "Tree mckine was not found in " << kineName << '\n';
    return;
  }

  std::vector<TParticle>* tracks = nullptr;
  kineTree->SetBranchAddress("tracks", &tracks);

  std::vector<double> parentP, parentPt, parentY;
  std::vector<double> daughterP, daughterPt;
  long nGoodPairs = 0;
  long nBadPairs = 0;
  double maxMomentumClosure = 0.;
  double maxEnergyClosure = 0.;

  const Long64_t nEvents = kineTree->GetEntries();
  for (Long64_t iev = 0; iev < nEvents; ++iev) {
    kineTree->GetEntry(iev);
    if (!tracks) continue;

    int parentIndex = -1;
    for (std::size_t i = 0; i < tracks->size(); ++i) {
      const auto& tr = (*tracks)[i];
      if (tr.GetPdgCode() == wantedParentPDG && tr.GetFirstMother() < 0) {
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
    }

    const double dpx = sumPx - parent.Px();
    const double dpy = sumPy - parent.Py();
    const double dpz = sumPz - parent.Pz();
    maxMomentumClosure = std::max(maxMomentumClosure,
                                  std::sqrt(dpx * dpx + dpy * dpy + dpz * dpz));
    maxEnergyClosure = std::max(maxEnergyClosure, std::abs(sumE - parent.Energy()));
  }

  auto hParentP = makeHistogram("hParentP",
    Form("%s parent |p|;|p| [GeV/c];parents", displayName(channel).c_str()), parentP);
  auto hParentPt = makeHistogram("hParentPt",
    Form("%s parent p_{T};p_{T} [GeV/c];parents", displayName(channel).c_str()), parentPt);
  auto hDaughterP = makeHistogram("hDaughterP",
    Form("%s daughter muons |p|;|p| [GeV/c];muons", displayName(channel).c_str()), daughterP);
  auto hDaughterPt = makeHistogram("hDaughterPt",
    Form("%s daughter muons p_{T};p_{T} [GeV/c];muons", displayName(channel).c_str()), daughterPt);
  auto hParentY = std::make_unique<TH1D>(
    "hParentY", Form("%s parent rapidity;rapidity y;parents", displayName(channel).c_str()),
    120, 0., 6.);
  hParentY->SetDirectory(nullptr);
  for (double y : parentY) hParentY->Fill(y);

  const std::string plotDir = std::string(runDir) + "/plots";
  gSystem->mkdir(plotDir.c_str(), true);
  gStyle->SetOptStat(0);

  TCanvas cKin("cKinematics", "", 1500, 900);
  cKin.Divide(3, 2);
  cKin.cd(1); gPad->SetLogy(); hParentP->Draw("HIST");
  cKin.cd(2); hParentPt->Draw("HIST");
  cKin.cd(3); hParentY->Draw("HIST");
  cKin.cd(4); gPad->SetLogy(); hDaughterP->Draw("HIST");
  cKin.cd(5); hDaughterPt->Draw("HIST");
  cKin.SaveAs(Form("%s/kinematics_%s.png", plotDir.c_str(), channel.c_str()));

  TFile out(Form("%s/kinematics_%s.root", plotDir.c_str(), channel.c_str()), "RECREATE");
  hParentP->Write();
  hParentPt->Write();
  hParentY->Write();
  hDaughterP->Write();
  hDaughterPt->Write();
  out.Close();

  std::cout << "\n========== DIMUON KINEMATICS: " << channel << " ==========\n";
  std::cout << "events                 = " << nEvents << '\n';
  std::cout << "good dimuon pairs      = " << nGoodPairs << '\n';
  std::cout << "events without 2 muons = " << nBadPairs << '\n';
  std::cout << "saved daughter muons   = " << daughterP.size() << '\n';
  std::cout << "max |sum p_mu - p_M|   = " << maxMomentumClosure << " GeV/c\n";
  std::cout << "max |sum E_mu - E_M|   = " << maxEnergyClosure << " GeV\n";
  std::cout << "plots: " << plotDir << '\n';
  std::cout << "================================================\n";
}
