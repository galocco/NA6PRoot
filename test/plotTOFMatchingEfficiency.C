#if !defined(__CINT__) || defined(__MAKECINT__)
#include <TTree.h>
#include <TFile.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TEfficiency.h>
#include <TParticle.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <array>
#include <cmath>
#include <string>
#include <vector>
#include "NA6PTrack.h"
#include "NA6PTOFCluster.h"
#include "NA6PLayoutParam.h"
#endif

namespace TOFMatchEffConst {
constexpr int minVerTelHits = 4;
constexpr double minMomentum = 0.;
constexpr double maxMomentum = 50.;
constexpr double minPt = 0.;
constexpr double maxPt = 5.;
constexpr double minEta = -3.;
constexpr double maxEta = 7.;
}

struct SpeciesConfig {
  const char* key;
  const char* label;
  int absPdg;
};

struct SpeciesHistograms {
  TH1F* allVsP = nullptr;
  TH1F* matchedVsP = nullptr;
  TH1F* trueMatchedVsP = nullptr;
  TH1F* allVsPt = nullptr;
  TH1F* matchedVsPt = nullptr;
  TH1F* trueMatchedVsPt = nullptr;
  TH1F* allVsEta = nullptr;
  TH1F* matchedVsEta = nullptr;
  TH1F* trueMatchedVsEta = nullptr;
  TH2F* allPtVsEta = nullptr;
  TH2F* matchedPtVsEta = nullptr;
  TH2F* trueMatchedPtVsEta = nullptr;
};

TTree* openFileAndTree(const char* path, const char* treeName)
{
  TFile* file = TFile::Open(path);
  if (!file || file->IsZombie()) {
    printf("ERROR: Cannot open %s\n", path);
    return nullptr;
  }
  TTree* tree = dynamic_cast<TTree*>(file->Get(treeName));
  if (!tree) {
    printf("ERROR: Cannot find tree '%s' in %s\n", treeName, path);
    return nullptr;
  }
  return tree;
}

double getTrackPt(const NA6PTrack& track)
{
  double pxyz[3];
  track.getPXYZ(pxyz);
  return std::hypot(pxyz[0], pxyz[1]);
}

double getTrackEta(const NA6PTrack& track)
{
  double pxyz[3];
  track.getPXYZ(pxyz);
  const double pt = std::hypot(pxyz[0], pxyz[1]);
  if (pt < 1.e-12) {
    return pxyz[2] >= 0. ? TOFMatchEffConst::maxEta : TOFMatchEffConst::minEta;
  }
  return std::asinh(pxyz[2] / pt);
}

TEfficiency* makeEfficiency1D(TH1* numerator, TH1* denominator, const char* name, const char* title)
{
  for (int bin = 1; bin <= numerator->GetNbinsX(); ++bin) {
    if (numerator->GetBinContent(bin) > denominator->GetBinContent(bin)) {
      numerator->SetBinContent(bin, denominator->GetBinContent(bin));
    }
  }
  if (!TEfficiency::CheckConsistency(*numerator, *denominator)) {
    printf("WARNING: TEfficiency consistency failed for %s\n", name);
    return nullptr;
  }
  auto* efficiency = new TEfficiency(*numerator, *denominator);
  efficiency->SetName(name);
  efficiency->SetTitle(title);
  return efficiency;
}

TEfficiency* makeEfficiency2D(TH2* numerator, TH2* denominator, const char* name, const char* title)
{
  for (int binX = 1; binX <= numerator->GetNbinsX(); ++binX) {
    for (int binY = 1; binY <= numerator->GetNbinsY(); ++binY) {
      if (numerator->GetBinContent(binX, binY) > denominator->GetBinContent(binX, binY)) {
        numerator->SetBinContent(binX, binY, denominator->GetBinContent(binX, binY));
      }
    }
  }
  if (!TEfficiency::CheckConsistency(*numerator, *denominator)) {
    printf("WARNING: 2D TEfficiency consistency failed for %s\n", name);
    return nullptr;
  }
  auto* efficiency = new TEfficiency(*numerator, *denominator);
  efficiency->SetName(name);
  efficiency->SetTitle(title);
  return efficiency;
}

void drawEfficiency1D(TEfficiency* efficiency, const char* canvasName, const char* canvasTitle, const char* outFile)
{
  if (!efficiency) {
    return;
  }
  auto* canvas = new TCanvas(canvasName, canvasTitle, 800, 600);
  canvas->SetLeftMargin(0.14);
  canvas->SetBottomMargin(0.13);
  canvas->SetRightMargin(0.05);
  canvas->SetTopMargin(0.07);
  efficiency->SetMarkerStyle(20);
  efficiency->SetMarkerSize(0.8);
  efficiency->SetLineWidth(2);
  efficiency->Draw("AP");
  gPad->Modified();
  gPad->Update();
  if (auto* graph = efficiency->GetPaintedGraph()) {
    graph->GetYaxis()->SetRangeUser(0., 1.05);
    graph->GetXaxis()->SetTitleSize(0.05);
    graph->GetYaxis()->SetTitleSize(0.05);
  }
  gPad->Modified();
  gPad->Update();
  canvas->SaveAs(outFile);
}

void drawEfficiency2D(TEfficiency* efficiency, const char* canvasName, const char* canvasTitle, const char* outFile)
{
  if (!efficiency) {
    return;
  }
  auto* canvas = new TCanvas(canvasName, canvasTitle, 1000, 800);
  canvas->SetLeftMargin(0.14);
  canvas->SetBottomMargin(0.12);
  canvas->SetRightMargin(0.15);
  canvas->SetTopMargin(0.07);
  efficiency->Draw("COLZ");
  canvas->SaveAs(outFile);
}

void plotTOFMatchingEfficiency(const char* dirSimu = ".")
{
  gStyle->SetOptStat(0);
  gStyle->SetTitleSize(0.05, "XY");
  gStyle->SetLabelSize(0.045, "XY");
  gStyle->SetTitleOffset(1.1, "X");
  gStyle->SetTitleOffset(1.0, "Y");

  auto* trackTree = openFileAndTree(Form("%s/TracksTOFMatching.root", dirSimu), "tracksTOFMatching");
  auto* mcTree = openFileAndTree(Form("%s/MCKine.root", dirSimu), "mckine");
  auto* tClus   = openFileAndTree(Form("%s/ClustersTOF.root",        dirSimu), "clustersTOF");
  if (!trackTree || !mcTree || !tClus) {
    return;
  }

  std::vector<NA6PTrack> tracks;
  std::vector<NA6PTrack>* trackPtr = &tracks;
  std::vector<TParticle>* mcArr = nullptr;
  std::vector<NA6PTOFCluster>  tofClus,     *tofClusPtr   = &tofClus;

  trackTree->SetBranchAddress("TOFMatching", &trackPtr);
  mcTree->SetBranchAddress("tracks", &mcArr);
  tClus->SetBranchAddress("TOF", &tofClusPtr);

  const int nEvents = std::min({(int)trackTree->GetEntries(), (int)mcTree->GetEntries(), (int)tClus->GetEntries()});
  printf("Processing %d events for TOF matching efficiency\n", nEvents);
  const int nVTLayers = NA6PLayoutParam::Instance().nVerTelPlanes;

  const int nMomBins = 100;
  const int nPtBins = 100;
  const int nEtaBins = 100;

  auto* hAllVsP = new TH1F("hAllVsP", ";p (GeV/c);tracks", nMomBins,
                           TOFMatchEffConst::minMomentum, TOFMatchEffConst::maxMomentum);
  auto* hMatchedVsP = new TH1F("hMatchedVsP", ";p (GeV/c);tracks", nMomBins,
                               TOFMatchEffConst::minMomentum, TOFMatchEffConst::maxMomentum);
  auto* hTrueMatchedVsP = new TH1F("hTrueMatchedVsP", ";p (GeV/c);tracks", nMomBins,
                                   TOFMatchEffConst::minMomentum, TOFMatchEffConst::maxMomentum);

  auto* hAllVsPt = new TH1F("hAllVsPt", ";p_{T} (GeV/c);tracks", nPtBins,
                            TOFMatchEffConst::minPt, TOFMatchEffConst::maxPt);
  auto* hMatchedVsPt = new TH1F("hMatchedVsPt", ";p_{T} (GeV/c);tracks", nPtBins,
                                TOFMatchEffConst::minPt, TOFMatchEffConst::maxPt);
  auto* hTrueMatchedVsPt = new TH1F("hTrueMatchedVsPt", ";p_{T} (GeV/c);tracks", nPtBins,
                                    TOFMatchEffConst::minPt, TOFMatchEffConst::maxPt);

  auto* hAllVsEta = new TH1F("hAllVsEta", ";#eta;tracks", nEtaBins,
                             TOFMatchEffConst::minEta, TOFMatchEffConst::maxEta);
  auto* hMatchedVsEta = new TH1F("hMatchedVsEta", ";#eta;tracks", nEtaBins,
                                 TOFMatchEffConst::minEta, TOFMatchEffConst::maxEta);
  auto* hTrueMatchedVsEta = new TH1F("hTrueMatchedVsEta", ";#eta;tracks", nEtaBins,
                                     TOFMatchEffConst::minEta, TOFMatchEffConst::maxEta);

  auto* hAllPtVsEta = new TH2F("hAllPtVsEta", ";p_{T} (GeV/c);#eta;tracks", nPtBins,
                               TOFMatchEffConst::minPt, TOFMatchEffConst::maxPt,
                               nEtaBins, TOFMatchEffConst::minEta, TOFMatchEffConst::maxEta);
  auto* hMatchedPtVsEta = new TH2F("hMatchedPtVsEta", ";p_{T} (GeV/c);#eta;tracks", nPtBins,
                                   TOFMatchEffConst::minPt, TOFMatchEffConst::maxPt,
                                   nEtaBins, TOFMatchEffConst::minEta, TOFMatchEffConst::maxEta);
  auto* hTrueMatchedPtVsEta = new TH2F("hTrueMatchedPtVsEta", ";p_{T} (GeV/c);#eta;tracks", nPtBins,
                                       TOFMatchEffConst::minPt, TOFMatchEffConst::maxPt,
                                       nEtaBins, TOFMatchEffConst::minEta, TOFMatchEffConst::maxEta);

  const std::array<SpeciesConfig, 3> speciesConfigs{{
    {"p", "p", 2212},
    {"K", "K", 321},
    {"pi", "#pi", 211},
  }};
  std::array<SpeciesHistograms, 3> speciesHistograms;

  auto bookSpeciesHistograms = [&](const SpeciesConfig& species) {
    SpeciesHistograms histograms;
    histograms.allVsP = new TH1F(Form("hAllVsP_%s", species.key),
                                 Form(";%s momentum p (GeV/c);tracks", species.label),
                                 nMomBins, TOFMatchEffConst::minMomentum, TOFMatchEffConst::maxMomentum);
    histograms.matchedVsP = new TH1F(Form("hMatchedVsP_%s", species.key),
                                     Form(";%s momentum p (GeV/c);tracks", species.label),
                                     nMomBins, TOFMatchEffConst::minMomentum, TOFMatchEffConst::maxMomentum);
    histograms.trueMatchedVsP = new TH1F(Form("hTrueMatchedVsP_%s", species.key),
                                         Form(";%s momentum p (GeV/c);tracks", species.label),
                                         nMomBins, TOFMatchEffConst::minMomentum, TOFMatchEffConst::maxMomentum);

    histograms.allVsPt = new TH1F(Form("hAllVsPt_%s", species.key),
                                  Form(";%s p_{T} (GeV/c);tracks", species.label),
                                  nPtBins, TOFMatchEffConst::minPt, TOFMatchEffConst::maxPt);
    histograms.matchedVsPt = new TH1F(Form("hMatchedVsPt_%s", species.key),
                                      Form(";%s p_{T} (GeV/c);tracks", species.label),
                                      nPtBins, TOFMatchEffConst::minPt, TOFMatchEffConst::maxPt);
    histograms.trueMatchedVsPt = new TH1F(Form("hTrueMatchedVsPt_%s", species.key),
                                          Form(";%s p_{T} (GeV/c);tracks", species.label),
                                          nPtBins, TOFMatchEffConst::minPt, TOFMatchEffConst::maxPt);

    histograms.allVsEta = new TH1F(Form("hAllVsEta_%s", species.key),
                                   Form(";%s #eta;tracks", species.label),
                                   nEtaBins, TOFMatchEffConst::minEta, TOFMatchEffConst::maxEta);
    histograms.matchedVsEta = new TH1F(Form("hMatchedVsEta_%s", species.key),
                                       Form(";%s #eta;tracks", species.label),
                                       nEtaBins, TOFMatchEffConst::minEta, TOFMatchEffConst::maxEta);
    histograms.trueMatchedVsEta = new TH1F(Form("hTrueMatchedVsEta_%s", species.key),
                                           Form(";%s #eta;tracks", species.label),
                                           nEtaBins, TOFMatchEffConst::minEta, TOFMatchEffConst::maxEta);

    histograms.allPtVsEta = new TH2F(Form("hAllPtVsEta_%s", species.key),
                                     Form(";%s p_{T} (GeV/c);#eta;tracks", species.label),
                                     nPtBins, TOFMatchEffConst::minPt, TOFMatchEffConst::maxPt,
                                     nEtaBins, TOFMatchEffConst::minEta, TOFMatchEffConst::maxEta);
    histograms.matchedPtVsEta = new TH2F(Form("hMatchedPtVsEta_%s", species.key),
                                         Form(";%s p_{T} (GeV/c);#eta;tracks", species.label),
                                         nPtBins, TOFMatchEffConst::minPt, TOFMatchEffConst::maxPt,
                                         nEtaBins, TOFMatchEffConst::minEta, TOFMatchEffConst::maxEta);
    histograms.trueMatchedPtVsEta = new TH2F(Form("hTrueMatchedPtVsEta_%s", species.key),
                                             Form(";%s p_{T} (GeV/c);#eta;tracks", species.label),
                                             nPtBins, TOFMatchEffConst::minPt, TOFMatchEffConst::maxPt,
                                             nEtaBins, TOFMatchEffConst::minEta, TOFMatchEffConst::maxEta);
    return histograms;
  };

  for (size_t iSpecies = 0; iSpecies < speciesConfigs.size(); ++iSpecies) {
    speciesHistograms[iSpecies] = bookSpeciesHistograms(speciesConfigs[iSpecies]);
  }

  auto getSpeciesIndex = [&](int absPdg) {
    for (size_t iSpecies = 0; iSpecies < speciesConfigs.size(); ++iSpecies) {
      if (speciesConfigs[iSpecies].absPdg == absPdg) {
        return static_cast<int>(iSpecies);
      }
    }
    return -1;
  };

  long long nAllTracks = 0;
  long long nMatchedTracks = 0;
  long long nTrueMatchedTracks = 0;
  long long nMatchedWithInvalidToFIndex = 0;
  long long nProtonTracks = 0;
  long long nProtonMatched = 0;
  long long nProtonTrueMatched = 0;

  for (int iEvent = 0; iEvent < nEvents; ++iEvent) {
    trackTree->GetEvent(iEvent);
    mcTree->GetEvent(iEvent);
    tClus  ->GetEvent(iEvent);

    for (const auto& track : tracks) {
      if (track.getNVTHits() < TOFMatchEffConst::minVerTelHits) {
        continue;
      }

      const double momentum = track.getP();
      const double pt = getTrackPt(track);
      const double eta = getTrackEta(track);
      const bool isMatched = track.getTOF() > 0.;
      ++nAllTracks;
      // retrive the particle ID of the tof cluster associated to the track, and check if it matches the track's particle ID (if available)
      int tofClusterPid = -1;
      if (isMatched) {
        ++nMatchedTracks;
        const int clusterIndex = track.getClusterIndex(nVTLayers);
        if (clusterIndex >= 0 && clusterIndex < static_cast<int>(tofClus.size())) {
          tofClusterPid = tofClus[clusterIndex].getParticleID();
        } else {
          ++nMatchedWithInvalidToFIndex;
        }
      }
      const bool isTrueMatched = isMatched && tofClusterPid == track.getParticleID();
      if (isTrueMatched) {
        ++nTrueMatchedTracks;
      }
      const int absTrackPid = std::abs(track.getParticleID());
      int speciesIndex = -1;
      if (mcArr && absTrackPid >= 0 && absTrackPid < static_cast<int>(mcArr->size())) {
        speciesIndex = getSpeciesIndex(std::abs(mcArr->at(absTrackPid).GetPdgCode()));
      }

      if (speciesIndex == 0) {
        ++nProtonTracks;
        if (isMatched) {
          ++nProtonMatched;
        }
        if (isTrueMatched) {
          ++nProtonTrueMatched;
        }
      }

      hAllVsP->Fill(momentum);
      hAllVsPt->Fill(pt);
      hAllVsEta->Fill(eta);
      hAllPtVsEta->Fill(pt, eta);

      if (isMatched) {
        hMatchedVsP->Fill(momentum);
        hMatchedVsPt->Fill(pt);
        hMatchedVsEta->Fill(eta);
        hMatchedPtVsEta->Fill(pt, eta);
      }

      if (isTrueMatched) {
        hTrueMatchedVsP->Fill(momentum);
        hTrueMatchedVsPt->Fill(pt);
        hTrueMatchedVsEta->Fill(eta);
        hTrueMatchedPtVsEta->Fill(pt, eta);
      }

      if (speciesIndex >= 0) {
        auto& species = speciesHistograms[speciesIndex];
        species.allVsP->Fill(momentum);
        species.allVsPt->Fill(pt);
        species.allVsEta->Fill(eta);
        species.allPtVsEta->Fill(pt, eta);
        if (isMatched) {
          species.matchedVsP->Fill(momentum);
          species.matchedVsPt->Fill(pt);
          species.matchedVsEta->Fill(eta);
          species.matchedPtVsEta->Fill(pt, eta);
        }
        if (isTrueMatched) {
          species.trueMatchedVsP->Fill(momentum);
          species.trueMatchedVsPt->Fill(pt);
          species.trueMatchedVsEta->Fill(eta);
          species.trueMatchedPtVsEta->Fill(pt, eta);
        }
      }
    }
  }

  auto* outputFile = TFile::Open(Form("%s/TOFMatchingEfficiency.root", dirSimu), "recreate");
  for (auto* hist : {hAllVsP, hMatchedVsP, hTrueMatchedVsP, hAllVsPt, hMatchedVsPt, hTrueMatchedVsPt,
                     hAllVsEta, hMatchedVsEta, hTrueMatchedVsEta}) {
    hist->Write();
  }
  hAllPtVsEta->Write();
  hMatchedPtVsEta->Write();
  hTrueMatchedPtVsEta->Write();
  for (const auto& species : speciesHistograms) {
    for (auto* hist : {species.allVsP, species.matchedVsP, species.trueMatchedVsP,
                       species.allVsPt, species.matchedVsPt, species.trueMatchedVsPt,
                       species.allVsEta, species.matchedVsEta, species.trueMatchedVsEta}) {
      hist->Write();
    }
    species.allPtVsEta->Write();
    species.matchedPtVsEta->Write();
    species.trueMatchedPtVsEta->Write();
  }

  auto cloneH1 = [](TH1F* hist, const char* newName) -> TH1F* {
    auto* clone = static_cast<TH1F*>(hist->Clone(newName));
    clone->SetDirectory(nullptr);
    return clone;
  };
  auto cloneH2 = [](TH2F* hist, const char* newName) -> TH2F* {
    auto* clone = static_cast<TH2F*>(hist->Clone(newName));
    clone->SetDirectory(nullptr);
    return clone;
  };

  auto* effMatchVsP = makeEfficiency1D(cloneH1(hMatchedVsP, "hMatchedVsP_cloned"), hAllVsP,
                                       "effTOFMatchVsP", "TOF matching efficiency;p (GeV/c);#varepsilon");
  auto* effTrueVsP = makeEfficiency1D(cloneH1(hTrueMatchedVsP, "hTrueMatchedVsP_cloned"), hAllVsP,
                                      "effTrueTOFMatchVsP", "True TOF matching efficiency;p (GeV/c);#varepsilon");
  auto* purityVsP = makeEfficiency1D(cloneH1(hTrueMatchedVsP, "hTrueMatchedVsP_forPurity"), hMatchedVsP,
                                     "purityTOFMatchVsP", "TOF matching purity;p (GeV/c);Purity");

  auto* effMatchVsPt = makeEfficiency1D(cloneH1(hMatchedVsPt, "hMatchedVsPt_cloned"), hAllVsPt,
                                        "effTOFMatchVsPt", "TOF matching efficiency;p_{T} (GeV/c);#varepsilon");
  auto* effTrueVsPt = makeEfficiency1D(cloneH1(hTrueMatchedVsPt, "hTrueMatchedVsPt_cloned"), hAllVsPt,
                                       "effTrueTOFMatchVsPt", "True TOF matching efficiency;p_{T} (GeV/c);#varepsilon");
  auto* purityVsPt = makeEfficiency1D(cloneH1(hTrueMatchedVsPt, "hTrueMatchedVsPt_forPurity"), hMatchedVsPt,
                                      "purityTOFMatchVsPt", "TOF matching purity;p_{T} (GeV/c);Purity");

  auto* effMatchVsEta = makeEfficiency1D(cloneH1(hMatchedVsEta, "hMatchedVsEta_cloned"), hAllVsEta,
                                         "effTOFMatchVsEta", "TOF matching efficiency;#eta;#varepsilon");
  auto* effTrueVsEta = makeEfficiency1D(cloneH1(hTrueMatchedVsEta, "hTrueMatchedVsEta_cloned"), hAllVsEta,
                                        "effTrueTOFMatchVsEta", "True TOF matching efficiency;#eta;#varepsilon");
  auto* purityVsEta = makeEfficiency1D(cloneH1(hTrueMatchedVsEta, "hTrueMatchedVsEta_forPurity"), hMatchedVsEta,
                                       "purityTOFMatchVsEta", "TOF matching purity;#eta;Purity");

  auto* effMatchPtVsEta = makeEfficiency2D(cloneH2(hMatchedPtVsEta, "hMatchedPtVsEta_cloned"), hAllPtVsEta,
                                           "effTOFMatchPtVsEta", "TOF matching efficiency;p_{T} (GeV/c);#eta");
  auto* effTruePtVsEta = makeEfficiency2D(cloneH2(hTrueMatchedPtVsEta, "hTrueMatchedPtVsEta_cloned"), hAllPtVsEta,
                                          "effTrueTOFMatchPtVsEta", "True TOF matching efficiency;p_{T} (GeV/c);#eta");
  auto* purityPtVsEta = makeEfficiency2D(cloneH2(hTrueMatchedPtVsEta, "hTrueMatchedPtVsEta_forPurity"), hMatchedPtVsEta,
                                         "purityTOFMatchPtVsEta", "TOF matching purity;p_{T} (GeV/c);#eta");

  struct SpeciesEfficiencies {
    TEfficiency* effMatchVsP = nullptr;
    TEfficiency* effTrueVsP = nullptr;
    TEfficiency* purityVsP = nullptr;
    TEfficiency* effMatchVsPt = nullptr;
    TEfficiency* effTrueVsPt = nullptr;
    TEfficiency* purityVsPt = nullptr;
    TEfficiency* effMatchVsEta = nullptr;
    TEfficiency* effTrueVsEta = nullptr;
    TEfficiency* purityVsEta = nullptr;
    TEfficiency* effMatchPtVsEta = nullptr;
    TEfficiency* effTruePtVsEta = nullptr;
    TEfficiency* purityPtVsEta = nullptr;
  };
  std::array<SpeciesEfficiencies, 3> speciesEfficiencies;

  for (size_t iSpecies = 0; iSpecies < speciesConfigs.size(); ++iSpecies) {
    const auto& cfg = speciesConfigs[iSpecies];
    const auto& hist = speciesHistograms[iSpecies];
    auto& eff = speciesEfficiencies[iSpecies];
    eff.effMatchVsP = makeEfficiency1D(cloneH1(hist.matchedVsP, Form("hMatchedVsP_%s_cloned", cfg.key)), hist.allVsP,
                                       Form("effTOFMatchVsP_%s", cfg.key),
                                       Form("%s TOF matching efficiency;p (GeV/c);#varepsilon", cfg.label));
    eff.effTrueVsP = makeEfficiency1D(cloneH1(hist.trueMatchedVsP, Form("hTrueMatchedVsP_%s_cloned", cfg.key)), hist.allVsP,
                                      Form("effTrueTOFMatchVsP_%s", cfg.key),
                                      Form("%s true TOF matching efficiency;p (GeV/c);#varepsilon", cfg.label));
    eff.purityVsP = makeEfficiency1D(cloneH1(hist.trueMatchedVsP, Form("hTrueMatchedVsP_%s_forPurity", cfg.key)), hist.matchedVsP,
                                     Form("purityTOFMatchVsP_%s", cfg.key),
                                     Form("%s TOF matching purity;p (GeV/c);Purity", cfg.label));

    eff.effMatchVsPt = makeEfficiency1D(cloneH1(hist.matchedVsPt, Form("hMatchedVsPt_%s_cloned", cfg.key)), hist.allVsPt,
                                        Form("effTOFMatchVsPt_%s", cfg.key),
                                        Form("%s TOF matching efficiency;p_{T} (GeV/c);#varepsilon", cfg.label));
    eff.effTrueVsPt = makeEfficiency1D(cloneH1(hist.trueMatchedVsPt, Form("hTrueMatchedVsPt_%s_cloned", cfg.key)), hist.allVsPt,
                                       Form("effTrueTOFMatchVsPt_%s", cfg.key),
                                       Form("%s true TOF matching efficiency;p_{T} (GeV/c);#varepsilon", cfg.label));
    eff.purityVsPt = makeEfficiency1D(cloneH1(hist.trueMatchedVsPt, Form("hTrueMatchedVsPt_%s_forPurity", cfg.key)), hist.matchedVsPt,
                                      Form("purityTOFMatchVsPt_%s", cfg.key),
                                      Form("%s TOF matching purity;p_{T} (GeV/c);Purity", cfg.label));

    eff.effMatchVsEta = makeEfficiency1D(cloneH1(hist.matchedVsEta, Form("hMatchedVsEta_%s_cloned", cfg.key)), hist.allVsEta,
                                         Form("effTOFMatchVsEta_%s", cfg.key),
                                         Form("%s TOF matching efficiency;#eta;#varepsilon", cfg.label));
    eff.effTrueVsEta = makeEfficiency1D(cloneH1(hist.trueMatchedVsEta, Form("hTrueMatchedVsEta_%s_cloned", cfg.key)), hist.allVsEta,
                                        Form("effTrueTOFMatchVsEta_%s", cfg.key),
                                        Form("%s true TOF matching efficiency;#eta;#varepsilon", cfg.label));
    eff.purityVsEta = makeEfficiency1D(cloneH1(hist.trueMatchedVsEta, Form("hTrueMatchedVsEta_%s_forPurity", cfg.key)), hist.matchedVsEta,
                                       Form("purityTOFMatchVsEta_%s", cfg.key),
                                       Form("%s TOF matching purity;#eta;Purity", cfg.label));

    eff.effMatchPtVsEta = makeEfficiency2D(cloneH2(hist.matchedPtVsEta, Form("hMatchedPtVsEta_%s_cloned", cfg.key)), hist.allPtVsEta,
                                           Form("effTOFMatchPtVsEta_%s", cfg.key),
                                           Form("%s TOF matching efficiency;p_{T} (GeV/c);#eta", cfg.label));
    eff.effTruePtVsEta = makeEfficiency2D(cloneH2(hist.trueMatchedPtVsEta, Form("hTrueMatchedPtVsEta_%s_cloned", cfg.key)), hist.allPtVsEta,
                                          Form("effTrueTOFMatchPtVsEta_%s", cfg.key),
                                          Form("%s true TOF matching efficiency;p_{T} (GeV/c);#eta", cfg.label));
    eff.purityPtVsEta = makeEfficiency2D(cloneH2(hist.trueMatchedPtVsEta, Form("hTrueMatchedPtVsEta_%s_forPurity", cfg.key)), hist.matchedPtVsEta,
                                         Form("purityTOFMatchPtVsEta_%s", cfg.key),
                                         Form("%s TOF matching purity;p_{T} (GeV/c);#eta", cfg.label));
  }

  for (auto* efficiency : {effMatchVsP, effTrueVsP, purityVsP, effMatchVsPt, effTrueVsPt, purityVsPt,
                           effMatchVsEta, effTrueVsEta, purityVsEta, effMatchPtVsEta, effTruePtVsEta, purityPtVsEta}) {
    if (efficiency) {
      efficiency->Write();
    }
  }
  for (const auto& eff : speciesEfficiencies) {
    for (auto* efficiency : {eff.effMatchVsP, eff.effTrueVsP, eff.purityVsP,
                             eff.effMatchVsPt, eff.effTrueVsPt, eff.purityVsPt,
                             eff.effMatchVsEta, eff.effTrueVsEta, eff.purityVsEta,
                             eff.effMatchPtVsEta, eff.effTruePtVsEta, eff.purityPtVsEta}) {
      if (efficiency) {
        efficiency->Write();
      }
    }
  }

  drawEfficiency1D(effMatchVsP, "cEffTOFMatchVsP", "TOF matching efficiency vs p",
                   Form("%s/tof_matching_eff_vs_p.png", dirSimu));
  drawEfficiency1D(effTrueVsP, "cEffTrueTOFMatchVsP", "True TOF matching efficiency vs p",
                   Form("%s/tof_true_matching_eff_vs_p.png", dirSimu));
  drawEfficiency1D(purityVsP, "cPurityTOFMatchVsP", "TOF matching purity vs p",
                   Form("%s/tof_matching_purity_vs_p.png", dirSimu));

  drawEfficiency1D(effMatchVsPt, "cEffTOFMatchVsPt", "TOF matching efficiency vs p_{T}",
                   Form("%s/tof_matching_eff_vs_pt.png", dirSimu));
  drawEfficiency1D(effTrueVsPt, "cEffTrueTOFMatchVsPt", "True TOF matching efficiency vs p_{T}",
                   Form("%s/tof_true_matching_eff_vs_pt.png", dirSimu));
  drawEfficiency1D(purityVsPt, "cPurityTOFMatchVsPt", "TOF matching purity vs p_{T}",
                   Form("%s/tof_matching_purity_vs_pt.png", dirSimu));

  drawEfficiency1D(effMatchVsEta, "cEffTOFMatchVsEta", "TOF matching efficiency vs #eta",
                   Form("%s/tof_matching_eff_vs_eta.png", dirSimu));
  drawEfficiency1D(effTrueVsEta, "cEffTrueTOFMatchVsEta", "True TOF matching efficiency vs #eta",
                   Form("%s/tof_true_matching_eff_vs_eta.png", dirSimu));
  drawEfficiency1D(purityVsEta, "cPurityTOFMatchVsEta", "TOF matching purity vs #eta",
                   Form("%s/tof_matching_purity_vs_eta.png", dirSimu));

  drawEfficiency2D(effMatchPtVsEta, "cEffTOFMatchPtVsEta", "TOF matching efficiency vs (p_{T}, #eta)",
                   Form("%s/tof_matching_eff_pt_vs_eta.png", dirSimu));
  drawEfficiency2D(effTruePtVsEta, "cEffTrueTOFMatchPtVsEta", "True TOF matching efficiency vs (p_{T}, #eta)",
                   Form("%s/tof_true_matching_eff_pt_vs_eta.png", dirSimu));
  drawEfficiency2D(purityPtVsEta, "cPurityTOFMatchPtVsEta", "TOF matching purity vs (p_{T}, #eta)",
                   Form("%s/tof_matching_purity_pt_vs_eta.png", dirSimu));

  for (size_t iSpecies = 0; iSpecies < speciesConfigs.size(); ++iSpecies) {
    const auto& cfg = speciesConfigs[iSpecies];
    const auto& eff = speciesEfficiencies[iSpecies];
    drawEfficiency1D(eff.effMatchVsP, Form("cEffTOFMatchVsP_%s", cfg.key),
                     Form("%s TOF matching efficiency vs p", cfg.label),
                     Form("%s/tof_matching_eff_vs_p_%s.png", dirSimu, cfg.key));
    drawEfficiency1D(eff.effTrueVsP, Form("cEffTrueTOFMatchVsP_%s", cfg.key),
                     Form("%s true TOF matching efficiency vs p", cfg.label),
                     Form("%s/tof_true_matching_eff_vs_p_%s.png", dirSimu, cfg.key));
    drawEfficiency1D(eff.purityVsP, Form("cPurityTOFMatchVsP_%s", cfg.key),
                     Form("%s TOF matching purity vs p", cfg.label),
                     Form("%s/tof_matching_purity_vs_p_%s.png", dirSimu, cfg.key));

    drawEfficiency1D(eff.effMatchVsPt, Form("cEffTOFMatchVsPt_%s", cfg.key),
                     Form("%s TOF matching efficiency vs p_{T}", cfg.label),
                     Form("%s/tof_matching_eff_vs_pt_%s.png", dirSimu, cfg.key));
    drawEfficiency1D(eff.effTrueVsPt, Form("cEffTrueTOFMatchVsPt_%s", cfg.key),
                     Form("%s true TOF matching efficiency vs p_{T}", cfg.label),
                     Form("%s/tof_true_matching_eff_vs_pt_%s.png", dirSimu, cfg.key));
    drawEfficiency1D(eff.purityVsPt, Form("cPurityTOFMatchVsPt_%s", cfg.key),
                     Form("%s TOF matching purity vs p_{T}", cfg.label),
                     Form("%s/tof_matching_purity_vs_pt_%s.png", dirSimu, cfg.key));

    drawEfficiency1D(eff.effMatchVsEta, Form("cEffTOFMatchVsEta_%s", cfg.key),
                     Form("%s TOF matching efficiency vs #eta", cfg.label),
                     Form("%s/tof_matching_eff_vs_eta_%s.png", dirSimu, cfg.key));
    drawEfficiency1D(eff.effTrueVsEta, Form("cEffTrueTOFMatchVsEta_%s", cfg.key),
                     Form("%s true TOF matching efficiency vs #eta", cfg.label),
                     Form("%s/tof_true_matching_eff_vs_eta_%s.png", dirSimu, cfg.key));
    drawEfficiency1D(eff.purityVsEta, Form("cPurityTOFMatchVsEta_%s", cfg.key),
                     Form("%s TOF matching purity vs #eta", cfg.label),
                     Form("%s/tof_matching_purity_vs_eta_%s.png", dirSimu, cfg.key));

    drawEfficiency2D(eff.effMatchPtVsEta, Form("cEffTOFMatchPtVsEta_%s", cfg.key),
                     Form("%s TOF matching efficiency vs (p_{T}, #eta)", cfg.label),
                     Form("%s/tof_matching_eff_pt_vs_eta_%s.png", dirSimu, cfg.key));
    drawEfficiency2D(eff.effTruePtVsEta, Form("cEffTrueTOFMatchPtVsEta_%s", cfg.key),
                     Form("%s true TOF matching efficiency vs (p_{T}, #eta)", cfg.label),
                     Form("%s/tof_true_matching_eff_pt_vs_eta_%s.png", dirSimu, cfg.key));
    drawEfficiency2D(eff.purityPtVsEta, Form("cPurityTOFMatchPtVsEta_%s", cfg.key),
                     Form("%s TOF matching purity vs (p_{T}, #eta)", cfg.label),
                     Form("%s/tof_matching_purity_pt_vs_eta_%s.png", dirSimu, cfg.key));
  }

  outputFile->Close();
    printf("TOF matching summary (selected tracks): all=%lld matched=%lld true=%lld invalidTOFidx=%lld\n",
      nAllTracks, nMatchedTracks, nTrueMatchedTracks, nMatchedWithInvalidToFIndex);
    if (nAllTracks > 0) {
      printf("  match eff = %.4f, true-match eff = %.4f, purity = %.4f\n",
        double(nMatchedTracks) / double(nAllTracks),
        double(nTrueMatchedTracks) / double(nAllTracks),
        nMatchedTracks > 0 ? double(nTrueMatchedTracks) / double(nMatchedTracks) : 0.0);
    }
    printf("Proton summary: all=%lld matched=%lld true=%lld\n",
      nProtonTracks, nProtonMatched, nProtonTrueMatched);
    if (nProtonTracks > 0) {
      printf("  proton match eff = %.4f, proton true-match eff = %.4f, proton purity = %.4f\n",
        double(nProtonMatched) / double(nProtonTracks),
        double(nProtonTrueMatched) / double(nProtonTracks),
        nProtonMatched > 0 ? double(nProtonTrueMatched) / double(nProtonMatched) : 0.0);
    }
  mcTree->GetCurrentFile()->Close();
    tClus->GetCurrentFile()->Close();
  trackTree->GetCurrentFile()->Close();
  printf("Done — TOF matching efficiency plots saved in %s/\n", dirSimu);
}