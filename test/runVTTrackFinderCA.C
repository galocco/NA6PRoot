/*
  CA Track Finder Macro

  Usage:
    root -l
    root [0] .L setup_macros.C     // Load NA6PRoot libraries and setup paths
    root [1] .L runTrackFinderCA.C+  // Compile the macro with ACLiC
    root [2] runTrackFinderCA()      // Run with default parameters

  Or with custom parameters:
    root [2] runTrackFinderCA(0, 100, "../fullgeo", 5)
*/

#if !defined(__CINT__) || defined(__MAKECINT__)
#include <TTree.h>
#include <TFile.h>
#include <TH1F.h>
#include <TVector3.h>
#include <TParticle.h>
#include <TCanvas.h>
#include <TLegend.h>
#include "NA6PVerTelCluster.h"
#include "NA6PTrack.h"
#include "NA6PVertex.h"
#include "NA6PTrackerCA.h"
#include "MagneticField.h"
#include "NA6PVerTelHit.h"
#include "NA6PMCTruthContainer.h"
#include "NA6PReconstruction.h"
#endif

const int maxIterationsCA = NA6PTrackerCA::kMaxIterationsCA;

void runVTTrackFinderCA(int firstEv = 0,
                        int lastEv = 9999,
                        const char* dirSimu = ".",
                        int nLayers = 5)
{

  auto magField = new MagneticField();
  if (!Propagator::loadField() || !Propagator::loadGeometry(Form("%s/geometry.root", dirSimu))) {
    return;
  }

  int nMomBins = 40;
  TH1F* hMomGen = new TH1F("hMomGen", ";p (GeV/c);counts", nMomBins, 0., 10.);
  TH1F* hEtaGen = new TH1F("hEtaGen", ";#eta;counts", 20, 1., 5.);
  TH1F* hMomReco = new TH1F("hMomReco", ";p (GeV/c);counts", nMomBins, 0., 10.);
  TH1F* hEtaReco = new TH1F("hEtaReco", ";#eta;counts", 20, 1., 5.);
  TH1F* hMomGoodReco = new TH1F("hMomGoodReco", ";p (GeV/c);counts", nMomBins, 0., 10.);
  TH1F* hEtaGoodReco = new TH1F("hEtaGoodReco", ";#eta;counts", 20, 1., 5.);
  TH1F** hMomRecoIterCA = new TH1F*[maxIterationsCA];
  TH1F** hEtaRecoIterCA = new TH1F*[maxIterationsCA];
  int colors[maxIterationsCA] = {kMagenta + 1, kBlue, kGreen + 1, kOrange + 1, kRed + 1, kRed, kRed - 9, kGray + 2, kGray + 1, kGray};
  for (int jIteration = 0; jIteration < maxIterationsCA; ++jIteration) {
    hMomRecoIterCA[jIteration] = new TH1F(Form("hMomRecoIterCA%d", jIteration), ";p (GeV/c);counts", nMomBins, 0., 10.);
    hEtaRecoIterCA[jIteration] = new TH1F(Form("hEtaRecoIterCA%d", jIteration), ";#eta;counts", 20, 1., 5.);
  }

  std::unique_ptr<NA6PTrackerCA> tracker = std::make_unique<NA6PTrackerCA>();
  tracker->configureFromRecoParamVT();
  //  tracker->setNumberOfIterations(1);
  //  tracker->setIterationParams(0, 0.04, 0.1, 4., 0.4, 0.02, 0.002, 5000, 5000, 5000, 5);
  //  tracker->setIterationParams(1, 0.10, 0.2, 9., 0.6, 0.05, 0.0075, 1000, 1000, 1000, 3);
  tracker->printConfiguration();

  //  tracker->setVerbosity(true);
  TFile* fk = new TFile(Form("%s/MCKine.root", dirSimu));
  TTree* mcTree = (TTree*)fk->Get("mckine");
  std::vector<TParticle>* mcArr = nullptr;
  mcTree->SetBranchAddress("tracks", &mcArr);

  TFile* fh = new TFile(Form("%s/HitsVerTel.root", dirSimu));
  TTree* th = (TTree*)fh->Get("hitsVerTel");
  std::vector<NA6PVerTelHit> vtHits, *vtHitsPtr = &vtHits;
  th->SetBranchAddress("VerTel", &vtHitsPtr);

  TFile* fc = new TFile(Form("%s/ClustersVerTel.root", dirSimu));
  printf("Open cluster file: %s\n", fc->GetName());
  TTree* tc = (TTree*)fc->Get("clustersVerTel");
  std::vector<NA6PVerTelCluster> vtClus, *vtClusPtr = &vtClus;
  NA6PMCTruthContainer vtCluMCLabels, *vtCluMCLabelsPtr = &vtCluMCLabels;
  tc->SetBranchAddress("VerTel", &vtClusPtr);
  tc->SetBranchAddress("VerTelMCTruth", &vtCluMCLabelsPtr);
  int nEv = tc->GetEntries();
  if (lastEv > nEv || lastEv < 0)
    lastEv = nEv;
  if (firstEv < 0)
    firstEv = 0;

  NA6PVertex primVert;
  NA6PReconstruction rec("dummy");
  std::vector<NA6PMCComposedLabel> trkMCLabs;

  int nIterationsCA = tracker->getNIterations();

  for (int jEv = firstEv; jEv < lastEv; jEv++) {
    mcTree->GetEvent(jEv);
    th->GetEvent(jEv);
    int nPart = mcArr->size();
    double zvert = 0;
    // get primary vertex position from the Kine Tree
    for (int jp = 0; jp < nPart; jp++) {
      auto curPart = mcArr->at(jp);
      if (curPart.IsPrimary()) {
        zvert = curPart.Vz();
        break;
      }
    }
    primVert.setXYZ(0., 0., zvert);
    uint nHits = vtHits.size();
    std::vector<bool> isTrackable(nPart, false);
    std::vector<int> hitMapParticle(nPart, 0);
    for (int jp = 0; jp < nPart; jp++) {
      auto curPart = mcArr->at(jp);
      int maskHits = 0;
      for (size_t jHit = 0; jHit < nHits; ++jHit) {
        const auto& hit = vtHits.at(jHit);
        int idPart = hit.getTrackID();
        if (idPart == jp) {
          int nLay = hit.getDetectorID() / 4;
          maskHits |= (1 << nLay);
        }
      }
      hitMapParticle[jp] = maskHits;
      if (maskHits == (1 << nLayers) - 1) {
        isTrackable[jp] = true;
        double pxPart = curPart.Px();
        double pyPart = curPart.Py();
        double pzPart = curPart.Pz();
        double momPart = curPart.P();
        double phiPart = curPart.Phi();
        double thetaPart = std::acos(pzPart / momPart);
        double etaPart = -std::log(std::tan(thetaPart / 2.));
        hMomGen->Fill(momPart);
        hEtaGen->Fill(etaPart);
      }
    }
    tc->GetEvent(jEv);
    for (size_t i = 0; i < vtClus.size(); ++i) {
      vtClus[i].setClusterIndex(static_cast<int>(i));
    }
    tracker->setClusterMCTruth(&vtCluMCLabels);
    tracker->findTracks(vtClus, &primVert);
    std::vector<NA6PTrack> trks = tracker->getTracks();
    rec.assignMCLabels(trks, trkMCLabs, vtCluMCLabels);
    int nTrks = trks.size();
    for (int jT = 0; jT < nTrks; jT++) {
      NA6PTrack tr = trks[jT];
      NA6PMCComposedLabel mcCompLabel = trkMCLabs[jT];
      int idPartTrack = mcCompLabel.getTrackID();
      int jIteration = tr.getCAIteration();
      if (tr.getNHits() == 5) {
        if (!isTrackable[idPartTrack] && !mcCompLabel.isFake()) {
          printf("Mismatch: track %d has 5 hits, label %d, but hit map %d\n", jT, idPartTrack, hitMapParticle[idPartTrack]);
        }
        auto curPart = mcArr->at(idPartTrack);
        double pxPart = curPart.Px();
        double pyPart = curPart.Py();
        double pzPart = curPart.Pz();
        double momPart = curPart.P();
        double phiPart = curPart.Phi();
        double thetaPart = std::acos(pzPart / momPart);
        double etaPart = -std::log(std::tan(thetaPart / 2.));
        hMomReco->Fill(momPart);
        hEtaReco->Fill(etaPart);
        if (jIteration >= 0 && jIteration < maxIterationsCA) {
          hMomRecoIterCA[jIteration]->Fill(momPart);
          hEtaRecoIterCA[jIteration]->Fill(etaPart);
        }
        if (!mcCompLabel.isFake()) {
          hMomGoodReco->Fill(momPart);
          hEtaGoodReco->Fill(etaPart);
        }
      }
    }
  }

  TH1F* hPurityMom = (TH1F*)hMomGoodReco->Clone("hPurityMom");
  hPurityMom->GetYaxis()->SetTitle("purity");
  for (int iBin = 1; iBin <= hMomGoodReco->GetNbinsX(); iBin++) {
    double cg = hMomGoodReco->GetBinContent(iBin);
    double ct = hMomReco->GetBinContent(iBin);
    if (ct == 0) {
      hPurityMom->SetBinContent(iBin, 0.);
      hPurityMom->SetBinError(iBin, 0.);
    } else {
      double p = cg / ct;
      double ep = std::sqrt(p * (1 - p) / ct);
      hPurityMom->SetBinContent(iBin, p);
      hPurityMom->SetBinError(iBin, ep);
    }
  }
  TH1F* hPurityEta = (TH1F*)hEtaGoodReco->Clone("hPurityEta");
  hPurityEta->GetYaxis()->SetTitle("purity");
  for (int iBin = 1; iBin <= hEtaGoodReco->GetNbinsX(); iBin++) {
    double cg = hEtaGoodReco->GetBinContent(iBin);
    double ct = hEtaReco->GetBinContent(iBin);
    if (ct == 0) {
      hPurityEta->SetBinContent(iBin, 0.);
      hPurityEta->SetBinError(iBin, 0.);
    } else {
      double p = cg / ct;
      double ep = std::sqrt(p * (1 - p) / ct);
      hPurityEta->SetBinContent(iBin, p);
      hPurityEta->SetBinError(iBin, ep);
    }
  }

  TCanvas* cef = new TCanvas("cef", "", 1400, 800);
  cef->Divide(2, 2);
  cef->cd(1);
  hMomGen->SetLineColor(kGray + 1);
  hMomGen->SetLineWidth(3);
  hMomGen->Draw();
  hMomReco->SetLineWidth(2);
  hMomReco->Draw("same");
  TLegend* leg = new TLegend(0.5, 0.6, 0.89, 0.8);
  leg->AddEntry(hMomGen, "Generated, 5 hits in VT");
  for (int jIteration = 0; jIteration < nIterationsCA; ++jIteration) {
    if (hMomRecoIterCA[jIteration]->GetEntries() > 0) {
      hMomRecoIterCA[jIteration]->SetLineColor(colors[jIteration]);
      hMomRecoIterCA[jIteration]->SetLineWidth(2);
      hMomRecoIterCA[jIteration]->Draw("same");
      leg->AddEntry(hMomRecoIterCA[jIteration], Form("Reconstructed Interation %d, 5 hits in VT", jIteration));
    }
  }
  leg->Draw();
  cef->cd(2);
  TH1F* hEffMom = (TH1F*)hMomReco->Clone("hEffMom");
  hEffMom->Divide(hMomReco, hMomGen, 1., 1., "B");
  hEffMom->GetYaxis()->SetTitle("Efficiency");
  hEffMom->SetStats(0);
  hEffMom->Draw();
  cef->cd(3);
  hEtaGen->SetLineColor(kGray + 1);
  hEtaGen->SetLineWidth(3);
  hEtaGen->Draw();
  hEtaReco->SetLineWidth(2);
  hEtaReco->Draw("same");
  for (int jIteration = 0; jIteration < nIterationsCA; ++jIteration) {
    if (hEtaRecoIterCA[jIteration]->GetEntries() > 0) {
      hEtaRecoIterCA[jIteration]->SetLineColor(colors[jIteration]);
      hEtaRecoIterCA[jIteration]->SetLineWidth(2);
      hEtaRecoIterCA[jIteration]->Draw("same");
    }
  }
  cef->cd(4);
  TH1F* hEffEta = (TH1F*)hEtaReco->Clone("hEffEta");
  hEffEta->Divide(hEtaReco, hEtaGen, 1., 1., "B");
  hEffEta->GetYaxis()->SetTitle("Efficiency");
  hEffEta->SetStats(0);
  hEffEta->Draw();
  cef->SaveAs("TrackingEffic-VT.png");

  TCanvas* cpu = new TCanvas("cpu", "", 1400, 800);
  cpu->Divide(2, 2);
  cpu->cd(1);
  hMomReco->Draw();
  hMomGoodReco->SetLineColor(kGreen + 1);
  hMomGoodReco->SetLineWidth(2);
  hMomGoodReco->Draw("same");
  cpu->cd(2);
  hPurityMom->GetYaxis()->SetTitle("Purity");
  hPurityMom->SetMinimum(0.8);
  hPurityMom->SetStats(0);
  hPurityMom->SetLineWidth(2);
  hPurityMom->Draw();
  cpu->cd(3);
  hEtaReco->Draw();
  hEtaGoodReco->SetLineColor(kGreen + 1);
  hEtaGoodReco->SetLineWidth(2);
  hEtaGoodReco->Draw("same");
  cpu->cd(4);
  hPurityEta->GetYaxis()->SetTitle("Purity");
  hPurityEta->SetMinimum(0.8);
  hPurityEta->SetStats(0);
  hPurityEta->SetLineWidth(2);
  hPurityEta->Draw();
  cpu->SaveAs("TrackingPurity-VT.png");
}
