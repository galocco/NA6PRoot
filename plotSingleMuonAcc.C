#if !defined(__CINT__) || defined(__MAKECINT__)
#include <TTree.h>
#include <TFile.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TF1.h>
#include <TParticle.h>
#include <TCanvas.h>
#include <TPaveStats.h>
#include <TLegend.h>
#include <TEfficiency.h>
#include <TGeoManager.h>
#include <TStyle.h>
#include <TLorentzVector.h>
#include "NA6PTrack.h"
#include "MagneticField.h"
#include "NA6PVerTelHit.h"
#include "NA6PMuonSpecModularHit.h"
#include "NA6PLayoutParam.h"
#include "NA6PFastTrackFitter.h"
#endif

// function to set title size and marker style for a histogram

void setStyle(TH1* hist)
{
  hist->SetMarkerSize(0.8);
  hist->SetTitleSize(0.05, "X");
  hist->SetTitleSize(0.05, "Y");
}

void setStyle(TEfficiency* eff)
{
  if (!eff) return;

  eff->SetMarkerSize(0.8);
  eff->SetMarkerStyle(20);
  eff->SetLineWidth(2);
}

void setStyleAxes(TEfficiency* eff)
{
  if (!eff || !gPad) return;

  gPad->Modified();
  gPad->Update(); // necessario per creare il painted graph
  auto* gr = eff->GetPaintedGraph();
  if (!gr) return;
  gr->SetMarkerSize(0.8);
  gr->GetXaxis()->SetTitleSize(0.06);
  gr->GetYaxis()->SetTitleSize(0.06);
  gr->GetXaxis()->SetLabelSize(0.045);
  gr->GetYaxis()->SetLabelSize(0.045);
  gr->GetXaxis()->SetTitleOffset(1.1);
  gr->GetYaxis()->SetTitleOffset(1.15);
  gPad->Modified();
  gPad->Update();
}

void plotSingleMuonAcc(const char* dirSimu = ".", int NRecoClusters = 4)
{
  // --- Initialization & Field Loading ---
  auto& param = NA6PLayoutParam::Instance();
  param.updateFromFile(Form("%s/na6pLayout.ini", dirSimu), "", true);

  auto magField = new MagneticField();
  magField->loadField();
  magField->setAsGlobalField();

  NA6PFastTrackFitter fitter;
  fitter.loadGeometry();
  fitter.enableMaterialCorrections();

  // --- File Loading ---
  TFile* fVT = TFile::Open(Form("%s/TracksVerTel.root", dirSimu));
  if (!fVT || fVT->IsZombie())
    return;
  TTree* trTreeVT = (TTree*)fVT->Get("tracksVerTel");
  if (!trTreeVT)
    return;
  std::vector<NA6PTrack>* trArrVT = nullptr;
  trTreeVT->SetBranchAddress("VerTel", &trArrVT);

  TFile* fMS = TFile::Open(Form("%s/TracksMuonSpec.root", dirSimu));
  if (!fMS || fMS->IsZombie())
    return;
  TTree* trTreeMS = (TTree*)fMS->Get("tracksMuonSpec");
  if (!trTreeMS)
    return;
  std::vector<NA6PTrack>* trArrMS = nullptr;
  trTreeMS->SetBranchAddress("MuonSpec", &trArrMS);

  TFile* fk = TFile::Open(Form("%s/MCKine.root", dirSimu));
  if (!fk || fk->IsZombie())
    return;
  TTree* mcTree = (TTree*)fk->Get("mckine");
  if (!mcTree)
    return;
  std::vector<TParticle>* mcArr = nullptr;
  mcTree->SetBranchAddress("tracks", &mcArr);

  float y0 = 2.2;
  float range = 5;
  // --- Histogram Definitions (Restored) ---
  TH1D* hRapidityVT = new TH1D("hRapidityVT", " Muons from J/#psi;#it{y}; Counts", 100, y0 - range, y0 + range);
  TH1D* hRapidityVTNeg = new TH1D("hRapidityVTNeg", " Muons from J/#psi;#it{y}; Counts", 100, y0 - range, y0 + range);
  TH1D* hRapidityVTPos = new TH1D("hRapidityVTPos", " Muons from J/#psi;#it{y}; Counts", 100, y0 - range, y0 + range);
  TH1D* hRapidityMS = new TH1D("hRapidityMS", " Muons from J/#psi;#it{y}; Counts", 100, y0 - range, y0 + range);
  TH1D* hRapidityMSNeg = new TH1D("hRapidityMSNeg", " Muons from J/#psi;#it{y}; Counts", 100, y0 - range, y0 + range);
  TH1D* hRapidityMSPos = new TH1D("hRapidityMSPos", " Muons from J/#psi;#it{y}; Counts", 100, y0 - range, y0 + range);
  TH1D* hRapidityMC = new TH1D("hRapidityMC", " Muons from J/#psi;#it{y}; Counts", 100, y0 - range, y0 + range);
  TH1D* hRapidityMC2 = new TH1D("hRapidityMC2", " Muons from J/#psi;#it{y}; Counts", 100, y0 - range, y0 + range);
  TH1D* hIsPIDGoodVT = new TH1D("hIsPIDGoodVT", " Is PID good?;Is PID good?; Counts", 2, -0.5, 1.5);
  TH1D* hIsPIDGoodMS = new TH1D("hIsPIDGoodMS", " Is PID good?;Is PID good?; Counts", 2, -0.5, 1.5);
  int nEv = std::min((int)trTreeVT->GetEntries(), std::min((int)trTreeMS->GetEntries(), (int)mcTree->GetEntries()));
  printf("Number of events = %d\n", nEv);

  for (int jEv = 0; jEv < nEv; jEv++) {
    mcTree->GetEvent(jEv);
    trTreeVT->GetEvent(jEv);
    trTreeMS->GetEvent(jEv);

    int nPart = (mcArr) ? mcArr->size() : 0;
    int nTracksVT = (trArrVT) ? trArrVT->size() : 0;
    int nTracksMS = (trArrMS) ? trArrMS->size() : 0;
    printf("Event %d particles = %d tracks VT = %d tracks MS = %d\n", jEv, nPart, nTracksVT, nTracksMS);

    // loop over the MCkine to fille the rapidity of muons
    for (int jp = 0; jp < nPart; jp++) {
      auto& curPart = mcArr->at(jp);
      if (TMath::Abs(curPart.GetPdgCode()) == 13) {
        double px = curPart.Px();
        double py = curPart.Py();
        double pz = curPart.Pz();
        double p = TMath::Sqrt(px * px + py * py + pz * pz);
        double rapidity = 0.5 * TMath::Log((p + pz) / (p - pz));
        hRapidityMC->Fill(rapidity);
        hRapidityMC2->Fill(rapidity);
      }
    }

    for (int jTr = 0; jTr < nTracksVT; ++jTr) {
      NA6PTrack& tr = trArrVT->at(jTr);
      int nClusters = tr.getNHits();
      if (TMath::Abs(tr.getParticleID()) >= nPart)
        hIsPIDGoodVT->Fill(1);
      else
        hIsPIDGoodVT->Fill(0);
      if (nClusters < NRecoClusters)
        continue;
      if (tr.getParticleID() < 0 || tr.getParticleID() >= nPart)
        continue; // check if the track is a muon

      auto& curPart = mcArr->at(tr.getParticleID());
      if (TMath::Abs(curPart.GetPdgCode()) == 13) {

        double px = curPart.Px();
        double py = curPart.Py();
        double pz = curPart.Pz();
        double p = TMath::Sqrt(px * px + py * py + pz * pz);
        double rapidity = 0.5 * TMath::Log((p + pz) / (p - pz));
        // check if the track is a muon
        hRapidityVT->Fill(rapidity);
        if (tr.getCharge() < 0)
          hRapidityVTNeg->Fill(rapidity);
        if (tr.getCharge() > 0)
          hRapidityVTPos->Fill(rapidity);
      }
    }

    for (int jTr = 0; jTr < nTracksMS; ++jTr) {
      NA6PTrack& tr = trArrMS->at(jTr);
      int nClusters = tr.getNHits();
      if (TMath::Abs(tr.getParticleID()) >= nPart)
        hIsPIDGoodMS->Fill(1);
      else
        hIsPIDGoodMS->Fill(0);

      if (nClusters < NRecoClusters)
        continue;
      if (tr.getParticleID() < 0 || tr.getParticleID() >= nPart)
        continue; // check if the track is a muon

      auto& curPart = mcArr->at(tr.getParticleID());
      if (TMath::Abs(curPart.GetPdgCode()) == 13) {

        double px = curPart.Px();
        double py = curPart.Py();
        double pz = curPart.Pz();
        double p = TMath::Sqrt(px * px + py * py + pz * pz);
        double rapidity = 0.5 * TMath::Log((p + pz) / (p - pz));
        // check if the track is a muon
        hRapidityMS->Fill(rapidity);
        if (tr.getCharge() < 0)
          hRapidityMSNeg->Fill(rapidity);
        if (tr.getCharge() > 0)
          hRapidityMSPos->Fill(rapidity);
      }
    }

  }

  hRapidityMC2->Scale(0.5); // account for both mu+ and mu-
  gStyle->SetOptFit(1111);

  TFile fOut(Form("%s/single_muon.root", dirSimu), "RECREATE");
  TEfficiency* hEffVsRapidityVT = nullptr;
  TEfficiency* hEffVsRapidityVTNeg = nullptr;
  TEfficiency* hEffVsRapidityVTPos = nullptr;
  TEfficiency* hEffVsRapidityMS = nullptr;
  TEfficiency* hEffVsRapidityMSNeg = nullptr;
  TEfficiency* hEffVsRapidityMSPos = nullptr;

  if (TEfficiency::CheckConsistency(*hRapidityVT, *hRapidityMC)) {
    hEffVsRapidityVT = new TEfficiency(*hRapidityVT, *hRapidityMC);
    hEffVsRapidityVT->SetName("hEffVsRapidityVT");
    hEffVsRapidityVT->SetTitle("Efficiency x Acceptance for muons from J/#psi;#it{y};#epsilon #times A");
  }
  if (TEfficiency::CheckConsistency(*hRapidityVTNeg, *hRapidityMC2)) {
    hEffVsRapidityVTNeg = new TEfficiency(*hRapidityVTNeg, *hRapidityMC2);
    hEffVsRapidityVTNeg->SetName("hEffVsRapidityVTNeg");
    hEffVsRapidityVTNeg->SetTitle("Efficiency x Acceptance for muons from J/#psi;#it{y};#epsilon #times A");
  }
  if (TEfficiency::CheckConsistency(*hRapidityVTPos, *hRapidityMC2)) {
    hEffVsRapidityVTPos = new TEfficiency(*hRapidityVTPos, *hRapidityMC2);
    hEffVsRapidityVTPos->SetName("hEffVsRapidityVTPos");
    hEffVsRapidityVTPos->SetTitle("Efficiency x Acceptance for muons from J/#psi;#it{y};#epsilon #times A");
  }
  if (TEfficiency::CheckConsistency(*hRapidityMS, *hRapidityMC)) {
    hEffVsRapidityMS = new TEfficiency(*hRapidityMS, *hRapidityMC);
    hEffVsRapidityMS->SetName("hEffVsRapidityMS");
    hEffVsRapidityMS->SetTitle("Efficiency x Acceptance for muons from J/#psi;#it{y};#epsilon #times A");
  }
  if (TEfficiency::CheckConsistency(*hRapidityMSNeg, *hRapidityMC2)) {
    hEffVsRapidityMSNeg = new TEfficiency(*hRapidityMSNeg, *hRapidityMC2);
    hEffVsRapidityMSNeg->SetName("hEffVsRapidityMSNeg");
    hEffVsRapidityMSNeg->SetTitle("Efficiency x Acceptance for muons from J/#psi;#it{y};#epsilon #times A");
  }
  if (TEfficiency::CheckConsistency(*hRapidityMSPos, *hRapidityMC2)) {
    hEffVsRapidityMSPos = new TEfficiency(*hRapidityMSPos, *hRapidityMC2);
    hEffVsRapidityMSPos->SetName("hEffVsRapidityMSPos");
    hEffVsRapidityMSPos->SetTitle("Efficiency x Acceptance for muons from J/#psi;#it{y};#epsilon #times A");
  }

  setStyle(hRapidityVT);
  setStyle(hRapidityVTNeg);
  setStyle(hRapidityVTPos);
  setStyle(hRapidityMS);
  setStyle(hRapidityMSNeg);
  setStyle(hRapidityMSPos);

  hIsPIDGoodVT->Write();
  hIsPIDGoodMS->Write();
  hRapidityVT->Write();
  hRapidityVTNeg->Write();
  hRapidityVTPos->Write();
  hRapidityMS->Write();
  hRapidityMSNeg->Write();
  hRapidityMSPos->Write();

  if (hEffVsRapidityVT)
    hEffVsRapidityVT->Write();
  if (hEffVsRapidityVTPos)
    hEffVsRapidityVTPos->Write();
  if (hEffVsRapidityVTNeg)
    hEffVsRapidityVTNeg->Write();
  if (hEffVsRapidityMS)
    hEffVsRapidityMS->Write();
  if (hEffVsRapidityMSPos)
    hEffVsRapidityMSPos->Write();
  if (hEffVsRapidityMSNeg)
    hEffVsRapidityMSNeg->Write();

  hRapidityMC->Write();

  gStyle->SetOptStat(0);
  TCanvas* c1 = new TCanvas("c1", "Single Muon Rapidity", 800, 600);
  // set margins
  c1->SetLeftMargin(0.15);
  c1->SetRightMargin(0.05);
  c1->SetBottomMargin(0.15);

  float effVT = hRapidityVT->GetEntries() / hRapidityMC->GetEntries();
  float effMS = hRapidityMS->GetEntries() / hRapidityMC->GetEntries();
  hRapidityVT->SetLineColor(kBlue);
  hRapidityVT->SetMarkerColor(kBlue);
  hRapidityVT->SetMarkerStyle(20);
  hRapidityVT->Draw("EP");
  hRapidityMS->SetLineColor(kRed);
  hRapidityMS->SetMarkerColor(kRed);
  hRapidityMS->SetMarkerStyle(20);
  hRapidityMS->Draw("EP SAME");
  hRapidityMC->SetLineColor(kBlack);
  hRapidityMC->SetMarkerColor(kBlack);
  hRapidityMC->SetMarkerStyle(20);
  hRapidityMC->Draw("EP SAME");
  TLegend* legend = new TLegend(0.65, 0.5, 0.95, 0.9);
  legend->AddEntry(hRapidityVT, Form("VT tracks: #epsilon #times A = %.1f%%", effVT * 100), "lep");
  legend->AddEntry(hRapidityMS, Form("MS tracks: #epsilon #times A = %.1f%%", effMS * 100), "lep");
  legend->AddEntry(hRapidityMC, "Generated", "lep");
  legend->Draw();
  c1->SaveAs(Form("%s/single_muon_rapidity_comparison.png", dirSimu));
  
  setStyle(hEffVsRapidityVT);
  setStyle(hEffVsRapidityMS);
  bool firstEffDrawn = false;
  if (hEffVsRapidityVT) {
    hEffVsRapidityVT->SetLineColor(kBlue);
    hEffVsRapidityVT->SetMarkerColor(kBlue);
    hEffVsRapidityVT->Draw("AP");
    setStyleAxes(hEffVsRapidityVT);
    firstEffDrawn = true;
  }
  if (hEffVsRapidityMS) {
    hEffVsRapidityMS->SetLineColor(kRed);
    hEffVsRapidityMS->SetMarkerColor(kRed);
    hEffVsRapidityMS->Draw(firstEffDrawn ? "P SAME" : "AP");
    if (!firstEffDrawn) {
      setStyleAxes(hEffVsRapidityMS);
      firstEffDrawn = true;
    }
  }
  gPad->Modified();
  gPad->Update();
  TLegend* legend2 = new TLegend(0.2, 0.7, 0.38, 0.88);
  if (hEffVsRapidityVT)
    legend2->AddEntry(hEffVsRapidityVT, "VT", "lep");
  if (hEffVsRapidityMS)
    legend2->AddEntry(hEffVsRapidityMS, "MS", "lep");
  legend2->Draw();
  c1->SaveAs(Form("%s/single_muon_efficiency_comparison.png", dirSimu));

  fOut.Close();
}