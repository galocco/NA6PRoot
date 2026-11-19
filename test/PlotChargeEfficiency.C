#if !defined(__CINT__) || defined(__MAKECINT__)
#include <TFile.h>
#include <TTree.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TH3F.h>
#include <TCanvas.h>
#include <TPad.h>
#include <TLegend.h>
#include <TEfficiency.h>
#include <TParticle.h>
#include <TDatabasePDG.h>
#include <TLatex.h>
#include <TStyle.h>
#include <TMath.h>
#include <cmath>
#include <iostream>
#include <vector>
#include <algorithm>
#include "NA6PTrack.h"
#include "NA6PVerTelHit.h"
#include "NA6PLayoutParam.h"
#include "ConfigurableParam.h"
#endif

void setAxisStyle(TH1F* hist){
  hist->GetXaxis()->SetTitleSize(0.06);
  hist->GetYaxis()->SetTitleSize(0.06);
  hist->GetXaxis()->SetLabelSize(0.05);
  hist->GetYaxis()->SetLabelSize(0.05);
  hist->GetXaxis()->SetTitleOffset(1.20f);
  hist->GetYaxis()->SetTitleOffset(1.30f);
}

void setAxisStyle(TH2F* hist){
  hist->GetXaxis()->SetTitleSize(0.06);
  hist->GetYaxis()->SetTitleSize(0.06);
  hist->GetXaxis()->SetLabelSize(0.05);
  hist->GetYaxis()->SetLabelSize(0.05);
  hist->GetXaxis()->SetTitleOffset(1.20f);
  hist->GetYaxis()->SetTitleOffset(1.30f);
}

int GetTargetIndexFromLayout(float z,
                             const std::vector<float>& zTarget,
                             const std::vector<float>& halfThick)
{
  int best = 0;
  float bestDz = 1.e9f;
  for (size_t it = 0; it < zTarget.size(); ++it) {
    float dz = std::abs(z - zTarget[it]);
    if (dz <= halfThick[it] + 1.e-6f)
      return (int)it;
    if (dz < bestDz) {
      bestDz = dz;
      best = (int)it;
    }
  }
  return best;
}

void PlotChargeEfficiency(const char* dirSimu = ".",
                          const char* outPrefix = "ChargeEfficiency",
                          const char* layoutFile = "na6pLayout.ini",
                          int nPtBins = 60,
                          float ptMin = 0.f,
                          float ptMax = 3.f,
                          int nEtaBins = 100,
                          float etaMin = -3.f,
                          float etaMax = 7.f)
{
  std::cout << "\n=== PlotChargeEfficiency macro ===" << std::endl;
  std::cout << "Reading from directory: " << dirSimu << std::endl;
  // Open files
  const char* trackFile = Form("%s/TracksVerTel.root", dirSimu);
  const char* mcFile = Form("%s/MCKine.root", dirSimu);
  
  std::cout << "Opening " << trackFile << std::endl;
  std::cout << "Opening " << mcFile << std::endl;
  
  TFile* fTracks = TFile::Open(trackFile, "READ");
  if (!fTracks || fTracks->IsZombie()) {
    std::cerr << "Cannot open " << trackFile << std::endl;
    return;
  }

  TFile* fMC = TFile::Open(mcFile, "READ");
  if (!fMC || fMC->IsZombie()) {
    std::cerr << "Cannot open " << mcFile << std::endl;
    fTracks->Close();
    return;
  }

  std::cout << "Files opened successfully" << std::endl;
  std::cout << "Listing contents of TracksVerTel.root:" << std::endl;
  fTracks->ls();
  std::cout << "\nListing contents of MCKine.root:" << std::endl;
  fMC->ls();
  
  std::cout << "\nLooking for trees 'tracksVerTel' and 'mckine'" << std::endl;
  
  TTree* trTree = (TTree*)fTracks->Get("tracksVerTel");
  TTree* mcTree = (TTree*)fMC->Get("mckine");

  if (!trTree) {
    std::cerr << "Cannot find tree 'tracksVerTel' in " << trackFile << std::endl;
    fTracks->Close();
    fMC->Close();
    return;
  }
  
  if (!mcTree) {
    std::cerr << "Cannot find tree 'mckine' in " << mcFile << std::endl;
    fTracks->Close();
    fMC->Close();
    return;
  }

  std::vector<NA6PTrack>* trArr = nullptr;
  std::vector<TParticle>* mcArr = nullptr;
  trTree->SetBranchAddress("VerTel", &trArr);
  mcTree->SetBranchAddress("tracks", &mcArr);

  // Check if trees have entries
  int nEvTracks = trTree->GetEntries();
  int nEvMC = mcTree->GetEntries();
  
  std::cout << "TracksVerTel tree has " << nEvTracks << " entries" << std::endl;
  std::cout << "MCKine tree has " << nEvMC << " entries" << std::endl;
  
  if (nEvTracks == 0 || nEvMC == 0) {
    std::cerr << "Error: One or both trees are empty!" << std::endl;
    if (nEvTracks == 0) std::cerr << "  - TracksVerTel is empty" << std::endl;
    if (nEvMC == 0) std::cerr << "  - MCKine is empty" << std::endl;
    fTracks->Close();
    fMC->Close();
    return;
  }

  // Load layout
  na6p::conf::ConfigurableParam::updateFromFile(layoutFile, "", true);
  if (TGeoGlobalMagField::Instance()->GetField() == nullptr) {
    auto magField = new MagneticField();
    magField->loadField();
    magField->setAsGlobalField();
  }
  auto& layoutPar = NA6PLayoutParam::Instance();
  int nTargets = std::max(1, std::min((int)layoutPar.nTargets, (int)NA6PLayoutParam::MaxTargets));
  std::vector<float> zTarget(nTargets, 0.f), halfThick(nTargets, 0.05f);
  for (int it = 0; it < nTargets; ++it) {
    zTarget[it] = layoutPar.posTargetZ[it] + layoutPar.shiftTargets[2];
    halfThick[it] = 0.5f * std::max(1.e-4f, layoutPar.thicknessTarget[it]);
  }

  gStyle->SetOptStat(0);

  // Histograms for positive particles
  TH1F* hPtGenPos = new TH1F("hPtGenPos", ";#it{p}_{T} (GeV/#it{c});generated", nPtBins, ptMin, ptMax);
  TH1F* hPtRecPos = new TH1F("hPtRecPos", ";#it{p}_{T} (GeV/#it{c});reconstructed", nPtBins, ptMin, ptMax);
  TH1F* hEtaGenPos = new TH1F("hEtaGenPos", ";#it{#eta};generated", nEtaBins, etaMin, etaMax);
  TH1F* hEtaRecPos = new TH1F("hEtaRecPos", ";#it{#eta};reconstructed", nEtaBins, etaMin, etaMax);

  TH2F* h2GenPos = new TH2F("h2GenPos", ";#it{#eta};#it{p}_{T} (GeV/#it{c})", nEtaBins, etaMin, etaMax, nPtBins, ptMin, ptMax);
  TH2F* h2RecPos = new TH2F("h2RecPos", ";#it{#eta};#it{p}_{T} (GeV/#it{c})", nEtaBins, etaMin, etaMax, nPtBins, ptMin, ptMax);

  TH3F* h3GenPos = new TH3F("h3GenPos", ";#it{#eta};#it{p}_{T} (GeV/#it{c});target index", nEtaBins, etaMin, etaMax, nPtBins, ptMin, ptMax, nTargets, -0.5, nTargets - 0.5);
  TH3F* h3RecPos = new TH3F("h3RecPos", ";#it{#eta};#it{p}_{T} (GeV/#it{c});target index", nEtaBins, etaMin, etaMax, nPtBins, ptMin, ptMax, nTargets, -0.5, nTargets - 0.5);

  // Histograms for negative particles
  TH1F* hPtGenNeg = new TH1F("hPtGenNeg", ";#it{p}_{T} (GeV/#it{c});generated", nPtBins, ptMin, ptMax);
  TH1F* hPtRecNeg = new TH1F("hPtRecNeg", ";#it{p}_{T} (GeV/#it{c});reconstructed", nPtBins, ptMin, ptMax);
  TH1F* hEtaGenNeg = new TH1F("hEtaGenNeg", ";#it{#eta};generated", nEtaBins, etaMin, etaMax);
  TH1F* hEtaRecNeg = new TH1F("hEtaRecNeg", ";#it{#eta};reconstructed", nEtaBins, etaMin, etaMax);

  TH2F* h2GenNeg = new TH2F("h2GenNeg", ";#it{#eta};#it{p}_{T} (GeV/#it{c})", nEtaBins, etaMin, etaMax, nPtBins, ptMin, ptMax);
  TH2F* h2RecNeg = new TH2F("h2RecNeg", ";#it{#eta};#it{p}_{T} (GeV/#it{c})", nEtaBins, etaMin, etaMax, nPtBins, ptMin, ptMax);

  TH3F* h3GenNeg = new TH3F("h3GenNeg", ";#it{#eta};#it{p}_{T} (GeV/#it{c});target index", nEtaBins, etaMin, etaMax, nPtBins, ptMin, ptMax, nTargets, -0.5, nTargets - 0.5);
  TH3F* h3RecNeg = new TH3F("h3RecNeg", ";#it{#eta};#it{p}_{T} (GeV/#it{c});target index", nEtaBins, etaMin, etaMax, nPtBins, ptMin, ptMax, nTargets, -0.5, nTargets - 0.5);

  // Impact parameter histograms: primary vs secondary
  const int nMomBins = 20;
  const double maxP = 10.;
  const double impMax = 1000.;  // in microns (expanded range)
  const int nImpBins = 200;  // bins for 1D histograms
  TH2F* hImpParXVsPPrimary = new TH2F("hImpParXVsPPrimary", ";p (GeV/#it{c});Track Imp. Par. X (#mum);counts", nMomBins, 0., maxP, nImpBins, -impMax, impMax);
  TH2F* hImpParYVsPPrimary = new TH2F("hImpParYVsPPrimary", ";p (GeV/#it{c});Track Imp. Par. Y (#mum);counts", nMomBins, 0., maxP, nImpBins, -impMax, impMax);
  TH2F* hImpParXVsPSecondary = new TH2F("hImpParXVsPSecondary", ";p (GeV/#it{c});Track Imp. Par. X (#mum);counts", nMomBins, 0., maxP, nImpBins, -impMax, impMax);
  TH2F* hImpParYVsPSecondary = new TH2F("hImpParYVsPSecondary", ";p (GeV/#it{c});Track Imp. Par. Y (#mum);counts", nMomBins, 0., maxP, nImpBins, -impMax, impMax);
  
  // 1D impact parameter histograms
  TH1F* hImpParXPrimary = new TH1F("hImpParXPrimary", ";Track Imp. Par. X (#mum);counts", nImpBins, -impMax, impMax);
  TH1F* hImpParYPrimary = new TH1F("hImpParYPrimary", ";Track Imp. Par. Y (#mum);counts", nImpBins, -impMax, impMax);
  TH1F* hImpParXSecondary = new TH1F("hImpParXSecondary", ";Track Imp. Par. X (#mum);counts", nImpBins, -impMax, impMax);
  TH1F* hImpParYSecondary = new TH1F("hImpParYSecondary", ";Track Imp. Par. Y (#mum);counts", nImpBins, -impMax, impMax);
  TH1F* hDxySignedPrimary = new TH1F("hDxySignedPrimary", ";#it{d}_{xy} (#mum);Counts", nImpBins, -impMax, impMax);
  TH1F* hDxySignedSecondary = new TH1F("hDxySignedSecondary", ";#it{d}_{xy} (#mum);Counts", nImpBins, -impMax, impMax);

  // Per-target histograms
  std::vector<TH1F*> hPtGenPosT(nTargets), hPtRecPosT(nTargets);
  std::vector<TH1F*> hEtaGenPosT(nTargets), hEtaRecPosT(nTargets);
  std::vector<TH2F*> h2GenPosT(nTargets), h2RecPosT(nTargets);
  std::vector<TH1F*> hPtGenNegT(nTargets), hPtRecNegT(nTargets);
  std::vector<TH1F*> hEtaGenNegT(nTargets), hEtaRecNegT(nTargets);
  std::vector<TH2F*> h2GenNegT(nTargets), h2RecNegT(nTargets);
  std::vector<double> nGenPosT(nTargets, 0.), nRecPosT(nTargets, 0.);
  std::vector<double> nGenNegT(nTargets, 0.), nRecNegT(nTargets, 0.);

  for (int it = 0; it < nTargets; ++it) {
    hPtGenPosT[it] = new TH1F(Form("hPtGenPos_t%d", it), ";#it{p}_{T} (GeV/#it{c});generated", nPtBins, ptMin, ptMax);
    hPtRecPosT[it] = new TH1F(Form("hPtRecPos_t%d", it), ";#it{p}_{T} (GeV/#it{c});reconstructed", nPtBins, ptMin, ptMax);
    hEtaGenPosT[it] = new TH1F(Form("hEtaGenPos_t%d", it), ";#it{#eta};generated", nEtaBins, etaMin, etaMax);
    hEtaRecPosT[it] = new TH1F(Form("hEtaRecPos_t%d", it), ";#it{#eta};reconstructed", nEtaBins, etaMin, etaMax);
    h2GenPosT[it] = new TH2F(Form("h2GenPos_t%d", it), ";#it{#eta};#it{p}_{T} (GeV/#it{c})", nEtaBins, etaMin, etaMax, nPtBins, ptMin, ptMax);
    h2RecPosT[it] = new TH2F(Form("h2RecPos_t%d", it), ";#it{#eta};#it{p}_{T} (GeV/#it{c})", nEtaBins, etaMin, etaMax, nPtBins, ptMin, ptMax);

    hPtGenNegT[it] = new TH1F(Form("hPtGenNeg_t%d", it), ";#it{p}_{T} (GeV/#it{c});generated", nPtBins, ptMin, ptMax);
    hPtRecNegT[it] = new TH1F(Form("hPtRecNeg_t%d", it), ";#it{p}_{T} (GeV/#it{c});reconstructed", nPtBins, ptMin, ptMax);
    hEtaGenNegT[it] = new TH1F(Form("hEtaGenNeg_t%d", it), ";#it{#eta};generated", nEtaBins, etaMin, etaMax);
    hEtaRecNegT[it] = new TH1F(Form("hEtaRecNeg_t%d", it), ";#it{#eta};reconstructed", nEtaBins, etaMin, etaMax);
    h2GenNegT[it] = new TH2F(Form("h2GenNeg_t%d", it), ";#it{#eta};#it{p}_{T} (GeV/#it{c})", nEtaBins, etaMin, etaMax, nPtBins, ptMin, ptMax);
    h2RecNegT[it] = new TH2F(Form("h2RecNeg_t%d", it), ";#it{#eta};#it{p}_{T} (GeV/#it{c})", nEtaBins, etaMin, etaMax, nPtBins, ptMin, ptMax);
  }

  // Event loop
  int nEv = trTree->GetEntries();
  std::cout << "Processing " << nEv << " events..." << std::endl;

  for (int iEv = 0; iEv < nEv; ++iEv) {
    if (iEv % 100 == 0) {
      std::cout << "Event " << iEv << "/" << nEv << std::endl;
    }
    mcTree->GetEntry(iEv);
    trTree->GetEntry(iEv);

    if (!mcArr || !trArr) {
      continue;
    }

    // Get primary vertex
    double xvert = 0, yvert = 0, zvert = 0;
    for (int jp = 0; jp < (int)mcArr->size(); ++jp) {
      auto curPart = mcArr->at(jp);
      if (curPart.IsPrimary()) {
        xvert = curPart.Vx();
        yvert = curPart.Vy();
        zvert = curPart.Vz();
        break;
      }
    }

    // Fill generated particles
    for (int jp = 0; jp < (int)mcArr->size(); ++jp) {
      auto curPart = mcArr->at(jp);
      if (!curPart.IsPrimary()) {
        continue;
      }

      int pdg = curPart.GetPdgCode();
      int charge = 0;
      // Get charge from PDG code
      TParticlePDG* pPDG = TDatabasePDG::Instance()->GetParticle(pdg);
      if (pPDG) {
        charge = (int)(pPDG->Charge() / 3); // Charge in units of e/3
      }

      double px = curPart.Px();
      double py = curPart.Py();
      double pz = curPart.Pz();
      double pt = std::sqrt(px * px + py * py);
      double p = curPart.P();
      double theta = std::acos(pz / p);
      double eta = -std::log(std::tan(theta / 2.));

      if (pt < ptMin || pt >= ptMax || eta < etaMin || eta >= etaMax) {
        continue;
      }

      int targetIdx = GetTargetIndexFromLayout(zvert, zTarget, halfThick);

      if (charge > 0) {
        hPtGenPos->Fill(pt);
        hEtaGenPos->Fill(eta);
        h2GenPos->Fill(eta, pt);
        h3GenPos->Fill(eta, pt, targetIdx);
        hPtGenPosT[targetIdx]->Fill(pt);
        hEtaGenPosT[targetIdx]->Fill(eta);
        h2GenPosT[targetIdx]->Fill(eta, pt);
        nGenPosT[targetIdx] += 1.;
      } else if (charge < 0) {
        hPtGenNeg->Fill(pt);
        hEtaGenNeg->Fill(eta);
        h2GenNeg->Fill(eta, pt);
        h3GenNeg->Fill(eta, pt, targetIdx);
        hPtGenNegT[targetIdx]->Fill(pt);
        hEtaGenNegT[targetIdx]->Fill(eta);
        h2GenNegT[targetIdx]->Fill(eta, pt);
        nGenNegT[targetIdx] += 1.;
      }
    }

    // Fill reconstructed tracks
    
    for (int jTr = 0; jTr < (int)trArr->size(); ++jTr) {
      NA6PTrack tr = trArr->at(jTr);
      int nClusters = tr.getNVTHits();
      if (nClusters < 5) {
        continue;
      }

      double pxyzReco[3];
      tr.getPXYZ(pxyzReco);
      double ptReco = std::sqrt(pxyzReco[0] * pxyzReco[0] + pxyzReco[1] * pxyzReco[1]);
      double pReco = tr.getP();
      double thetaReco = std::acos(pxyzReco[2] / pReco);
      double etaReco = -std::log(std::tan(thetaReco / 2.));
      
      // Calculate impact parameters at the primary vertex z position
      tr.propagateToZBxByBz(zvert);
      double impparX = tr.getXLab() - xvert;  // Impact parameter in X: distance from vertex X
      double impparY = tr.getYLab() - yvert;  // Impact parameter in Y: distance from vertex Y
      double dxySigned = 0.;
      if (ptReco > 0.) {
        dxySigned = (impparX * pxyzReco[1] - impparY * pxyzReco[0]) / ptReco;
      }

      int mcLabel = tr.getParticleID();
      if (mcLabel < 0) {
        continue; // Skip fake tracks
      }

      auto mcPart = mcArr->at(mcLabel);
      
      // Fill impact parameter histograms for all tracks
      if (mcPart.IsPrimary()) {
        hImpParXVsPPrimary->Fill(pReco, impparX * 1e4);
        hImpParYVsPPrimary->Fill(pReco, impparY * 1e4);
        hImpParXPrimary->Fill(impparX * 1e4);
        hImpParYPrimary->Fill(impparY * 1e4);
        hDxySignedPrimary->Fill(dxySigned * 1e4);
      } else {
        hImpParXVsPSecondary->Fill(pReco, impparX * 1e4);
        hImpParYVsPSecondary->Fill(pReco, impparY * 1e4);
        hImpParXSecondary->Fill(impparX * 1e4);
        hImpParYSecondary->Fill(impparY * 1e4);
        hDxySignedSecondary->Fill(dxySigned * 1e4);
      }

      // Efficiency calculation only for primary particles
      if (!mcPart.IsPrimary()) {
        continue;
      }


      int pdg = mcPart.GetPdgCode();
      int charge = 0;
      TParticlePDG* pPDG = TDatabasePDG::Instance()->GetParticle(pdg);
      if (pPDG) {
        charge = (int)(pPDG->Charge() / 3);
      }

      int targetIdx = GetTargetIndexFromLayout(zvert, zTarget, halfThick);
      float ptGen = mcPart.Pt();
      float etaGen = mcPart.Eta();

      if (ptGen < ptMin || ptGen >= ptMax || etaGen < etaMin || etaGen >= etaMax) {
        continue;
      }
      if (charge > 0) {
        hPtRecPos->Fill(ptGen);
        hEtaRecPos->Fill(etaGen);
        h2RecPos->Fill(etaGen, ptGen);
        h3RecPos->Fill(etaGen, ptGen, targetIdx);
        hPtRecPosT[targetIdx]->Fill(ptGen);
        hEtaRecPosT[targetIdx]->Fill(etaGen);
        h2RecPosT[targetIdx]->Fill(etaGen, ptGen);
        nRecPosT[targetIdx] += 1.;
      } else if (charge < 0) {
        hPtRecNeg->Fill(ptGen);
        hEtaRecNeg->Fill(etaGen);
        h2RecNeg->Fill(etaGen, ptGen);
        h3RecNeg->Fill(etaGen, ptGen, targetIdx);
        hPtRecNegT[targetIdx]->Fill(ptGen);
        hEtaRecNegT[targetIdx]->Fill(etaGen);
        h2RecNegT[targetIdx]->Fill(etaGen, ptGen);
        nRecNegT[targetIdx] += 1.;
      }
    }
  }

  // Create efficiency histograms
  TEfficiency* hEffPtPos = new TEfficiency(*hPtRecPos, *hPtGenPos);
  TEfficiency* hEffEtaPos = new TEfficiency(*hEtaRecPos, *hEtaGenPos);
  TEfficiency* hEffPtNeg = new TEfficiency(*hPtRecNeg, *hPtGenNeg);
  TEfficiency* hEffEtaNeg = new TEfficiency(*hEtaRecNeg, *hEtaGenNeg);
  std::cout << "hPtGenNeg entries: " << hPtGenNeg->GetEntries() << std::endl;
  std::cout << "hPtRecNeg entries: " << hPtRecNeg->GetEntries() << std::endl;
  // 2D efficiency
  TH2F* hEff2DPos = (TH2F*)h2RecPos->Clone("hEff2DPos");
  hEff2DPos->SetTitle("Efficiency x Acceptance (Positive);#it{#eta};#it{p}_{T} (GeV/#it{c})");
  hEff2DPos->Divide(h2RecPos, h2GenPos, 1., 1., "B");

  TH2F* hEff2DNeg = (TH2F*)h2RecNeg->Clone("hEff2DNeg");
  hEff2DNeg->SetTitle("Efficiency x Acceptance (Negative);#it{#eta};#it{p}_{T} (GeV/#it{c})");
  hEff2DNeg->Divide(h2RecNeg, h2GenNeg, 1., 1., "B");

  // 3D efficiency
  TH3F* hEff3DPos = (TH3F*)h3RecPos->Clone("hEff3DPos");
  hEff3DPos->SetTitle("Efficiency x Acceptance (Positive);#it{#eta};#it{p}_{T} (GeV/#it{c});target index");
  hEff3DPos->Divide(h3RecPos, h3GenPos, 1., 1., "B");

  TH3F* hEff3DNeg = (TH3F*)h3RecNeg->Clone("hEff3DNeg");
  hEff3DNeg->SetTitle("Efficiency x Acceptance (Negative);#it{#eta};#it{p}_{T} (GeV/#it{c});target index");
  hEff3DNeg->Divide(h3RecNeg, h3GenNeg, 1., 1., "B");

  double effPosGlobal = (hPtGenPos->GetEntries() > 0.) ? (hPtRecPos->GetEntries() / hPtGenPos->GetEntries()) : 0.;
  double effNegGlobal = (hPtGenNeg->GetEntries() > 0.) ? (hPtRecNeg->GetEntries() / hPtGenNeg->GetEntries()) : 0.;

  // Plot 1D efficiencies
  TCanvas* c1D = new TCanvas("c1D", "Efficiency 1D", 1600, 700);
  c1D->Divide(2, 1);

  auto drawEff = [](TEfficiency* eff, Color_t color, Style_t markerStyle) {
    eff->SetLineColor(color);
    eff->SetMarkerColor(color);
    eff->SetMarkerStyle(markerStyle);
    eff->SetMarkerSize(1.25);
    eff->Draw("AP");
    gPad->Update();
    if (auto* gr = eff->GetPaintedGraph()) {
      gr->SetMinimum(0.0);
      gr->SetMaximum(1);
      gr->GetXaxis()->SetTitleSize(0.060f);
      gr->GetYaxis()->SetTitleSize(0.060f);
      gr->GetXaxis()->SetLabelSize(0.050f);
      gr->GetYaxis()->SetLabelSize(0.050f);
      gr->GetXaxis()->SetTitleOffset(1.15f);
      gr->GetYaxis()->SetTitleOffset(1.30f);
      gr->GetXaxis()->SetNdivisions(510);
      gr->GetYaxis()->SetNdivisions(510);
    }
    gPad->Modified();
    gPad->Update();
  };

  c1D->cd(1);
  gPad->SetLeftMargin(0.18);
  gPad->SetBottomMargin(0.14);
  hEffPtPos->SetTitle("Efficiency vs #it{p}_{T} (Positive);#it{p}_{T} (GeV/#it{c});Efficiency");
  drawEff(hEffPtPos, kRed + 1, 20);

  c1D->cd(2);
  gPad->SetLeftMargin(0.18);
  gPad->SetBottomMargin(0.14);
  hEffPtNeg->SetTitle("Efficiency vs #it{p}_{T} (Negative);#it{p}_{T} (GeV/#it{c});Efficiency");
  drawEff(hEffPtNeg, kBlue + 1, 21);

  c1D->SaveAs(Form("%s-1D-Pt.png", outPrefix));

  TCanvas* c1DEta = new TCanvas("c1DEta", "Efficiency 1D Eta", 1600, 700);
  c1DEta->Divide(2, 1);

  c1DEta->cd(1);
  gPad->SetLeftMargin(0.18);
  gPad->SetBottomMargin(0.14);
  hEffEtaPos->SetTitle("Efficiency vs #it{#eta} (Positive);#it{#eta};Efficiency");
  drawEff(hEffEtaPos, kRed + 1, 20);

  c1DEta->cd(2);
  gPad->SetLeftMargin(0.18);
  gPad->SetBottomMargin(0.14);
  hEffEtaNeg->SetTitle("Efficiency vs #it{#eta} (Negative);#it{#eta};Efficiency");
  drawEff(hEffEtaNeg, kBlue + 1, 21);

  c1DEta->SaveAs(Form("%s-1D-Eta.png", outPrefix));

  // Plot 2D efficiencies
  TCanvas* c2D = new TCanvas("c2D", "Efficiency 2D", 1600, 700);
  c2D->Divide(2, 1);

  c2D->cd(1);
  gPad->SetLeftMargin(0.16);
  gPad->SetRightMargin(0.12);
  gPad->SetBottomMargin(0.18);
  gPad->SetTopMargin(0.07);
  gStyle->SetPaintTextFormat(".2f");
  gStyle->SetPalette(kBird);
  hEff2DPos->SetMinimum(0.);
  hEff2DPos->SetMaximum(1.);
  hEff2DPos->GetXaxis()->SetRangeUser(1, 5);
  setAxisStyle(hEff2DPos);
  hEff2DPos->GetXaxis()->SetNdivisions(510);
  hEff2DPos->GetYaxis()->SetNdivisions(510);
  hEff2DPos->Draw("colz");

  c2D->cd(2);
  gPad->SetLeftMargin(0.16);
  gPad->SetRightMargin(0.12);
  gPad->SetBottomMargin(0.18);
  gPad->SetTopMargin(0.07);
  gStyle->SetPaintTextFormat(".2f");
  gStyle->SetPalette(kBird);
  hEff2DNeg->SetMinimum(0.);
  hEff2DNeg->SetMaximum(1.);
  hEff2DNeg->GetXaxis()->SetRangeUser(1, 5);
  setAxisStyle(hEff2DNeg);
  hEff2DNeg->GetXaxis()->SetNdivisions(510);
  hEff2DNeg->GetYaxis()->SetNdivisions(510);
  hEff2DNeg->Draw("colz");

  c2D->SaveAs(Form("%s-2D.png", outPrefix));

  // Impact parameter canvas
  TCanvas* cImpPar = new TCanvas("cImpPar", "Impact Parameters Primary vs Secondary", 1200, 500);
  cImpPar->Divide(2, 1);

  cImpPar->cd(1);
  gPad->SetLeftMargin(0.14);
  gPad->SetRightMargin(0.16);
  gPad->SetBottomMargin(0.14);
  gPad->SetTopMargin(0.07);
  hImpParXVsPPrimary->SetTitle("Impact Parameter X - Primary;p (GeV/#it{c});Imp. Par. X (#mum)");
  hImpParXVsPPrimary->Draw("colz");

  cImpPar->cd(2);
  gPad->SetLeftMargin(0.14);
  gPad->SetRightMargin(0.16);
  gPad->SetBottomMargin(0.14);
  gPad->SetTopMargin(0.07);
  hImpParXVsPSecondary->SetTitle("Impact Parameter X - Secondary;p (GeV/#it{c});Imp. Par. X (#mum)");
  hImpParXVsPSecondary->Draw("colz");

  cImpPar->SaveAs(Form("%s-ImpParX.png", outPrefix));

  TCanvas* cImpParY = new TCanvas("cImpParY", "Impact Parameters Y", 1200, 500);
  cImpParY->Divide(2, 1);

  cImpParY->cd(1);
  gPad->SetLeftMargin(0.14);
  gPad->SetRightMargin(0.16);
  gPad->SetBottomMargin(0.14);
  gPad->SetTopMargin(0.07);
  hImpParYVsPPrimary->SetTitle("Impact Parameter Y - Primary;p (GeV/#it{c});Imp. Par. Y (#mum)");
  hImpParYVsPPrimary->Draw("colz");

  cImpParY->cd(2);
  gPad->SetLeftMargin(0.14);
  gPad->SetRightMargin(0.16);
  gPad->SetBottomMargin(0.14);
  gPad->SetTopMargin(0.07);
  hImpParYVsPSecondary->SetTitle("Impact Parameter Y - Secondary;p (GeV/#it{c});Imp. Par. Y (#mum)");
  hImpParYVsPSecondary->Draw("colz");

  cImpParY->SaveAs(Form("%s-ImpParY.png", outPrefix));

  TCanvas* cGenPtEta = new TCanvas("cGenPtEta", "", 1600, 700);
  cGenPtEta->Divide(2, 1);

  // Pad sinistro
  cGenPtEta->cd(1);
  gPad->SetLeftMargin(0.18);
  gPad->SetBottomMargin(0.18);
  gPad->SetTopMargin(0.07);
  gPad->SetRightMargin(0.10);

  // Pad destro
  cGenPtEta->cd(2);
  gPad->SetLeftMargin(0.18);
  gPad->SetBottomMargin(0.18);
  gPad->SetTopMargin(0.07);
  gPad->SetRightMargin(0.10);



  cGenPtEta->cd(1);
  hPtGenPos->GetYaxis()->SetTitle("#frac{dN}{d#it{p}_{T}} (GeV/#it{c})^{-1}");

  hPtGenPos->Scale(1./nEv);
  hPtGenNeg->Scale(1./nEv);
  setAxisStyle(hPtGenPos);
  hPtGenPos->SetLineColor(kRed + 1);
  hPtGenPos->SetLineWidth(2);
  hPtGenPos->Draw("hist");
  hPtGenNeg->SetLineColor(kBlue + 1);
  hPtGenNeg->SetLineWidth(2);
  hPtGenNeg->Draw("hist same");
  TLegend *legGenPt = new TLegend(0.65, 0.75, 0.88, 0.88);
  legGenPt->SetBorderSize(0);
  legGenPt->SetFillStyle(0);
  legGenPt->AddEntry(hPtGenPos, "Positive", "l");
  legGenPt->AddEntry(hPtGenNeg, "Negative", "l");
  legGenPt->Draw();
  cGenPtEta->cd(2);
  hEtaGenPos->Scale(1./nEv);
  hEtaGenNeg->Scale(1./nEv);
  hEtaGenPos->GetYaxis()->SetTitle("#frac{dN}{d#it{#eta}}");
  setAxisStyle(hEtaGenPos);
  hEtaGenPos->SetLineColor(kRed + 1);
  hEtaGenPos->SetLineWidth(2);
  hEtaGenPos->Draw("hist");
  hEtaGenNeg->SetLineColor(kBlue + 1);
  hEtaGenNeg->SetLineWidth(2);
  hEtaGenNeg->Draw("hist same");
  cGenPtEta->SaveAs(Form("%s-GenPtEtaPos.png", outPrefix));

  TCanvas* cDxySigned = new TCanvas("cDxySigned", "Signed dxy Primary vs Secondary", 800, 700);

  gPad->SetLeftMargin(0.18);
  gPad->SetBottomMargin(0.18);
  gPad->SetTopMargin(0.07);
  gPad->SetRightMargin(0.10);
  setAxisStyle(hDxySignedPrimary);
  setAxisStyle(hDxySignedSecondary);

  hDxySignedPrimary->Scale(1./nEv);
  hDxySignedPrimary->GetXaxis()->SetNdivisions(805);
  hDxySignedSecondary->Scale(1./nEv);
  hDxySignedPrimary->Draw("hist");
  hDxySignedPrimary->SetLineColor(kRed + 1);
  hDxySignedPrimary->SetLineWidth(2);
  hDxySignedSecondary->Draw("hist same");
  hDxySignedSecondary->SetLineColor(kBlue + 1);
  hDxySignedSecondary->SetLineWidth(2);
  TLegend *legDxy = new TLegend(0.65, 0.75, 0.88, 0.88);
  legDxy->SetBorderSize(0);
  legDxy->SetFillStyle(0);
  legDxy->AddEntry(hDxySignedPrimary, "Primary", "l");
  legDxy->AddEntry(hDxySignedSecondary, "Secondary", "l");
  legDxy->Draw();
  cDxySigned->SetLogy();
  cDxySigned->SaveAs(Form("%s-DxySigned.png", outPrefix));

  // Save to ROOT file
  TFile* fout = new TFile(Form("%s.root", outPrefix), "RECREATE");
  hPtGenPos->Write();
  hPtRecPos->Write();
  hEtaGenPos->Write();
  hEtaRecPos->Write();
  hPtGenNeg->Write();
  hPtRecNeg->Write();
  hEtaGenNeg->Write();
  hEtaRecNeg->Write();
  h2GenPos->Write();
  h2RecPos->Write();
  h2GenNeg->Write();
  h2RecNeg->Write();
  h3GenPos->Write();
  h3RecPos->Write();
  h3GenNeg->Write();
  h3RecNeg->Write();
  hEffPtPos->Write();
  hEffEtaPos->Write();
  hEffPtNeg->Write();
  hEffEtaNeg->Write();
  hEff2DPos->Write();
  hEff2DNeg->Write();
  hEff3DPos->Write();
  hEff3DNeg->Write();
  hImpParXVsPPrimary->Write();
  hImpParYVsPPrimary->Write();
  hImpParXVsPSecondary->Write();
  hImpParYVsPSecondary->Write();
  hImpParXPrimary->Write();
  hImpParYPrimary->Write();
  hImpParXSecondary->Write();
  hImpParYSecondary->Write();
  hDxySignedPrimary->Write();
  hDxySignedSecondary->Write();
  c1D->Write();
  c1DEta->Write();
  c2D->Write();
  cImpPar->Write();
  cImpParY->Write();
  cDxySigned->Write();

  for (int it = 0; it < nTargets; ++it) {
    hPtGenPosT[it]->Write();
    hPtRecPosT[it]->Write();
    hEtaGenPosT[it]->Write();
    hEtaRecPosT[it]->Write();
    h2GenPosT[it]->Write();
    h2RecPosT[it]->Write();
    hPtGenNegT[it]->Write();
    hPtRecNegT[it]->Write();
    hEtaGenNegT[it]->Write();
    hEtaRecNegT[it]->Write();
    h2GenNegT[it]->Write();
    h2RecNegT[it]->Write();
  }

  fout->Close();

  // Summary
  std::cout << "\n=== Charge Efficiency Summary ===\n";
  std::cout << "Input directory: " << dirSimu << "\n";
  std::cout << "Layout file: " << layoutFile << "\n";
  std::cout << "\nPositive particles:\n";
  std::cout << "  Generated: " << hPtGenPos->GetEntries() << "\n";
  std::cout << "  Reconstructed: " << hPtRecPos->GetEntries() << "\n";
  std::cout << "  Global efficiency: " << effPosGlobal << "\n";
  std::cout << "\nNegative particles:\n";
  std::cout << "  Generated: " << hPtGenNeg->GetEntries() << "\n";
  std::cout << "  Reconstructed: " << hPtRecNeg->GetEntries() << "\n";
  std::cout << "  Global efficiency: " << effNegGlobal << "\n";
  std::cout << "\nTargets: " << nTargets << "\n";
  std::cout << "Output: " << outPrefix << ".root and " << outPrefix << "-*.png\n\n";

  fTracks->Close();
  fMC->Close();
}
