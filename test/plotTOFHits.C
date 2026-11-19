#if !defined(__CINT__) || defined(__MAKECINT__)
#include <TTree.h>
#include <TFile.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TH3D.h>
#include <TEfficiency.h>
#include <TParticle.h>
#include <TCanvas.h>
#include <TGraph.h>
#include <TStyle.h>
#include <TDatabasePDG.h>
#include <TMath.h>
#include <cmath>
#include <algorithm>
#include <unordered_map>
#include "NA6PTOFCluster.h"
#include "NA6PVerTelCluster.h"
#endif

namespace ConstHits {
  constexpr double cLight      = 29.9792458;
  constexpr double mPion       = 0.13957;
  constexpr double mKaon       = 0.49368;
  constexpr double nsPerSec    = 1.e9;

  constexpr int    nTargets    = 5;
  constexpr double targetZ[nTargets] = {-5.3, -4.1, -2.9, -1.7, -0.5};
  constexpr int minVerTelHits = 4;
}

TTree* openFileAndTreeHits(const char* path, const char* treeName)
{
  TFile* f = TFile::Open(path);
  if (!f || f->IsZombie()) {
    printf("ERROR: Cannot open %s\n", path);
    return nullptr;
  }
  TTree* t = dynamic_cast<TTree*>(f->Get(treeName));
  if (!t) {
    printf("ERROR: Cannot find tree '%s' in %s\n", treeName, path);
    return nullptr;
  }
  return t;
}

void drawBetaCurveHits(double mass, int color)
{
  constexpr int np = 200;
  double xp[np], yp[np];
  for (int i = 0; i < np; ++i) {
    xp[i] = 0.1 + i * 19.9 / np;
    yp[i] = xp[i] / std::sqrt(xp[i] * xp[i] + mass * mass);
  }
  auto* gr = new TGraph(np, xp, yp);
  gr->SetLineColor(color);
  gr->SetLineWidth(2);
  gr->SetLineStyle(2);
  gr->Draw("L same");
}

TEfficiency* makeEfficiencyHits(TH1* hNum, TH1* hDen,
                                const char* name, const char* title)
{
  for (int b = 1; b <= hNum->GetNbinsX(); ++b) {
    if (hNum->GetBinContent(b) > hDen->GetBinContent(b)) {
      hNum->SetBinContent(b, hDen->GetBinContent(b));
    }
  }
  if (!TEfficiency::CheckConsistency(*hNum, *hDen)) {
    printf("  WARNING [%s]: TEfficiency consistency check failed — skipping.\n", name);
    return nullptr;
  }
  auto* eff = new TEfficiency(*hNum, *hDen);
  eff->SetName(name);
  eff->SetTitle(title);
  return eff;
}

void plotTOFHits(const char* dirSimu = ".", float pMaxIDProton = 3)
{
  gStyle->SetOptStat(0);
  gStyle->SetTitleSize(0.05, "XY");
  gStyle->SetLabelSize(0.045, "XY");
  gStyle->SetTitleOffset(1.1, "X");
  gStyle->SetTitleOffset(1.0, "Y");

  auto* tMC     = openFileAndTreeHits(Form("%s/MCKine.root",         dirSimu), "mckine");
  auto* tClus   = openFileAndTreeHits(Form("%s/ClustersTOF.root",    dirSimu), "clustersTOF");
  auto* tVerTel = openFileAndTreeHits(Form("%s/ClustersVerTel.root", dirSimu), "clustersVerTel");
  if (!tMC || !tClus || !tVerTel) return;

  std::vector<TParticle>* mcArr = nullptr;
  std::vector<NA6PTOFCluster> tofClus, *tofClusPtr = &tofClus;
  std::vector<NA6PVerTelCluster> verTelClus, *verTelClusPtr = &verTelClus;

  tMC->SetBranchAddress("tracks", &mcArr);
  tClus->SetBranchAddress("TOF", &tofClusPtr);
  tVerTel->SetBranchAddress("VerTel", &verTelClusPtr);

  const int nEv = std::min({(int)tMC->GetEntries(), (int)tClus->GetEntries(), (int)tVerTel->GetEntries()});
  printf("Processing %d events (hits-only)\n", nEv);

  auto* hBetaVsP = new TH2F("hBetaVsP", ";#it{p} (GeV/#it{c});#beta", 50, 0., 5., 300, 0.4, 1.10);
  auto* hBetaVsPprotons = new TH2F("hBetaVsPprotons", ";#it{p} (GeV/#it{c});#beta", 50, 0., 5., 300, 0.4, 1.10);

  TH2F* hBetaVsPtgt[ConstHits::nTargets];
  for (int it = 0; it < ConstHits::nTargets; ++it) {
    hBetaVsPtgt[it] = new TH2F(
      Form("hBetaVsP_tgt%d", it),
      Form("target %d (z=%.1f cm);#it{p} (GeV/#it{c});#beta", it, ConstHits::targetZ[it]),
      50, 0., 5., 300, 0.4, 1.10);
  }

  auto* hGenPVsYVsPt = new TH3D("hGenPVsYVsPt", ";p (GeV/c);y;p_{T} (GeV/c)", 50, 0., 50., 50, -2., 6., 50, 0., 5.);

  const int nRapBins = 100;
  const double rapMin = -3., rapMax = 7.;
  const int nEtaBins = 100;
  const double etaMin = -3., etaMax = 7.;
  const int nPBins = 100;
  const double pMinH = 0., pMaxH = 50.;
  const int nPtBins = 100;
  const double ptMin = 0., ptMax = 5.;

  auto* hProtonGenVsY = new TH1F("hProtonGenVsY", ";y;counts", nRapBins, rapMin, rapMax);
  auto* hProtonAccVsY = new TH1F("hProtonAccVsY", ";y;counts", nRapBins, rapMin, rapMax);
  auto* hProtonWithTOFVsY = new TH1F("hProtonWithTOFVsY", ";y;counts", nRapBins, rapMin, rapMax);
  auto* hProtonGenIDVsY = new TH1F("hProtonGenIDVsY", ";y;counts", nRapBins, rapMin, rapMax);
  auto* hProtonWithTOFIDVsY = new TH1F("hProtonWithTOFIDVsY", ";y;counts", nRapBins, rapMin, rapMax);

  auto* hProtonGenVsEta = new TH1F("hProtonGenVsEta", ";#eta;counts", nEtaBins, etaMin, etaMax);
  auto* hProtonAccVsEta = new TH1F("hProtonAccVsEta", ";#eta;counts", nEtaBins, etaMin, etaMax);
  auto* hProtonWithTOFVsEta = new TH1F("hProtonWithTOFVsEta", ";#eta;counts", nEtaBins, etaMin, etaMax);
  auto* hProtonGenIDVsEta = new TH1F("hProtonGenIDVsEta", ";#eta;counts", nEtaBins, etaMin, etaMax);
  auto* hProtonWithTOFIDVsEta = new TH1F("hProtonWithTOFIDVsEta", ";#eta;counts", nEtaBins, etaMin, etaMax);

  auto* hProtonGenPVsY = new TH2F("hProtonGenPVsY", ";p (GeV/c);y;counts", nPBins, pMinH, pMaxH, nRapBins, rapMin, rapMax);
  auto* hProtonAccPVsY = new TH2F("hProtonAccPVsY", ";p (GeV/c);y;counts", nPBins, pMinH, pMaxH, nRapBins, rapMin, rapMax);
  auto* hProtonRecPVsY = new TH2F("hProtonRecPVsY", ";p (GeV/c);y;counts", nPBins, pMinH, pMaxH, nRapBins, rapMin, rapMax);

  auto* hProtonGenPtVsEta = new TH2F("hProtonGenPtVsEta", ";p_{T} (GeV/c);#eta;counts", nPtBins, ptMin, ptMax, nEtaBins, etaMin, etaMax);
  auto* hProtonRecPtVsEta = new TH2F("hProtonRecPtVsEta", ";p_{T} (GeV/c);#eta;counts", nPtBins, ptMin, ptMax, nEtaBins, etaMin, etaMax);
  auto* hProtonAccPtVsEta = new TH2F("hProtonAccPtVsEta", ";p_{T} (GeV/c);#eta;counts", nPtBins, ptMin, ptMax, nEtaBins, etaMin, etaMax);

  auto* hProtonGenVsPt = new TH1F("hProtonGenVsPt", ";p_{T} (GeV/c);counts", nPtBins, ptMin, ptMax);
  auto* hProtonAccVsPt = new TH1F("hProtonAccVsPt", ";p_{T} (GeV/c);counts", nPtBins, ptMin, ptMax);
  auto* hProtonWithTOFVsPt = new TH1F("hProtonWithTOFVsPt", ";p_{T} (GeV/c);counts", nPtBins, ptMin, ptMax);
  auto* hProtonGenIDVsPt = new TH1F("hProtonGenIDVsPt", ";p_{T} (GeV/c);counts", nPtBins, ptMin, ptMax);
  auto* hProtonWithTOFIDVsPt = new TH1F("hProtonWithTOFIDVsPt", ";p_{T} (GeV/c);counts", nPtBins, ptMin, ptMax);

  auto* hProtonGenVsP = new TH1F("hProtonGenVsP", ";p (GeV/c);counts", nPBins, pMinH, pMaxH);
  auto* hProtonWithTOFVsP = new TH1F("hProtonWithTOFVsP", ";p (GeV/c);counts", nPBins, pMinH, pMaxH);

  const double mProton = TDatabasePDG::Instance()->GetParticle(2212)->Mass();

  for (int jEv = 0; jEv < nEv; ++jEv) {
    tMC->GetEvent(jEv);
    tClus->GetEvent(jEv);
    tVerTel->GetEvent(jEv);

    std::unordered_map<int, int> verTelHits;
    for (const auto& clu : verTelClus) {
      verTelHits[clu.getParticleID()]++;
    }

    std::unordered_map<int, int> tofHits;
    for (const auto& clu : tofClus) {
      tofHits[clu.getParticleID()]++;
    }

    double xVertex = 0., yVertex = 0., zVertex = 0.;
    for (const auto& p : *mcArr) {
      if (p.IsPrimary()) {
        xVertex = p.Vx();
        yVertex = p.Vy();
        zVertex = p.Vz();
        break;
      }
    }

    int iTgt = 0;
    double bestDz = std::abs(zVertex - ConstHits::targetZ[0]);
    for (int it = 1; it < ConstHits::nTargets; ++it) {
      double dz = std::abs(zVertex - ConstHits::targetZ[it]);
      if (dz < bestDz) { bestDz = dz; iTgt = it; }
    }

    int nPart = -1;
    for (const auto& p : *mcArr) {
      nPart++;
      if (!p.IsPrimary()) continue;

      const int pdg = p.GetPdgCode();
      if (pdg != 2212 && pdg != -2212) continue;

      const double rap = p.Y();
      const double eta = p.Eta();
      const double pt = p.Pt();
      const double mom = p.P();
      const double pz  = p.Pz();

      hProtonGenVsY->Fill(rap);
      hProtonGenVsEta->Fill(eta);
      hProtonGenVsPt->Fill(pt);
      hProtonGenVsP->Fill(mom);
      hProtonGenPVsY->Fill(mom, rap);
      hProtonGenPtVsEta->Fill(pt, eta);
      hGenPVsYVsPt->Fill(mom, rap, pt);

      const bool reconstructible =
        (verTelHits.count(nPart) && verTelHits.at(nPart) >= ConstHits::minVerTelHits) &&
        (tofHits.count(nPart) && tofHits.at(nPart) >= 1);

      if (reconstructible) {
        hProtonRecPVsY->Fill(mom, rap);
        hProtonRecPtVsEta->Fill(pt, eta);
      }

      if (mom < pMaxIDProton) {
        hProtonGenIDVsY->Fill(rap);
        hProtonGenIDVsEta->Fill(eta);
        hProtonGenIDVsPt->Fill(pt);

        if (reconstructible) {
          if (eta < 0){
            printf("Proton pz = % .2f GeV/c, pt = %.2f GeV/c, eta = %.2f, rapidity = %.2f, nHits = %d, nHitsTOF = %d\n", pz, pt, eta, rap, verTelHits.at(nPart), tofHits.at(nPart));
            printf("  Vertex: x=%.1f cm, y=%.1f cm, z=%.1f cm\n", xVertex, yVertex, zVertex);
            printf(" Vx = %.1f cm, Vy = %.1f cm, Vz = %.1f cm\n", p.Vx(), p.Vy(), p.Vz());
            for(const auto& clu : tofClus) {
              if (clu.getParticleID() == nPart) {
                printf("  VerTel hit: x=%.1f cm, y=%.1f cm, z=%.1f cm, time=%.2f ns\n", clu.getXLab(), clu.getYLab(), clu.getZLab(), clu.getTime() * ConstHits::nsPerSec);
              }
            }
          }
          hProtonAccVsY->Fill(rap);
          hProtonAccVsEta->Fill(eta);
          hProtonAccVsPt->Fill(pt);
          hProtonAccPVsY->Fill(mom, rap);
          hProtonAccPtVsEta->Fill(pt, eta);
        }
      }
    }

    for (const auto& clu : tofClus) {
      const int pid = clu.getParticleID();
      if (pid < 0 || pid >= (int)mcArr->size()) continue;
      const TParticle& part = (*mcArr)[pid];
      if (!part.IsPrimary()) continue;

      const int pdg = part.GetPdgCode();
      if (pdg != 2212 && pdg != -2212) continue;

      const double mom = part.P();
      const double pt  = part.Pt();
      const double rap = part.Y();
      const double eta = part.Eta();

      const double dx = clu.getXLab() - part.Vx();
      const double dy = clu.getYLab() - part.Vy();
      const double dz = clu.getZLab() - part.Vz();
      const double pathLen = std::sqrt(dx * dx + dy * dy + dz * dz);
      const double tof_ns = clu.getTime() * ConstHits::nsPerSec;
      if (tof_ns <= 0.) continue;

      const double beta = pathLen / (ConstHits::cLight * tof_ns);
      hBetaVsP->Fill(mom, beta);
      hBetaVsPtgt[iTgt]->Fill(mom, beta);
      hBetaVsPprotons->Fill(mom, beta);

      hProtonWithTOFVsY->Fill(rap);
      hProtonWithTOFVsEta->Fill(eta);
      hProtonWithTOFVsPt->Fill(pt);
      hProtonWithTOFVsP->Fill(mom);

      if (mom < pMaxIDProton && verTelHits.count(pid) && verTelHits.at(pid) >= ConstHits::minVerTelHits) {
        hProtonWithTOFIDVsY->Fill(rap);
        hProtonWithTOFIDVsEta->Fill(eta);
        hProtonWithTOFIDVsPt->Fill(pt);
      }
    }
  }

  auto* fOut = TFile::Open(Form("%s/TOFplots_hits.root", dirSimu), "recreate");

  for (auto* h : {hBetaVsP, hBetaVsPprotons}) h->Write();
  for (auto* h : {hProtonGenVsY, hProtonAccVsY, hProtonWithTOFVsY, hProtonGenIDVsY, hProtonWithTOFIDVsY}) h->Write();
  for (auto* h : {hProtonGenVsEta, hProtonAccVsEta, hProtonWithTOFVsEta, hProtonGenIDVsEta, hProtonWithTOFIDVsEta}) h->Write();
  for (auto* h : {hProtonGenVsPt, hProtonAccVsPt, hProtonWithTOFVsPt, hProtonGenIDVsPt, hProtonWithTOFIDVsPt}) h->Write();
  for (auto* h : {hProtonGenVsP, hProtonWithTOFVsP}) h->Write();

  hProtonGenPVsY->Write();
  hProtonRecPVsY->Write();
  hProtonAccPVsY->Write();
  hProtonGenPtVsEta->Write();
  hProtonRecPtVsEta->Write();
  hProtonAccPtVsEta->Write();
  hGenPVsYVsPt->Write();
  for (int it = 0; it < ConstHits::nTargets; ++it) hBetaVsPtgt[it]->Write();

  auto cloneH1 = [](TH1F* h, const char* newName) -> TH1F* {
    auto* c = (TH1F*)h->Clone(newName); c->SetDirectory(nullptr); return c;
  };

  auto* hTOFIDrap_clamped = cloneH1(hProtonWithTOFIDVsY, "hProtonWithTOFID_clamped");
  auto* hTOFIDeta_clamped = cloneH1(hProtonWithTOFIDVsEta, "hProtonWithTOFIDEta_clamped");
  auto* hTOFIDpt_clamped = cloneH1(hProtonWithTOFIDVsPt, "hProtonWithTOFIDPt_clamped");

  auto* effRap = makeEfficiencyHits(hProtonWithTOFVsY, hProtonGenVsY,
    "effProtonVsRapidity", "Proton TOF acceptance vs rapidity;y;#varepsilon");
  auto* effEta = makeEfficiencyHits(hProtonWithTOFVsEta, hProtonGenVsEta,
    "effProtonVsEta", "Proton TOF acceptance vs #eta;#eta;#varepsilon");
  auto* effPt = makeEfficiencyHits(hProtonWithTOFVsPt, hProtonGenVsPt,
    "effProtonVsPt", "Proton TOF acceptance vs p_{T};p_{T} (GeV/c);#varepsilon");

  auto* acceptanceRap = makeEfficiencyHits(hProtonAccVsY, hProtonGenVsY,
    "accProtonVsRapidity",
    Form("Proton acceptance (4 VerTel hits, p < %.0f GeV/c) vs y;y;Acceptance", (double)pMaxIDProton));
  auto* acceptanceEta = makeEfficiencyHits(hProtonAccVsEta, hProtonGenVsEta,
    "accProtonVsEta",
    Form("Proton acceptance (4 VerTel hits, p < %.0f GeV/c) vs #eta;#eta;Acceptance", (double)pMaxIDProton));
  auto* acceptancePt = makeEfficiencyHits(hProtonAccVsPt, hProtonGenVsPt,
    "accProtonVsPt",
    Form("Proton acceptance (4 VerTel hits, p < %.0f GeV/c) vs p_{T};p_{T} (GeV/c);Acceptance", (double)pMaxIDProton));

  auto* effTrkRap = makeEfficiencyHits(hTOFIDrap_clamped, hProtonAccVsY,
    "effTrackingVsRapidity",
    Form("Proton tracking efficiency (p < %.0f GeV/c) vs y;y;Tracking efficiency", (double)pMaxIDProton));
  auto* effTrkEta = makeEfficiencyHits(hTOFIDeta_clamped, hProtonAccVsEta,
    "effTrackingVsEta",
    Form("Proton tracking efficiency (p < %.0f GeV/c) vs #eta;#eta;Tracking efficiency", (double)pMaxIDProton));
  auto* effTrkPt = makeEfficiencyHits(hTOFIDpt_clamped, hProtonAccVsPt,
    "effTrackingVsPt",
    Form("Proton tracking efficiency (p < %.0f GeV/c) vs p_{T};p_{T} (GeV/c);Tracking efficiency", (double)pMaxIDProton));

  // 2D efficiencies: pT vs eta
  auto* hProtonWithTOFID_clamped_2D = (TH2F*)hProtonRecPtVsEta->Clone("hProtonWithTOFID_clamped_2D");
  hProtonWithTOFID_clamped_2D->SetDirectory(nullptr);
  
  auto* eff2D_TOFAcceptance = makeEfficiencyHits(hProtonWithTOFID_clamped_2D, hProtonGenPtVsEta,
    "eff2D_TOFAcceptance_PtVsEta",
    "Proton TOF acceptance vs (p_{T}, #eta);p_{T} (GeV/c);#eta");
  
  auto* hProtonAcc_clamped_2D = (TH2F*)hProtonAccPtVsEta->Clone("hProtonAcc_clamped_2D");
  hProtonAcc_clamped_2D->SetDirectory(nullptr);
  
  auto* eff2D_GeometricAcc = makeEfficiencyHits(hProtonAcc_clamped_2D, hProtonGenPtVsEta,
    "eff2D_GeometricAcceptance_PtVsEta",
    Form("Proton geometric acceptance (4 VerTel hits, p < %.0f GeV/c) vs (p_{T}, #eta);p_{T} (GeV/c);#eta",
         (double)pMaxIDProton));

  for (auto* eff : {effRap, effEta, effPt, acceptanceRap, acceptanceEta, acceptancePt, effTrkRap, effTrkEta, effTrkPt,
                    eff2D_TOFAcceptance, eff2D_GeometricAcc}) {
    if (eff) eff->Write();
  }

  auto* c1 = new TCanvas("cBetaVsP_hits", "#beta vs p (hits-only)", 1800, 1400);
  c1->SetLeftMargin(0.12); c1->SetBottomMargin(0.12);
  c1->SetRightMargin(0.14); c1->SetTopMargin(0.08);
  c1->SetLogz();
  hBetaVsP->GetYaxis()->SetTitleOffset(0.9);
  hBetaVsP->GetYaxis()->SetRangeUser(0., 1.1);
  hBetaVsP->Draw("colz");
  drawBetaCurveHits(ConstHits::mPion,  kRed);
  drawBetaCurveHits(ConstHits::mKaon,  kBlue);
  drawBetaCurveHits(mProton,           kGreen + 2);
  c1->SaveAs(Form("%s/tof_hits_beta_vs_p.png", dirSimu));

  auto drawEff = [&](TEfficiency* eff, const char* cName, const char* cTitle, const char* outFile) {
    if (!eff) return;
    auto* c = new TCanvas(cName, cTitle, 800, 600);
    c->SetLeftMargin(0.14); c->SetBottomMargin(0.13);
    c->SetRightMargin(0.05); c->SetTopMargin(0.07);
    eff->SetMarkerStyle(20);
    eff->SetMarkerSize(0.8);
    eff->SetLineWidth(2);
    eff->Draw("AP");
    gPad->Modified(); gPad->Update();
    if (auto* gr = eff->GetPaintedGraph()) {
      gr->GetYaxis()->SetRangeUser(0., 1.05);
      gr->GetXaxis()->SetTitleSize(0.05);
      gr->GetYaxis()->SetTitleSize(0.05);
    }
    gPad->Modified(); gPad->Update();
    c->SaveAs(Form("%s/%s", dirSimu, outFile));
  };

  drawEff(effRap, "cEffRapHits", "Proton acceptance vs y (hits-only)", "tof_hits_proton_acc_rapidity.png");
  drawEff(effEta, "cEffEtaHits", "Proton acceptance vs #eta (hits-only)", "tof_hits_proton_acc_eta.png");
  drawEff(effPt, "cEffPtHits", "Proton acceptance vs p_{T} (hits-only)", "tof_hits_proton_acc_pt.png");
  drawEff(acceptanceRap, "cAccRapHits", "Proton geometric acceptance vs y (hits-only)", "tof_hits_proton_geoacc_rapidity.png");
  drawEff(acceptanceEta, "cAccEtaHits", "Proton geometric acceptance vs #eta (hits-only)", "tof_hits_proton_geoacc_eta.png");
  drawEff(acceptancePt, "cAccPtHits", "Proton geometric acceptance vs p_{T} (hits-only)", "tof_hits_proton_geoacc_pt.png");
  drawEff(effTrkRap, "cEffTrkRapHits", "Proton tracking efficiency vs y (hits-only)", "tof_hits_proton_efftrk_rapidity.png");
  drawEff(effTrkEta, "cEffTrkEtaHits", "Proton tracking efficiency vs #eta (hits-only)", "tof_hits_proton_efftrk_eta.png");
  drawEff(effTrkPt, "cEffTrkPtHits", "Proton tracking efficiency vs p_{T} (hits-only)", "tof_hits_proton_efftrk_pt.png");

  auto drawH2 = [&](TH2F* h, const char* cName, const char* cTitle, const char* outFile) {
    auto* c = new TCanvas(cName, cTitle, 900, 700);
    c->SetLeftMargin(0.14); c->SetBottomMargin(0.12);
    c->SetRightMargin(0.14); c->SetTopMargin(0.07);
    c->SetLogz();
    h->Draw("colz");
    c->SaveAs(Form("%s/%s", dirSimu, outFile));
  };

  auto drawH1 = [&](TH1F* h, const char* cName, const char* cTitle, const char* outFile) {
    auto* c = new TCanvas(cName, cTitle, 900, 700);
    c->SetLeftMargin(0.14); c->SetBottomMargin(0.12);
    c->SetRightMargin(0.05); c->SetTopMargin(0.07);
    h->SetLineWidth(2);
    h->SetMarkerStyle(20);
    h->SetMarkerSize(0.8);
    h->Draw("E P");
    c->SaveAs(Form("%s/%s", dirSimu, outFile));
  };

  drawH1(hProtonGenVsEta, "cProtonGenVsEtaHits", "Generated protons vs #eta (hits-only)", "tof_hits_proton_generated_eta.png");
  drawH1(hProtonWithTOFVsEta, "cProtonRecVsEtaHits", "Reconstructed protons vs #eta (hits-only)", "tof_hits_proton_reconstructed_eta.png");
  drawH1(hProtonAccVsEta, "cProtonAccVsEtaHits", "Accepted protons vs #eta (hits-only)", "tof_hits_proton_accepted_eta.png");

  drawH2(hProtonGenPtVsEta, "cProtonGenPtVsEtaHits", "Generated protons: p_{T} vs #eta (hits-only)", "tof_hits_proton_generated_pt_vs_eta.png");
  drawH2(hProtonRecPtVsEta, "cProtonRecPtVsEtaHits", "Reconstructed protons: p_{T} vs #eta (hits-only)", "tof_hits_proton_reconstructed_pt_vs_eta.png");
  drawH2(hProtonAccPtVsEta, "cProtonAccPtVsEtaHits", "Accepted protons: p_{T} vs #eta (hits-only)", "tof_hits_proton_accepted_pt_vs_eta.png");

  // Helper: draw a 2D TEfficiency canvas and save
  auto drawEff2D = [&](TEfficiency* eff, const char* cName, const char* cTitle, const char* outFile) {
    if (!eff) return;
    auto* c = new TCanvas(cName, cTitle, 1000, 800);
    c->SetLeftMargin(0.14); c->SetBottomMargin(0.12);
    c->SetRightMargin(0.15); c->SetTopMargin(0.07);
    c->SetLogz(0);
    eff->Draw("colz");
    c->SaveAs(Form("%s/%s", dirSimu, outFile));
  };

  drawEff2D(eff2D_TOFAcceptance, "cEff2D_TOFAcceptanceHits", "Proton TOF acceptance vs (p_{T}, #eta) (hits-only)",
            "tof_hits_proton_eff2d_tof_acceptance_pt_vs_eta.png");
  drawEff2D(eff2D_GeometricAcc, "cEff2D_GeometricAccHits", "Proton geometric acceptance vs (p_{T}, #eta) (hits-only)",
            "tof_hits_proton_eff2d_geometric_acceptance_pt_vs_eta.png");

  fOut->Close();
  tMC->GetCurrentFile()->Close();
  tClus->GetCurrentFile()->Close();
  tVerTel->GetCurrentFile()->Close();

  printf("Done (hits-only) — plots saved in %s/\n", dirSimu);
}
