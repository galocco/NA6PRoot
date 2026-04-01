#if !defined(__CINT__) || defined(__MAKECINT__)
#include <TTree.h>
#include <TFile.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TEfficiency.h>
#include <TParticle.h>
#include <TCanvas.h>
#include <TLine.h>
#include <TLegend.h>
#include <TLatex.h>
#include <TStyle.h>
#include <TDatabasePDG.h>
#include <TMath.h>
#include <cmath>
#include "NA6PTrack.h"
#include "NA6PLayoutParam.h"
#include "NA6PTOFCluster.h"
#include "NA6PVerTelCluster.h"
#endif

/// Compute 3D path length [cm] of a charged track in a uniform By field.
/// Coordinates and momenta are in the lab frame (cm, GeV/c).
/// By is the magnetic field y-component in Tesla.
double getPathLength(
  double x0, double y0, double z0, // start point [cm]
  double x1, double y1, double z1, // end   point [cm]
  double px, double py, double pz, // momentum at start [GeV/c]
  int q,                           // charge sign (+1 or -1)
  double By = -1.47                // Tesla
)
{
  const double k = 0.299792458; // GeV/(c·T·m)

  // pT in the bending plane (x-z, perpendicular to By)
  double pT = std::sqrt(px * px + pz * pz);
  double p = std::sqrt(px * px + py * py + pz * pz);

  if (pT < 1e-12 || std::abs(By) < 1e-12) {
    // straight-line fallback
    double dx = x1 - x0;
    double dy = y1 - y0;
    double dz = z1 - z0;
    return std::sqrt(dx * dx + dy * dy + dz * dz); // cm
  }

  // bending radius in cm  (R[m] = pT / (k·|q|·|B|) → R[cm] = 100·R[m])
  double R = 100.0 * pT / (k * std::abs(q) * std::abs(By)); // cm

  // unit tangent in the x-z plane
  double tx = px / pT;
  double tz = pz / pT;

  // inward normal (centre direction)
  double nx = -tz;
  double nz = tx;
  double sgn = (q * By > 0) ? +1.0 : -1.0;

  double xc = x0 + sgn * R * nx; // cm
  double zc = z0 + sgn * R * nz;

  double phi0 = std::atan2(z0 - zc, x0 - xc);
  double phi1 = std::atan2(z1 - zc, x1 - xc);
  double dphi = phi1 - phi0;
  while (dphi > M_PI)
    dphi -= 2.0 * M_PI;
  while (dphi < -M_PI)
    dphi += 2.0 * M_PI;

  // 3D arc length = R·|Δφ| · (p/pT)   [cm]
  double L = std::abs(R * dphi) * (p / pT);
  return L; // cm
}

/// Save projections of h2 onto the beta axis in momentum slices of width deltaP.
/// Creates a subdirectory inside the current TFile directory, named after the histogram.
void saveBetaSlices(TH2F* h2, double deltaP)
{
  TDirectory* curDir = gDirectory;
  TDirectory* sliceDir = curDir->mkdir(Form("BetaSlices_%s", h2->GetName()));
  sliceDir->cd();

  const TAxis* ax = h2->GetXaxis();
  double pMin = ax->GetXmin();
  double pMax = ax->GetXmax();

  for (double pLo = pMin; pLo + 1.e-9 < pMax; pLo += deltaP) {
    double pHi = std::min(pLo + deltaP, pMax);
    int binLo = ax->FindBin(pLo + 1.e-9);
    int binHi = ax->FindBin(pHi - 1.e-9);
    if (binLo < 1)
      binLo = 1;
    if (binHi > ax->GetNbins())
      binHi = ax->GetNbins();

    TH1D* hSlice = h2->ProjectionY(Form("%s_p%.0f_%.0f", h2->GetName(), pLo * 1000, pHi * 1000),
                                   binLo, binHi);
    hSlice->SetTitle(Form("#beta   (%.2f < #it{p} < %.2f GeV/#it{c});#beta;counts",
                          pLo, pHi));
    hSlice->Write();
  }

  curDir->cd();
}

void plotTOF(const char* dirSimu = ".", float pMaxIDProton = 3)
{
  gStyle->SetOptStat(0);
  gStyle->SetTitleSize(0.05, "XY");
  gStyle->SetLabelSize(0.045, "XY");
  gStyle->SetTitleOffset(1.1, "X");
  gStyle->SetTitleOffset(1.0, "Y");

  const double cLight = 29.9792458; // cm/ns

  // ── open input files ──────────────────────────────────────────────
  auto* fTracks = TFile::Open(Form("%s/TracksTOFMatching.root", dirSimu));
  if (!fTracks || fTracks->IsZombie()) {
    printf("Cannot open TracksTOFMatching.root\n");
    return;
  }
  auto* tTracks = (TTree*)fTracks->Get("tracksTOFMatching");
  if (!tTracks) {
    printf("Cannot find tree 'tracksTOFMatching'\n");
    return;
  }
  std::vector<NA6PTrack> tracks, *tracksPtr = &tracks;
  tTracks->SetBranchAddress("TOFMatching", &tracksPtr);

  auto* fMC = TFile::Open(Form("%s/MCKine.root", dirSimu));
  if (!fMC || fMC->IsZombie()) {
    printf("Cannot open MCKine.root\n");
    return;
  }
  auto* tMC = (TTree*)fMC->Get("mckine");
  if (!tMC) {
    printf("Cannot find tree 'mckine'\n");
    return;
  }
  std::vector<TParticle>* mcArr = nullptr;
  tMC->SetBranchAddress("tracks", &mcArr);

  auto* fClus = TFile::Open(Form("%s/ClustersTOF.root", dirSimu));
  if (!fClus || fClus->IsZombie()) {
    printf("Cannot open ClustersTOF.root\n");
    return;
  }
  auto* tClus = (TTree*)fClus->Get("clustersTOF");
  if (!tClus) {
    printf("Cannot find tree 'clustersTOF'\n");
    return;
  }
  std::vector<NA6PTOFCluster> tofClus, *tofClusPtr = &tofClus;
  tClus->SetBranchAddress("TOF", &tofClusPtr);

  // VerTel clusters
  auto* fVerTel = TFile::Open(Form("%s/ClustersVerTel.root", dirSimu));
  if (!fVerTel || fVerTel->IsZombie()) {
    printf("Cannot open ClustersVerTel.root\n");
    return;
  }
  auto* tVerTel = (TTree*)fVerTel->Get("clustersVerTel");
  if (!tVerTel) {
    printf("Cannot find tree 'clustersVerTel'\n");
    return;
  }
  std::vector<NA6PVerTelCluster> verTelClus, *verTelClusPtr = &verTelClus;
  tVerTel->SetBranchAddress("VerTel", &verTelClusPtr);

  int nEv = std::min({(int)tTracks->GetEntries(), (int)tMC->GetEntries(), (int)tClus->GetEntries()});
  printf("Processing %d events\n", nEv);

  // ── histograms ────────────────────────────────────────────────────
  // beta vs p  (all tracks with TOF)
  auto* hBetaVsP = new TH2F("hBetaVsP", ";#it{p} (GeV/#it{c});#beta",
                            50, 0., 5., 300, 0.4, 1.10);

  // beta vs p  (protons only)
  auto* hBetaVsPprotons = new TH2F("hBetaVsPprotons", ";#it{p} (GeV/#it{c});#beta",
                                   50, 0., 5., 300, 0.4, 1.10);

  // per-target beta vs p
  const int nTargets = 5;
  const double targetZ[nTargets] = {-5.3, -4.1, -2.9, -1.7, -0.5};
  TH2F* hBetaVsPtgt[nTargets];
  for (int it = 0; it < nTargets; ++it)
    hBetaVsPtgt[it] = new TH2F(Form("hBetaVsP_tgt%d", it),
                               Form("target %d (z=%.1f cm);#it{p} (GeV/#it{c});#beta", it, targetZ[it]),
                               50, 0., 5., 300, 0.4, 1.10);

  // mass² vs p  (m² = p²(1/β² − 1))
  auto* hMass2VsP = new TH2F("hMass2VsP", "m^{2} vs #it{p};#it{p} (GeV/#it{c});m^{2} (GeV^{2}/#it{c}^{4})",
                             200, 0., 20., 400, -1., 4.);

  // proton acceptance: generated vs matched in rapidity
  const int nRapBins = 50;
  const double rapMin = -2., rapMax = 6.;
  auto* hProtonGen = new TH1F("hProtonGen", ";y;counts", nRapBins, rapMin, rapMax);
  auto* hProtonTOF = new TH1F("hProtonTOF", ";y;counts", nRapBins, rapMin, rapMax);
  auto* hProtonAcc = new TH1F("hProtonAcc", ";y;counts", nRapBins, rapMin, rapMax);

  auto* hProtonGenID = new TH1F("hProtonGenID", ";y;counts", nRapBins, rapMin, rapMax);
  auto* hProtonTOFID = new TH1F("hProtonTOFID", ";y;counts", nRapBins, rapMin, rapMax);

  // proton acceptance in pT
  const int nPtBins = 50;
  const double ptMin = 0., ptMax = 5.;
  auto* hProtonGenPt = new TH1F("hProtonGenPt", ";p_{T} (GeV/c);counts", nPtBins, ptMin, ptMax);
  auto* hProtonTOFPt = new TH1F("hProtonTOFPt", ";p_{T} (GeV/c);counts", nPtBins, ptMin, ptMax);
  auto* hProtonAccPt = new TH1F("hProtonAccPt", ";p_{T} (GeV/c);counts", nPtBins, ptMin, ptMax);
  auto* hProtonGenIDPt = new TH1F("hProtonGenIDPt", ";p_{T} (GeV/c);counts", nPtBins, ptMin, ptMax);
  auto* hProtonTOFIDPt = new TH1F("hProtonTOFIDPt", ";p_{T} (GeV/c);counts", nPtBins, ptMin, ptMax);
  // proton acceptance in p
  const int nPBins = 50;
  const double pMin = 0., pMax = 50.;
  auto* hProtonGenP = new TH1F("hProtonGenP", ";p (GeV/c);counts", nPBins, pMin, pMax);
  auto* hProtonTOFP = new TH1F("hProtonTOFP", ";p (GeV/c);counts", nPBins, pMin, pMax);

  const double mProton = TDatabasePDG::Instance()->GetParticle(2212)->Mass();

  // ── event loop ────────────────────────────────────────────────────
  for (int jEv = 0; jEv < nEv; ++jEv) {
    tMC->GetEvent(jEv);
    tTracks->GetEvent(jEv);
    tClus->GetEvent(jEv);
    tVerTel->GetEvent(jEv);

    // create a map of VerTel clusters by their particle ID and number of hits (to check for reconstructibility)
    std::unordered_map<int, int> verTelMap;
    for (const auto& clu : verTelClus) {
      int pid = clu.getParticleID();
      if (verTelMap.find(pid) == verTelMap.end())
        verTelMap[pid] = 1;
      else
        verTelMap[pid]++;
    }
    std::unordered_map<int, int> tofMap;
    for (const auto& clu : tofClus) {
      int pid = clu.getParticleID();
      if (tofMap.find(pid) == tofMap.end())
        tofMap[pid] = 1;
      else
        tofMap[pid]++;
    }

    // find primary vertex Z from MC
    double zVertex = 0.;
    double xVertex = 0.;
    double yVertex = 0.;
    for (const auto& p : *mcArr) {
      if (p.IsPrimary()) {
        zVertex = p.Vz();
        xVertex = p.Vx();
        yVertex = p.Vy();
        break;
      }
    }

    // identify closest target
    int iTgt = 0;
    double bestDz = std::abs(zVertex - targetZ[0]);
    for (int it = 1; it < nTargets; ++it) {
      double dz = std::abs(zVertex - targetZ[it]);
      if (dz < bestDz) {
        bestDz = dz;
        iTgt = it;
      }
    }

    // fill generated proton denominators
    int nPart = -1;
    for (const auto& p : *mcArr) {
      nPart++;
      if (!p.IsPrimary())
        continue;
      int pdg = p.GetPdgCode();
      if (pdg == 2212 || pdg == -2212) {
        double rap = p.Y();
        double pt = p.Pt();
        hProtonGen->Fill(rap);
        hProtonGenPt->Fill(pt);
        hProtonGenP->Fill(p.P());
        if (p.P() < pMaxIDProton) {
          hProtonGenID->Fill(rap);
          hProtonGenIDPt->Fill(pt);

          
          if (verTelMap.find(nPart) != verTelMap.end() && verTelMap[nPart] >= 4 &&
              tofMap.find(nPart) != tofMap.end() && tofMap[nPart] >= 1) {
            hProtonAcc->Fill(rap);
            hProtonAccPt->Fill(pt);
          }
        }
      }
    }

    for (const auto& trk : tracks) {
      double tof = trk.getTOF(); // ns
      if (tof <= 0.)
        continue; // no TOF match

      double p = trk.getP();
      if (p < 1.e-6)
        continue;

      // Find the TOF cluster matched to this track (by matching time)
      double xTOF = 0., yTOF = 0., zTOF = 0.;
      bool foundCluster = false;
      for (const auto& clu : tofClus) {
        if (std::abs(clu.getTime() - tof) < 1.e-15) { // exact match on stored time
          xTOF = clu.getXLab();
          yTOF = clu.getYLab();
          zTOF = clu.getZLab();
          foundCluster = true;
          break;
        }
      }
      if (!foundCluster)
        continue; // skip if no cluster found

      // path length: straight-line from vertex to TOF cluster position [cm]
      double dx = xTOF - xVertex;
      double dy = yTOF - yVertex;
      double dz = zTOF - zVertex;
      double pathLen = std::sqrt(dx * dx + dy * dy + dz * dz); // cm

      // TOF from Geant4 is in seconds — convert to ns
      double tof_ns = tof * 1.e9;

      double beta = pathLen / (cLight * tof_ns); // cLight = 29.9792 cm/ns

      hBetaVsP->Fill(p, beta);
      hBetaVsPtgt[iTgt]->Fill(p, beta);

      if (beta > 0. && beta < 10.) { // sanity cut
        double mass2 = p * p * (1. / (beta * beta) - 1.);
        hMass2VsP->Fill(p, mass2);
      }

      // MC-truth proton plots
      int pid = trk.getParticleID();
      if (pid >= 0 && pid < (int)mcArr->size()) {
        int pdg = (*mcArr)[pid].GetPdgCode();
        if (pdg == 2212 || pdg == -2212) {
          double rap = (*mcArr)[pid].Y();
          double pt = (*mcArr)[pid].Pt();
          hProtonTOF->Fill(rap);

          hProtonTOFPt->Fill(pt);
          hProtonTOFP->Fill(p);
          hBetaVsPprotons->Fill(p, beta);

          if (p < pMaxIDProton) {
            hProtonTOFID->Fill(rap);
            hProtonTOFIDPt->Fill(pt);
          }
        }
      }
    }
  }

  // ── output file ───────────────────────────────────────────────────
  auto* fOut = TFile::Open(Form("%s/TOFplots.root", dirSimu), "recreate");
  hBetaVsP->Write();
  hBetaVsPprotons->Write();
  hMass2VsP->Write();
  hProtonGen->Write();
  hProtonAcc->Write();
  hProtonAccPt->Write();
  hProtonTOF->Write();
  hProtonGenPt->Write();
  hProtonTOFPt->Write();
  hProtonGenP->Write();
  hProtonTOFP->Write();
  hProtonGenID->Write();
  hProtonTOFID->Write();

  // ── per-target histograms ──────────────────────────────────────────
  for (int it = 0; it < nTargets; ++it)
    hBetaVsPtgt[it]->Write();

  // ── beta slices in p intervals ────────────────────────────────────
  saveBetaSlices(hBetaVsP, 0.2);
  saveBetaSlices(hBetaVsPprotons, 0.2);
  for (int it = 0; it < nTargets; ++it)
    saveBetaSlices(hBetaVsPtgt[it], 0.2);

  // ── TEfficiency for proton acceptance ─────────────────────────────
  TEfficiency* effRap = nullptr;
  if (TEfficiency::CheckConsistency(*hProtonTOF, *hProtonGen)) {
    effRap = new TEfficiency(*hProtonTOF, *hProtonGen);
    effRap->SetName("effProtonVsRapidity");
    effRap->SetTitle("Proton TOF acceptance vs rapidity;y;#varepsilon");
    effRap->Write();
  }
  TEfficiency* effPt = nullptr;
  if (TEfficiency::CheckConsistency(*hProtonTOFPt, *hProtonGenPt)) {
    effPt = new TEfficiency(*hProtonTOFPt, *hProtonGenPt);
    effPt->SetName("effProtonVsPt");
    effPt->SetTitle("Proton TOF acceptance vs p_{T};p_{T} (GeV/c);#varepsilon");
    effPt->Write();
  }
  TEfficiency* acceptanceRap = nullptr;
  if (TEfficiency::CheckConsistency(*hProtonAcc, *hProtonGen)) {
    printf("Creating TEfficiency for proton acceptance vs rapidity...\n");
    acceptanceRap = new TEfficiency(*hProtonAcc, *hProtonGen);
    acceptanceRap->SetName("accProtonVsRapidity");
    acceptanceRap->SetTitle(Form("Proton acceptance (with 4 VerTel hits and p < %g) vs rapidity;y;Acceptance", pMaxIDProton));
    acceptanceRap->Write();
  }
  TEfficiency* acceptancePt = nullptr;
  if (TEfficiency::CheckConsistency(*hProtonAccPt, *hProtonGenPt)) {
    printf("Creating TEfficiency for proton acceptance vs pT...\n");
    acceptancePt = new TEfficiency(*hProtonAccPt, *hProtonGenPt);
    acceptancePt->SetName("accProtonVsPt");
    acceptancePt->SetTitle(Form("Proton acceptance (with 4 VerTel hits and p < %g) vs p_{T};p_{T} (GeV/c);Acceptance", pMaxIDProton));
    acceptancePt->Write();
  }

  TEfficiency* effTrkRap = nullptr;
  for (int bin = 1; bin <= hProtonTOFID->GetNbinsX(); ++bin) {
    double genCount = hProtonAcc->GetBinContent(bin);
    double accCount = hProtonTOFID->GetBinContent(bin);
    if (accCount > genCount) {
      printf("Warning: in bin %d, TOF ID count (%.0f) exceeds generated count (%.0f). Setting TOF ID count to generated count for TEfficiency.\n",
             bin, accCount, genCount);
      hProtonTOFID->SetBinContent(bin, genCount);
    }
  }

  if (TEfficiency::CheckConsistency(*hProtonTOFID, *hProtonAcc)) {
    printf("Creating TEfficiency for proton tracking efficiency vs rapidity...\n");
    effTrkRap = new TEfficiency(*hProtonTOFID, *hProtonAcc);
    effTrkRap->SetName("effTrackingVsRapidity");
    effTrkRap->SetTitle(Form("Proton TOF effTrk (with 4 VerTel hits and p < %g) vs rapidity;y;Tracking efficiency", pMaxIDProton));
    effTrkRap->Write();
  }
  TEfficiency* effTrkPt = nullptr;
  if (TEfficiency::CheckConsistency(*hProtonTOFIDPt, *hProtonGenPt)) {
    printf("Creating TEfficiency for proton tracking efficiency vs pT...\n");
    effTrkPt = new TEfficiency(*hProtonTOFIDPt, *hProtonGenPt);
    effTrkPt->SetName("effTrackingVsPt");
    effTrkPt->SetTitle(Form("Proton TOF effTrk (with 4 VerTel hits and p < %g) vs p_{T};p_{T} (GeV/c);Tracking efficiency", pMaxIDProton));
    effTrkPt->Write();
  }
  // ── draw ──────────────────────────────────────────────────────────

  // 1) beta vs p
  auto* c1 = new TCanvas("cBetaVsP", "#beta vs p", 1800, 1400);
  c1->SetLeftMargin(0.12);
  c1->SetBottomMargin(0.12);
  c1->SetRightMargin(0.14);
  c1->SetTopMargin(0.08);
  c1->SetLogz();
  hBetaVsP->GetYaxis()->SetTitleOffset(0.9);
  hBetaVsP->GetYaxis()->SetRangeUser(0., 1.1);
  hBetaVsP->Draw("colz");
  // expected curves
  auto drawExpected = [&](double mass, int color, const char* label) {
    const int np = 200;
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
  };
  drawExpected(0.13957, kRed, "#pi");     // pion
  drawExpected(0.49368, kBlue, "K");      // kaon
  drawExpected(mProton, kGreen + 2, "p"); // proton
  c1->SaveAs(Form("%s/tof_beta_vs_p.png", dirSimu));

  // 3) mass² vs p
  auto* c2 = new TCanvas("cMass2VsP", "m^{2} vs p", 1800, 1400);
  c2->SetLeftMargin(0.12);
  c2->SetBottomMargin(0.12);
  c2->SetRightMargin(0.14);
  c2->SetTopMargin(0.08);
  c2->SetLogz();
  hMass2VsP->GetYaxis()->SetTitleOffset(0.9);
  hMass2VsP->Draw("colz");
  auto* linePi = new TLine(0., 0.13957 * 0.13957, 20., 0.13957 * 0.13957);
  linePi->SetLineColor(kRed);
  linePi->SetLineStyle(2);
  linePi->Draw();
  auto* lineK = new TLine(0., 0.49368 * 0.49368, 20., 0.49368 * 0.49368);
  lineK->SetLineColor(kBlue);
  lineK->SetLineStyle(2);
  lineK->Draw();
  auto* lineP = new TLine(0., mProton * mProton, 20., mProton * mProton);
  lineP->SetLineColor(kGreen + 2);
  lineP->SetLineStyle(2);
  lineP->Draw();
  c2->SaveAs(Form("%s/tof_mass2_vs_p.png", dirSimu));

  // 5) proton acceptance vs rapidity
  if (effRap) {
    auto* c4 = new TCanvas("cEffRap", "Proton acceptance vs y", 800, 600);
    c4->SetLeftMargin(0.14);
    c4->SetBottomMargin(0.13);
    c4->SetRightMargin(0.05);
    c4->SetTopMargin(0.07);
    effRap->SetMarkerStyle(20);
    effRap->SetMarkerSize(0.8);
    effRap->SetLineWidth(2);
    effRap->Draw("AP");
    gPad->Modified();
    gPad->Update();
    if (auto* gr = effRap->GetPaintedGraph()) {
      gr->GetYaxis()->SetRangeUser(0., 1.05);
      gr->GetXaxis()->SetTitleSize(0.05);
      gr->GetYaxis()->SetTitleSize(0.05);
    }
    gPad->Modified();
    gPad->Update();
    c4->SaveAs(Form("%s/tof_proton_acc_rapidity.png", dirSimu));
  }

  // 6) proton acceptance vs pT
  if (effPt) {
    auto* c5 = new TCanvas("cEffPt", "Proton acceptance vs p_{T}", 800, 600);
    c5->SetLeftMargin(0.14);
    c5->SetBottomMargin(0.13);
    c5->SetRightMargin(0.05);
    c5->SetTopMargin(0.07);
    effPt->SetMarkerStyle(20);
    effPt->SetMarkerSize(0.8);
    effPt->SetLineWidth(2);
    effPt->Draw("AP");
    gPad->Modified();
    gPad->Update();
    if (auto* gr = effPt->GetPaintedGraph()) {
      gr->GetYaxis()->SetRangeUser(0., 1.05);
      gr->GetXaxis()->SetTitleSize(0.05);
      gr->GetYaxis()->SetTitleSize(0.05);
    }
    gPad->Modified();
    gPad->Update();
    c5->SaveAs(Form("%s/tof_proton_acc_pt.png", dirSimu));
  }

  fOut->Close();
  fTracks->Close();
  fMC->Close();
  fClus->Close();

  printf("Done — plots saved in %s/\n", dirSimu);
}
