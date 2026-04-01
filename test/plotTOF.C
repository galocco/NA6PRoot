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

// ============================================================
//  Physical / detector constants (all in cm, ns, GeV/c, T)
// ============================================================
namespace Const {
  constexpr double cLight      = 29.9792458; // cm/ns
  constexpr double kMagFactor  = 0.299792458; // GeV/(c·T·m) — for radius formula
  constexpr double By          = -1.47;       // magnetic field [T], y-component
  constexpr double mPion       = 0.13957;     // GeV/c²
  constexpr double mKaon       = 0.49368;     // GeV/c²
  constexpr double nsPerSec    = 1.e9;        // unit conversion: s → ns
  constexpr double tofTimeEps  = 1.e-15;      // [ns] tolerance for cluster time matching
  constexpr double pTmin       = 1.e-12;      // [GeV/c] guard against zero-division
  constexpr double Bmin        = 1.e-12;      // [T]     guard against zero-division

  // Target positions along z [cm]
  constexpr int    nTargets    = 5;
  constexpr double targetZ[nTargets] = {-5.3, -4.1, -2.9, -1.7, -0.5};

  // Minimum VerTel hits for a track to be considered reconstructible
  constexpr int minVerTelHits = 4;
}

// ============================================================
/// Compute the 3-D arc length [cm] of a charged track in a
/// uniform By field (bending in the x-z plane).
///
/// All positions in cm, momenta in GeV/c, field in Tesla.
/// Falls back to a straight-line distance when pT ≈ 0 or B ≈ 0.
// ============================================================
double getPathLength(
  double x0, double y0, double z0,  // start point [cm]
  double x1, double y1, double z1,  // end   point [cm]
  double px, double py, double pz,  // momentum at start [GeV/c]
  int    q,                         // charge sign: +1 or -1
  double By = Const::By             // field [T]
)
{
  const double pT = std::sqrt(px * px + pz * pz); // transverse momentum [GeV/c]
  const double p  = std::sqrt(px * px + py * py + pz * pz);

  // ── straight-line fallback ──────────────────────────────
  if (pT < Const::pTmin || std::abs(By) < Const::Bmin) {
    const double dx = x1 - x0, dy = y1 - y0, dz = z1 - z0;
    return std::sqrt(dx * dx + dy * dy + dz * dz);
  }

  // ── helix geometry ──────────────────────────────────────
  // Bending radius in cm:  R[m] = pT / (k·|q|·|B|)
  const double R = 100.0 * pT / (Const::kMagFactor * std::abs(q) * std::abs(By));

  // Unit tangent in the x-z bending plane
  const double tx = px / pT;
  const double tz = pz / pT;

  // Centre of the helix circle (inward normal, sign from Lorentz force)
  const double sgn = (q * By > 0) ? +1.0 : -1.0;
  const double xc  = x0 + sgn * R * (-tz);
  const double zc  = z0 + sgn * R * ( tx);

  // Azimuthal angles of start and end points around the circle centre
  double phi0 = std::atan2(z0 - zc, x0 - xc);
  double phi1 = std::atan2(z1 - zc, x1 - xc);

  // Wrap Δφ to (-π, π]
  double dphi = phi1 - phi0;
  while (dphi >  M_PI) dphi -= 2.0 * M_PI;
  while (dphi < -M_PI) dphi += 2.0 * M_PI;

  // 3-D arc length = R·|Δφ| · (p / pT)
  return std::abs(R * dphi) * (p / pT); // cm
}

// ============================================================
/// Project h2 onto the β axis in momentum slices of width deltaP
/// and write each slice into a dedicated subdirectory.
// ============================================================
void saveBetaSlices(TH2F* h2, double deltaP)
{
  TDirectory* curDir  = gDirectory;
  TDirectory* sliceDir = curDir->mkdir(Form("BetaSlices_%s", h2->GetName()));
  sliceDir->cd();

  const TAxis* ax   = h2->GetXaxis();
  const double pMin = ax->GetXmin();
  const double pMax = ax->GetXmax();

  for (double pLo = pMin; pLo + 1.e-9 < pMax; pLo += deltaP) {
    const double pHi   = std::min(pLo + deltaP, pMax);
    const int    binLo = std::max(ax->FindBin(pLo + 1.e-9), 1);
    const int    binHi = std::min(ax->FindBin(pHi - 1.e-9), ax->GetNbins());

    TH1D* hSlice = h2->ProjectionY(
      Form("%s_p%.0f_%.0f", h2->GetName(), pLo * 1000, pHi * 1000),
      binLo, binHi);
    hSlice->SetTitle(Form("#beta   (%.2f < #it{p} < %.2f GeV/#it{c});#beta;counts",
                          pLo, pHi));
    hSlice->Write();
  }

  curDir->cd();
}

// ============================================================
/// Open a ROOT file and retrieve a named tree; abort on failure.
// ============================================================
TTree* openFileAndTree(const char* path, const char* treeName)
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

// ============================================================
/// Draw the expected β(p) curve for a particle of given mass.
// ============================================================
void drawBetaCurve(double mass, int color)
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

// ============================================================
/// Safely create a TEfficiency, checking consistency first.
/// Clamps numerator bins to the denominator where needed and
/// warns the user.
// ============================================================
TEfficiency* makeEfficiency(TH1* hNum, TH1* hDen,
                            const char* name, const char* title)
{
  // Clamp numerator to denominator bin-by-bin (required by TEfficiency)
  for (int b = 1; b <= hNum->GetNbinsX(); ++b) {
    if (hNum->GetBinContent(b) > hDen->GetBinContent(b)) {
      printf("  WARNING [%s] bin %d: numerator (%.0f) > denominator (%.0f) — clamping.\n",
             name, b, hNum->GetBinContent(b), hDen->GetBinContent(b));
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
  printf("  Created TEfficiency '%s'\n", name);
  return eff;
}

// ============================================================
//  Main macro
// ============================================================
void plotTOF(const char* dirSimu = ".", float pMaxIDProton = 3)
{
  gStyle->SetOptStat(0);
  gStyle->SetTitleSize(0.05, "XY");
  gStyle->SetLabelSize(0.045, "XY");
  gStyle->SetTitleOffset(1.1, "X");
  gStyle->SetTitleOffset(1.0, "Y");

  // ── open input files ────────────────────────────────────
  auto* tTracks = openFileAndTree(Form("%s/TracksTOFMatching.root", dirSimu), "tracksTOFMatching");
  auto* tMC     = openFileAndTree(Form("%s/MCKine.root",             dirSimu), "mckine");
  auto* tClus   = openFileAndTree(Form("%s/ClustersTOF.root",        dirSimu), "clustersTOF");
  auto* tVerTel = openFileAndTree(Form("%s/ClustersVerTel.root",     dirSimu), "clustersVerTel");
  if (!tTracks || !tMC || !tClus || !tVerTel) return;

  // ── branch wiring ───────────────────────────────────────
  std::vector<NA6PTrack>       tracks,      *tracksPtr    = &tracks;
  std::vector<TParticle>       *mcArr       = nullptr;
  std::vector<NA6PTOFCluster>  tofClus,     *tofClusPtr   = &tofClus;
  std::vector<NA6PVerTelCluster> verTelClus, *verTelClusPtr = &verTelClus;

  tTracks->SetBranchAddress("TOFMatching", &tracksPtr);
  tMC    ->SetBranchAddress("tracks",      &mcArr);
  tClus  ->SetBranchAddress("TOF",         &tofClusPtr);
  tVerTel->SetBranchAddress("VerTel",      &verTelClusPtr);

  const int nEv = std::min({(int)tTracks->GetEntries(),
                             (int)tMC    ->GetEntries(),
                             (int)tClus  ->GetEntries(),
                             (int)tVerTel->GetEntries()});
  printf("Processing %d events\n", nEv);

  // ── histogram declarations ──────────────────────────────
  // β vs p — all tracks with a TOF match
  auto* hBetaVsP = new TH2F("hBetaVsP",
    ";#it{p} (GeV/#it{c});#beta", 50, 0., 5., 300, 0.4, 1.10);

  // β vs p — MC-truth protons only
  auto* hBetaVsPprotons = new TH2F("hBetaVsPprotons",
    ";#it{p} (GeV/#it{c});#beta", 50, 0., 5., 300, 0.4, 1.10);

  // β vs p — per target
  TH2F* hBetaVsPtgt[Const::nTargets];
  for (int it = 0; it < Const::nTargets; ++it)
    hBetaVsPtgt[it] = new TH2F(
      Form("hBetaVsP_tgt%d", it),
      Form("target %d (z=%.1f cm);#it{p} (GeV/#it{c});#beta", it, Const::targetZ[it]),
      50, 0., 5., 300, 0.4, 1.10);

  // Generated proton phase space
  auto* hGenPVsYVsPt = new TH3D("hGenPVsYVsPt",
    ";p (GeV/c);y;p_{T} (GeV/c)", 50, 0., 50., 50, -2., 6., 50, 0., 5.);

  // Proton acceptance / efficiency numerators and denominators
  const int    nRapBins = 50;
  const double rapMin = -2., rapMax = 6.;
  auto* hProtonGen    = new TH1F("hProtonGen",    ";y;counts", nRapBins, rapMin, rapMax);
  auto* hProtonAcc    = new TH1F("hProtonAcc",    ";y;counts", nRapBins, rapMin, rapMax);
  auto* hProtonTOF    = new TH1F("hProtonTOF",    ";y;counts", nRapBins, rapMin, rapMax);
  auto* hProtonGenID  = new TH1F("hProtonGenID",  ";y;counts", nRapBins, rapMin, rapMax);
  auto* hProtonTOFID  = new TH1F("hProtonTOFID",  ";y;counts", nRapBins, rapMin, rapMax);

  const int    nPtBins = 50;
  const double ptMin = 0., ptMax = 5.;
  auto* hProtonGenPt   = new TH1F("hProtonGenPt",   ";p_{T} (GeV/c);counts", nPtBins, ptMin, ptMax);
  auto* hProtonAccPt   = new TH1F("hProtonAccPt",   ";p_{T} (GeV/c);counts", nPtBins, ptMin, ptMax);
  auto* hProtonTOFPt   = new TH1F("hProtonTOFPt",   ";p_{T} (GeV/c);counts", nPtBins, ptMin, ptMax);
  auto* hProtonGenIDPt = new TH1F("hProtonGenIDPt", ";p_{T} (GeV/c);counts", nPtBins, ptMin, ptMax);
  auto* hProtonTOFIDPt = new TH1F("hProtonTOFIDPt", ";p_{T} (GeV/c);counts", nPtBins, ptMin, ptMax);

  const int    nPBins = 50;
  const double pMinH = 0., pMaxH = 50.;
  auto* hProtonGenP = new TH1F("hProtonGenP", ";p (GeV/c);counts", nPBins, pMinH, pMaxH);
  auto* hProtonTOFP = new TH1F("hProtonTOFP", ";p (GeV/c);counts", nPBins, pMinH, pMaxH);

  const double mProton = TDatabasePDG::Instance()->GetParticle(2212)->Mass();

  // ── event loop ──────────────────────────────────────────
  for (int jEv = 0; jEv < nEv; ++jEv) {
    tMC    ->GetEvent(jEv);
    tTracks->GetEvent(jEv);
    tClus  ->GetEvent(jEv);
    tVerTel->GetEvent(jEv);

    // ── build lookup tables: particle ID → hit counts ──
    // VerTel: count hits per MC particle (need ≥ 4 for reconstructibility)
    std::unordered_map<int, int> verTelHits;
    for (const auto& clu : verTelClus)
      verTelHits[clu.getParticleID()]++;

    // TOF: count matched clusters per MC particle
    std::unordered_map<int, int> tofHits;
    for (const auto& clu : tofClus)
      tofHits[clu.getParticleID()]++;

    // ── primary vertex from first primary MC particle ──
    double xVertex = 0., yVertex = 0., zVertex = 0.;
    for (const auto& p : *mcArr) {
      if (p.IsPrimary()) {
        xVertex = p.Vx();
        yVertex = p.Vy();
        zVertex = p.Vz();
        break;
      }
    }

    // ── identify the closest target ───────────────────
    int iTgt = 0;
    double bestDz = std::abs(zVertex - Const::targetZ[0]);
    for (int it = 1; it < Const::nTargets; ++it) {
      double dz = std::abs(zVertex - Const::targetZ[it]);
      if (dz < bestDz) { bestDz = dz; iTgt = it; }
    }

    // ── generated proton denominators ─────────────────
    int nPart = -1;
    for (const auto& p : *mcArr) {
      nPart++;
      if (!p.IsPrimary()) continue;

      const int pdg = p.GetPdgCode();
      if (pdg != 2212 && pdg != -2212) continue;

      const double rap = p.Y();
      const double pt  = p.Pt();
      const double mom = p.P();

      hProtonGen  ->Fill(rap);
      hProtonGenPt->Fill(pt);
      hProtonGenP ->Fill(mom);
      hGenPVsYVsPt->Fill(mom, rap, pt);

      // IDentifiable protons: p below the PID momentum threshold
      if (mom < pMaxIDProton) {
        hProtonGenID  ->Fill(rap);
        hProtonGenIDPt->Fill(pt);

        // Geometrically reconstructible: ≥4 VerTel hits AND ≥1 TOF hit
        const bool reconstructible =
          (verTelHits.count(nPart) && verTelHits.at(nPart) >= Const::minVerTelHits) &&
          (tofHits.count(nPart)   && tofHits.at(nPart)    >= 1);

        if (reconstructible) {
          hProtonAcc  ->Fill(rap);
          hProtonAccPt->Fill(pt);
        }
      }
    }

    // ── reconstructed track loop ───────────────────────
    for (const auto& trk : tracks) {
      const double tof_s = trk.getTOF(); // raw value from tree [seconds, Geant4 convention]
      // MC-truth proton plots
      const int pid = trk.getParticleID();
      if (tof_s <= 0.) continue;         // no TOF match
      if (trk.getNVTHits() < Const::minVerTelHits) continue; // not reconstructible by geometry
      if (!(*mcArr)[pid].IsPrimary()) continue; // only consider primaries for this plot
      // Find the TOF cluster matched to this track by exact time comparison
      double xTOF = 0., yTOF = 0., zTOF = 0.;
      bool foundCluster = false;
      for (const auto& clu : tofClus) {
        if (std::abs(clu.getTime() - tof_s) < Const::tofTimeEps) {
          xTOF = clu.getXLab();
          yTOF = clu.getYLab();
          zTOF = clu.getZLab();
          foundCluster = true;
          break;
        }
      }
      if (!foundCluster) continue;

      // Retrieve momentum vector and charge for the helical path length
      double pxyz[3];
      trk.getPXYZ(pxyz);
      const int charge = trk.getCharge(); // must be +1 or -1

      // 3-D helical arc length from vertex to TOF cluster [cm]
      const double pathLen = getPathLength(
        xVertex, yVertex, zVertex,
        xTOF,    yTOF,    zTOF,
        pxyz[0], pxyz[1], pxyz[2],
        charge);

      // Convert Geant4 flight time from seconds to nanoseconds
      const double tof_ns = tof_s * Const::nsPerSec;

      const double beta = pathLen / (Const::cLight * tof_ns);

      const double rap = (*mcArr)[pid].Y();
      const double pt  = (*mcArr)[pid].Pt();
      const double mom = (*mcArr)[pid].P();

      hBetaVsP        ->Fill(mom, beta);
      hBetaVsPtgt[iTgt]->Fill(mom, beta);

      if (pid < 0 || pid >= (int)mcArr->size()) continue;

      const int pdg = (*mcArr)[pid].GetPdgCode();
      if (pdg != 2212 && pdg != -2212) continue;

      hProtonTOF  ->Fill(rap);
      hProtonTOFPt->Fill(pt);
      hProtonTOFP ->Fill(mom);
      hBetaVsPprotons->Fill(mom, beta);
      if (mom < pMaxIDProton) {
        // Mirror the denominator's geometric condition
        if (verTelHits.count(pid) && verTelHits.at(pid) >= Const::minVerTelHits) {
          hProtonTOFID  ->Fill(rap);
          hProtonTOFIDPt->Fill(pt);
        }
      }
    }
  }

  // ── output file ─────────────────────────────────────────
  auto* fOut = TFile::Open(Form("%s/TOFplots.root", dirSimu), "recreate");

  // Raw histograms
  for (auto* h : {hBetaVsP, hBetaVsPprotons})         h->Write();
  for (auto* h : {hProtonGen, hProtonAcc, hProtonTOF,
                  hProtonGenID, hProtonTOFID})          h->Write();
  for (auto* h : {hProtonGenPt, hProtonAccPt, hProtonTOFPt,
                  hProtonGenIDPt, hProtonTOFIDPt})      h->Write();
  for (auto* h : {hProtonGenP, hProtonTOFP})            h->Write();
  hGenPVsYVsPt->Write();
  for (int it = 0; it < Const::nTargets; ++it)
    hBetaVsPtgt[it]->Write();

  // β slices in momentum intervals
  saveBetaSlices(hBetaVsP, 0.2);
  saveBetaSlices(hBetaVsPprotons, 0.2);
  for (int it = 0; it < Const::nTargets; ++it)
    saveBetaSlices(hBetaVsPtgt[it], 0.2);

  // ── TEfficiency objects ──────────────────────────────────
  // Note: makeEfficiency() clamps numerator bins before creating the object,
  // and works on a copy so the original histograms are not modified.
  // Build copies for cases where clamping is needed.
  auto cloneH1 = [](TH1F* h, const char* newName) -> TH1F* {
    auto* c = (TH1F*)h->Clone(newName); c->SetDirectory(nullptr); return c;
  };

  auto* hTOFIDrap_clamped = cloneH1(hProtonTOFID,   "hProtonTOFID_clamped");
  auto* hTOFIDpt_clamped  = cloneH1(hProtonTOFIDPt, "hProtonTOFIDPt_clamped");

  auto* effRap        = makeEfficiency(hProtonTOF,        hProtonGen,
    "effProtonVsRapidity",
    "Proton TOF acceptance vs rapidity;y;#varepsilon");

  auto* effPt         = makeEfficiency(hProtonTOFPt,      hProtonGenPt,
    "effProtonVsPt",
    "Proton TOF acceptance vs p_{T};p_{T} (GeV/c);#varepsilon");

  auto* acceptanceRap = makeEfficiency(hProtonAcc,        hProtonGen,
    "accProtonVsRapidity",
    Form("Proton acceptance (4 VerTel hits, p < %.0f GeV/c) vs y;y;Acceptance",
         (double)pMaxIDProton));

  auto* acceptancePt  = makeEfficiency(hProtonAccPt,      hProtonGenPt,
    "accProtonVsPt",
    Form("Proton acceptance (4 VerTel hits, p < %.0f GeV/c) vs p_{T};p_{T} (GeV/c);Acceptance",
         (double)pMaxIDProton));

  auto* effTrkRap     = makeEfficiency(hTOFIDrap_clamped, hProtonAcc,
    "effTrackingVsRapidity",
    Form("Proton tracking efficiency (p < %.0f GeV/c) vs y;y;Tracking efficiency",
         (double)pMaxIDProton));

  auto* effTrkPt      = makeEfficiency(hTOFIDpt_clamped,  hProtonAccPt,
    "effTrackingVsPt",
    Form("Proton tracking efficiency (p < %.0f GeV/c) vs p_{T};p_{T} (GeV/c);Tracking efficiency",
         (double)pMaxIDProton));

  for (auto* eff : {effRap, effPt, acceptanceRap, acceptancePt, effTrkRap, effTrkPt})
    if (eff) eff->Write();

  // ── draw canvases ────────────────────────────────────────
  auto styleCanvas = [](TCanvas* c) {
    c->SetLeftMargin(0.12); c->SetBottomMargin(0.12);
    c->SetRightMargin(0.14); c->SetTopMargin(0.08);
  };

  // 1) β vs p — all tracks
  {
    auto* c1 = new TCanvas("cBetaVsP", "#beta vs p", 1800, 1400);
    styleCanvas(c1);
    c1->SetLogz();
    hBetaVsP->GetYaxis()->SetTitleOffset(0.9);
    hBetaVsP->GetYaxis()->SetRangeUser(0., 1.1);
    hBetaVsP->Draw("colz");
    drawBetaCurve(Const::mPion,  kRed);
    drawBetaCurve(Const::mKaon,  kBlue);
    drawBetaCurve(mProton,       kGreen + 2);
    c1->SaveAs(Form("%s/tof_beta_vs_p.png", dirSimu));
  }

  // Helper: draw a TEfficiency canvas and save it
  auto drawEff = [&](TEfficiency* eff, const char* cName, const char* cTitle,
                     const char* outFile) {
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

  drawEff(effRap,        "cEffRap",        "Proton acceptance vs y",
          "tof_proton_acc_rapidity.png");
  drawEff(effPt,         "cEffPt",         "Proton acceptance vs p_{T}",
          "tof_proton_acc_pt.png");
  drawEff(acceptanceRap, "cAccRap",        "Proton geometric acceptance vs y",
          "tof_proton_geoacc_rapidity.png");
  drawEff(acceptancePt,  "cAccPt",         "Proton geometric acceptance vs p_{T}",
          "tof_proton_geoacc_pt.png");
  drawEff(effTrkRap,     "cEffTrkRap",     "Proton tracking efficiency vs y",
          "tof_proton_efftrk_rapidity.png");
  drawEff(effTrkPt,      "cEffTrkPt",      "Proton tracking efficiency vs p_{T}",
          "tof_proton_efftrk_pt.png");

  // ── tidy up ──────────────────────────────────────────────
  fOut->Close();
  tTracks->GetCurrentFile()->Close();
  tMC    ->GetCurrentFile()->Close();
  tClus  ->GetCurrentFile()->Close();
  tVerTel->GetCurrentFile()->Close();

  printf("Done — plots saved in %s/\n", dirSimu);
}