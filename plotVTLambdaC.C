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
#include <TGeoManager.h>
#include <TStyle.h>
#include <TLorentzVector.h>
#include "NA6PTrack.h"
#include "MagneticField.h"
#include "NA6PVerTelHit.h"
#include "NA6PLayoutParam.h"
#include "NA6PFastTrackFitter.h"
#endif

// ============================================================
// MC truth helpers for Lambda_c+ -> p K- pi+  (PDG 4122)
// Decay products: proton (2212), K- (321, negative), pi+ (211, positive)
// ============================================================

bool isFromLambdaC(const TParticle& part, const std::vector<TParticle>& mcArr)
{
  int motherIndex = part.GetFirstMother();
  if (motherIndex < 0 || motherIndex >= static_cast<int>(mcArr.size())) return false;
  return TMath::Abs(mcArr[motherIndex].GetPdgCode()) == 4122;
}

bool isProtonFromLambdaC(const NA6PTrack& trk, const std::vector<TParticle>& mcArr)
{
  int pid = trk.getParticleID();
  if (pid < 0 || pid >= static_cast<int>(mcArr.size())) return false;
  const TParticle& part = mcArr[pid];
  if (TMath::Abs(part.GetPdgCode()) != 2212) return false; // not a proton
  return isFromLambdaC(part, mcArr);
}

bool isKaonFromLambdaC(const NA6PTrack& trk, const std::vector<TParticle>& mcArr)
{
  int pid = trk.getParticleID();
  if (pid < 0 || pid >= static_cast<int>(mcArr.size())) return false;
  const TParticle& part = mcArr[pid];
  if (TMath::Abs(part.GetPdgCode()) != 321) return false; // not a kaon
  return isFromLambdaC(part, mcArr);
}

bool isPionFromLambdaC(const NA6PTrack& trk, const std::vector<TParticle>& mcArr)
{
  int pid = trk.getParticleID();
  if (pid < 0 || pid >= static_cast<int>(mcArr.size())) return false;
  const TParticle& part = mcArr[pid];
  if (TMath::Abs(part.GetPdgCode()) != 211) return false; // not a pion
  return isFromLambdaC(part, mcArr);
}

// ============================================================
// Fit the invariant-mass sigma in each rapidity slice of a TH2D
// ============================================================
TH1D* getMassSigmaVsRapidity(TH2D* hMassVsRapidity,
                              const char* outName = "hMassSigmaVsRapidity",
                              int minEntriesPerBin = 20,
                              std::vector<TH1D*>* outProjections = nullptr)
{
  if (!hMassVsRapidity) return nullptr;

  const int    nRapidityBins = hMassVsRapidity->GetYaxis()->GetNbins();
  const double rapidityMin   = hMassVsRapidity->GetYaxis()->GetXmin();
  const double rapidityMax   = hMassVsRapidity->GetYaxis()->GetXmax();

  TH1D* hSigmaVsRapidity = new TH1D(outName,
    "#Lambda_{c} mass #sigma from Gaussian fit vs rapidity;Rapidity;#sigma_{M} (GeV/c^{2})",
    nRapidityBins, rapidityMin, rapidityMax);

  for (int iRap = 1; iRap <= nRapidityBins; ++iRap) {
    TH1D* hMassProj = hMassVsRapidity->ProjectionX(
      Form("%s_projYbin%d", outName, iRap), iRap, iRap);
    if (!hMassProj) continue;

    if (hMassProj->GetEntries() < minEntriesPerBin || hMassProj->GetMaximum() <= 0.) {
      if (outProjections) {
        TH1D* hProjClone = static_cast<TH1D*>(hMassProj->Clone(Form("%s_saved_projYbin%d", outName, iRap)));
        hProjClone->SetDirectory(nullptr);
        outProjections->push_back(hProjClone);
      }
      delete hMassProj; continue;
    }

    const int    maxBin     = hMassProj->GetMaximumBin();
    const double peak       = hMassProj->GetMean();
    const double rms        = hMassProj->GetRMS();
    double       halfWindow = (rms > 0.) ? 1.5 * rms : 0.03;
    if (halfWindow < 0.01) halfWindow = 0.01;

    double fitMin = TMath::Max(peak - halfWindow, hMassProj->GetXaxis()->GetXmin());
    double fitMax = TMath::Min(peak + halfWindow, hMassProj->GetXaxis()->GetXmax());

    TF1 fitFunc(Form("%s_fitYbin%d", outName, iRap), "gaus", fitMin, fitMax);
    fitFunc.SetParameters(hMassProj->GetMaximum(), peak, (rms > 0.) ? rms : 0.02);

    if (hMassProj->Fit(&fitFunc, "QR") == 0) {
      hSigmaVsRapidity->SetBinContent(iRap, TMath::Abs(fitFunc.GetParameter(2)));
      hSigmaVsRapidity->SetBinError(iRap, fitFunc.GetParError(2));
    }

    if (outProjections) {
      TH1D* hProjClone = static_cast<TH1D*>(hMassProj->Clone(Form("%s_saved_projYbin%d", outName, iRap)));
      hProjClone->SetDirectory(nullptr);
      outProjections->push_back(hProjClone);
    }
    delete hMassProj;
  }
  return hSigmaVsRapidity;
}

// ============================================================
// Main macro
// ============================================================
void plotVTLambdaC(const char* dirSimu = ".", int NRecoClusters = 4)
{
  // --- Initialization & field loading ---
  auto& param = NA6PLayoutParam::Instance();
  param.updateFromFile(Form("%s/na6pLayout.ini", dirSimu), "", true);

  auto magField = new MagneticField();
  magField->loadField();
  magField->setAsGlobalField();

  NA6PFastTrackFitter fitter;
  fitter.loadGeometry();
  fitter.enableMaterialCorrections();

  // --- File loading ---
  TFile* ft = TFile::Open(Form("%s/TracksVerTel.root", dirSimu));
  if (!ft || ft->IsZombie()) { printf("ERROR: cannot open TracksVerTel.root\n"); return; }
  TTree* trTree = (TTree*)ft->Get("tracksVerTel");
  std::vector<NA6PTrack>* trArr = nullptr;
  trTree->SetBranchAddress("VerTel", &trArr);

  TFile* fk = TFile::Open(Form("%s/MCKine.root", dirSimu));
  if (!fk || fk->IsZombie()) { printf("ERROR: cannot open MCKine.root\n"); return; }
  TTree* mcTree = (TTree*)fk->Get("mckine");
  std::vector<TParticle>* mcArr = nullptr;
  mcTree->SetBranchAddress("tracks", &mcArr);

  TFile* fh = TFile::Open(Form("%s/HitsVerTel.root", dirSimu));
  TTree* th  = (fh) ? (TTree*)fh->Get("hitsVerTel") : nullptr;

  // PDG masses (GeV/c^2)
  const float mProton = 0.93827; // proton
  const float mKaon   = 0.49367; // K
  const float mPion   = 0.13957; // pi
  const float mLambdaC = 2.28646; // Lambda_c+

  const int   nbinsM = 100;
  const int   nbinsY = 20;
  const double massWin = 0.20; // +/- window around Lambda_c mass

  // --- Histogram definitions ---
  // Reconstructed
  TH1D* hLcMass      = new TH1D("hLcMass",      "#Lambda_{c} invariant mass;M (GeV/c^{2});Counts",
                                 nbinsM, mLambdaC - massWin, mLambdaC + massWin);
  TH1D* hLcRapidity  = new TH1D("hLcRapidity",  "#Lambda_{c} rapidity;y;Counts",                   100, 0., 5.);
  TH1D* hLcPt        = new TH1D("hLcPt",        "#Lambda_{c} p_{T};p_{T} (GeV/c);Counts",          100, 0., 8.);
  TH2D* hLcMassVsPt  = new TH2D("hLcMassVsPt",  "#Lambda_{c} mass vs p_{T};M (GeV/c^{2});p_{T} (GeV/c)",
                                  nbinsM, mLambdaC - massWin, mLambdaC + massWin, 100, 0., 8.);
  TH2D* hLcMassVsY   = new TH2D("hLcMassVsY",   "#Lambda_{c} mass vs rapidity;M (GeV/c^{2});y",
                                  nbinsM, mLambdaC - massWin, mLambdaC + massWin, nbinsY, 0., 5.);

  // Gen-level / resolution
  TH1D* hLcMassGen       = new TH1D("hLcMassGen",       "#Lambda_{c} inv. mass (gen);M (GeV/c^{2});Counts",
                                     100, mLambdaC - massWin, mLambdaC + massWin);
  TH1D* hLcRapidityMC    = new TH1D("hLcRapidityMC",    "#Lambda_{c} rapidity (MC truth);y;Counts",    100, 0., 5.);
  TH2D* hLcYRecVsGen     = new TH2D("hLcYRecVsGen",     "y rec vs gen;y_{rec};y_{gen}",                100, 0., 5., 100, 0., 5.);
  TH2D* hLcPtRecVsGen    = new TH2D("hLcPtRecVsGen",    "p_{T} rec vs gen;p_{T,rec};p_{T,gen}",        100, 0., 8., 100, 0., 8.);
  TH1D* hLcDeltaY        = new TH1D("hLcDeltaY",        "#Delta y (rec-gen);#Delta y;Counts",           100, -1., 1.);
  TH1D* hLcDeltaPt       = new TH1D("hLcDeltaPt",       "#Delta p_{T} (rec-gen);#Delta p_{T} (GeV/c);Counts", 100, -6., 6.);

  // Track-quality monitoring
  TH1D* hRapidityAllTrk  = new TH1D("hRapidityAllTrk", "all-track rapidity;y;Counts", 100, -2., 7.);
  TH1D* hRapidityNeg     = new TH1D("hRapidityNeg",    "negative-track rapidity;y;Counts", 100, -2., 7.);
  TH1D* hRapidityPos     = new TH1D("hRapidityPos",    "positive-track rapidity;y;Counts", 100, -2., 7.);

  int nEv = trTree->GetEntries();
  printf("Number of events = %d\n", nEv);

  for (int jEv = 0; jEv < nEv; jEv++) {
    mcTree->GetEvent(jEv);
    trTree->GetEvent(jEv);
    if (th) th->GetEvent(jEv);

    int nPart   = (mcArr) ? (int)mcArr->size() : 0;
    int nTracks = (trArr) ? (int)trArr->size() : 0;
    printf("\nEvent %d: particles=%d  tracks=%d\n", jEv, nPart, nTracks);

    // Find decay vertex z from any Lambda_c decay product
    double zOrig = 0.;
    for (int jp = 0; jp < nPart; jp++) {
      const TParticle& curPart = mcArr->at(jp);
      if (isFromLambdaC(curPart, *mcArr)) {
        zOrig = curPart.Vz();
        printf("  zOrig (Lambda_c vertex): %.2f cm\n", zOrig);
        break;
      }
    }

    // Fill MC-truth Lambda_c rapidity
    for (int jp = 0; jp < nPart; jp++) {
      const TParticle& curPart = mcArr->at(jp);
      if (TMath::Abs(curPart.GetPdgCode()) == 4122) {
        double px = curPart.Px(), py = curPart.Py(), pz = curPart.Pz();
        double p  = TMath::Sqrt(px*px + py*py + pz*pz);
        if (p > TMath::Abs(pz)) // guard against log(0)
          hLcRapidityMC->Fill(curPart.Y());
      }
    }

    // Pre-filter tracks: require minimum cluster count and valid MC label
    std::vector<int> filteredIndices;
    for (int jTr = 0; jTr < nTracks; ++jTr) {
      NA6PTrack& tr = trArr->at(jTr);
      if (tr.getNHits() < NRecoClusters)               continue;
      if (tr.getParticleID() < 0 ||
          tr.getParticleID() >= nPart)                 continue;

      double pxyz[3]; tr.getPXYZ(pxyz);
      double p = TMath::Sqrt(pxyz[0]*pxyz[0] + pxyz[1]*pxyz[1] + pxyz[2]*pxyz[2]);
      if (p > TMath::Abs(pxyz[2])) {
        double rap = 0.5 * TMath::Log((p + pxyz[2]) / (p - pxyz[2]));
        hRapidityAllTrk->Fill(rap);
        if (tr.getCharge() < 0) hRapidityNeg->Fill(rap);
        if (tr.getCharge() > 0) hRapidityPos->Fill(rap);
      }
      filteredIndices.push_back(jTr);
    }

    // -------------------------------------------------------
    // 3-body combinatorics: proton (p>0) + kaon (q<0) + pion (q>0)
    // Lambda_c+ -> p  K-  pi+
    //              2212 -321 211
    // -------------------------------------------------------
    const int nFilt = (int)filteredIndices.size();
    if (nFilt < 3) continue;

    for (int i = 0; i < nFilt; ++i) {
      NA6PTrack trP = trArr->at(filteredIndices[i]);
      if (!isProtonFromLambdaC(trP, *mcArr)) continue;
      // proton must be positive
      if (trP.getCharge() <= 0) continue;

      fitter.propagateToZ(&trP, zOrig);
      double momP[3]; trP.getPXYZ(momP);
      double RecPp = trP.getP();

      int mcLabelP  = trP.getParticleID();
      bool hasMC_P  = (mcLabelP >= 0 && mcLabelP < nPart);
      TLorentzVector genP;
      if (hasMC_P) {
        TParticle& pp = mcArr->at(mcLabelP);
        genP.SetPxPyPzE(pp.Px(), pp.Py(), pp.Pz(), pp.Energy());
      }

      for (int j = 0; j < nFilt; ++j) {
        if (j == i) continue;
        NA6PTrack trK = trArr->at(filteredIndices[j]);
        if (!isKaonFromLambdaC(trK, *mcArr)) continue;
        // kaon must be negative
        if (trK.getCharge() >= 0) continue;

        fitter.propagateToZ(&trK, zOrig);
        double momK[3]; trK.getPXYZ(momK);
        double RecPk = trK.getP();

        int mcLabelK = trK.getParticleID();
        bool hasMC_K = (mcLabelK >= 0 && mcLabelK < nPart);
        TLorentzVector genK;
        if (hasMC_K) {
          TParticle& pk = mcArr->at(mcLabelK);
          genK.SetPxPyPzE(pk.Px(), pk.Py(), pk.Pz(), pk.Energy());
        }

        for (int k = 0; k < nFilt; ++k) {
          if (k == i || k == j) continue;
          NA6PTrack trPi = trArr->at(filteredIndices[k]);
          if (!isPionFromLambdaC(trPi, *mcArr)) continue;
          // pion must be positive
          if (trPi.getCharge() <= 0) continue;

          fitter.propagateToZ(&trPi, zOrig);
          double momPi[3]; trPi.getPXYZ(momPi);
          double RecPpi = trPi.getP();

          // Build Lorentz vectors with correct mass hypotheses
          TLorentzVector lv_p, lv_K, lv_pi;
          lv_p .SetPxPyPzE(momP [0], momP [1], momP [2],
                            TMath::Sqrt(RecPp  * RecPp  + mProton * mProton));
          lv_K .SetPxPyPzE(momK [0], momK [1], momK [2],
                            TMath::Sqrt(RecPk  * RecPk  + mKaon   * mKaon));
          lv_pi.SetPxPyPzE(momPi[0], momPi[1], momPi[2],
                            TMath::Sqrt(RecPpi * RecPpi + mPion   * mPion));

          TLorentzVector lc = lv_p + lv_K + lv_pi;

          hLcMass    ->Fill(lc.M());
          hLcRapidity->Fill(lc.Rapidity());
          hLcPt      ->Fill(lc.Pt());
          hLcMassVsPt->Fill(lc.M(), lc.Pt());
          hLcMassVsY ->Fill(lc.M(), lc.Rapidity());

          // Fill gen-level quantities if all three prongs have MC truth
          int mcLabelPi = trPi.getParticleID();
          bool hasMC_Pi = (mcLabelPi >= 0 && mcLabelPi < nPart);
          if (hasMC_P && hasMC_K && hasMC_Pi) {
            TParticle& ppi = mcArr->at(mcLabelPi);
            TLorentzVector genPi;
            genPi.SetPxPyPzE(ppi.Px(), ppi.Py(), ppi.Pz(), ppi.Energy());
            TLorentzVector lcGen = genP + genK + genPi;

            hLcMassGen  ->Fill(lcGen.M());
            hLcYRecVsGen ->Fill(lc.Rapidity(),   lcGen.Rapidity());
            hLcPtRecVsGen->Fill(lc.Pt(),         lcGen.Pt());
            hLcDeltaY   ->Fill(lc.Rapidity()  - lcGen.Rapidity());
            hLcDeltaPt  ->Fill(lc.Pt()        - lcGen.Pt());
          }
        } // pi loop
      }   // K loop
    }     // p loop
  }       // event loop

  // --- Plotting ---
  gStyle->SetOptFit(1111);

  // Invariant mass with Gaussian fit
  TCanvas* cMass = new TCanvas("cMass", "Lambda_c mass", 1200, 800);
  hLcMass->SetMarkerColor(kAzure + 8);
  hLcMass->SetLineColor(kBlue + 1);
  hLcMass->SetMarkerStyle(20);
  hLcMass->SetTitle("#Lambda_{c}^{+} Invariant Mass (p K^{-} #pi^{+});M (GeV/c^{2});Counts");
  hLcMass->Fit("gaus", "", "", mLambdaC - 0.05, mLambdaC + 0.05);
  hLcMass->Draw("EP");
  cMass->SaveAs(Form("%s/LambdaC_mass.png", dirSimu));

  // Mass sigma vs rapidity
  std::vector<TH1D*> massProjectionsVsY;
  TH1D* hMassSigmaVsY = getMassSigmaVsRapidity(hLcMassVsY, "hLcMassSigmaVsY", 20, &massProjectionsVsY);
  TCanvas* cSigmaY = new TCanvas("cSigmaY", "Mass sigma vs rapidity", 1200, 800);
  if (hMassSigmaVsY) {
    hMassSigmaVsY->SetMarkerColor(kRed + 1);
    hMassSigmaVsY->SetLineColor(kRed + 1);
    hMassSigmaVsY->GetYaxis()->SetRangeUser(0., 0.1);
    hMassSigmaVsY->Draw("EP");
  }
  cSigmaY->SaveAs(Form("%s/LambdaC_mass_sigma_vs_rapidity.png", dirSimu));

  // Rapidity and pT
  TCanvas* cKin = new TCanvas("cKin", "Lambda_c kinematics", 1600, 800);
  cKin->Divide(2, 1);
  cKin->cd(1); hLcRapidity->Draw(); hLcRapidityMC->SetLineColor(kRed); hLcRapidityMC->Draw("SAME");
  TLegend* leg = new TLegend(0.6, 0.7, 0.9, 0.9);
  leg->AddEntry(hLcRapidity,   "Reconstructed", "l");
  leg->AddEntry(hLcRapidityMC, "MC truth",       "l");
  leg->Draw();
  cKin->cd(2); hLcPt->Draw();
  cKin->SaveAs(Form("%s/LambdaC_kinematics.png", dirSimu));

  // Resolution
  TCanvas* cRes = new TCanvas("cRes", "Resolution", 1600, 800);
  cRes->Divide(2, 1);
  cRes->cd(1); hLcDeltaY ->Draw();
  cRes->cd(2); hLcDeltaPt->Draw();
  cRes->SaveAs(Form("%s/LambdaC_resolution.png", dirSimu));

  // --- Save to ROOT file ---
  TFile fOut(Form("%s/LambdaC_analysis.root", dirSimu), "RECREATE");
  hLcMass          ->Write();
  hLcRapidity      ->Write();
  hLcPt            ->Write();
  hLcMassVsPt      ->Write();
  hLcMassVsY       ->Write();
  hLcMassGen       ->Write();
  hLcRapidityMC    ->Write();
  hLcYRecVsGen     ->Write();
  hLcPtRecVsGen    ->Write();
  hLcDeltaY        ->Write();
  hLcDeltaPt       ->Write();
  hRapidityAllTrk  ->Write();
  hRapidityNeg     ->Write();
  hRapidityPos     ->Write();
  if (hMassSigmaVsY) hMassSigmaVsY->Write();
  if (!massProjectionsVsY.empty()) {
    TDirectory* dProj = fOut.mkdir("MassProjectionsVsRapidity");
    if (dProj) {
      dProj->cd();
      for (auto* hProj : massProjectionsVsY) {
        if (hProj) hProj->Write();
      }
      fOut.cd();
    }
  }
  fOut.Close();

  for (auto* hProj : massProjectionsVsY) {
    delete hProj;
  }
  printf("\nDone. Output: %s/LambdaC_analysis.root\n", dirSimu);
}