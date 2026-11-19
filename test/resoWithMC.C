// resoWithMC.C - Resolution analysis with MC truth matching
#include <unordered_map>
#include <algorithm>
#include <TFile.h>
#include <TTree.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TParticle.h>
#include "NA6PTrack.h"

void resoWithMC(
  bool useMC = true,
  const char* dir = ".",
  const char* dirMC = "../fullev")
{
  gStyle->SetOptStat(0);
  gStyle->SetPalette(kRainBow);

  // ---- Open files with error checking
  auto openFile = [](const char* name) -> TFile* {
    TFile* f = TFile::Open(name);
    if (!f || f->IsZombie()) {
      printf("Error: Cannot open file %s\n", name);
      return nullptr;
    }
    return f;
  };

  TFile* fVT = openFile(Form("%s/TracksVerTelMatching%s.root", dir, useMC ? "_MC" : ""));
  TFile* fMS = openFile(Form("%s/TracksMuonSpecMatching%s.root", dir, useMC ? "_MC" : ""));
  TFile* fM = openFile(Form("%s/TracksMatching%s.root", dir, useMC ? "_MC" : ""));
  TFile* fMC = openFile(Form("%s/MCKine.root", dirMC));

  if (!fVT || !fMS || !fM || !fMC) {
    printf("Failed to open required files. Exiting.\n");
    return;
  }

  // ---- Get trees with error checking
  auto getTree = [](TFile* f, const char* name) -> TTree* {
    TTree* t = (TTree*)f->Get(name);
    if (!t) {
      printf("Error: Cannot find tree %s\n", name);
    }
    return t;
  };

  TTree* tVT = getTree(fVT, "tracksMatching");
  TTree* tMS = getTree(fMS, "tracksMatching");
  TTree* tM = getTree(fM, "tracksMatching");
  TTree* tMC = getTree(fMC, "mckine");

  if (!tVT || !tMS || !tM || !tMC) {
    printf("Failed to retrieve required trees. Exiting.\n");
    return;
  }

  std::vector<NA6PTrack>* vVT = nullptr;
  std::vector<NA6PTrack>* vMS = nullptr;
  std::vector<NA6PTrack>* vM = nullptr;
  std::vector<TParticle>* vMC = nullptr;

  tVT->SetBranchAddress("Matching", &vVT);
  tMS->SetBranchAddress("Matching", &vMS);
  tM->SetBranchAddress("Matching", &vM);
  tMC->SetBranchAddress("tracks", &vMC);

  float maxP = 20.0; // max momentum for histograms
  float minP = 0.0;  // min momentum for histograms
  int nPBins = 10;    // number of momentum bins

  // Fill rejection histogram labels
  auto setRejLabels = [](TH1F* hRej) {
    hRej->GetXaxis()->SetBinLabel(1, "Fake track");
    hRej->GetXaxis()->SetBinLabel(2, "No MC match");
    hRej->GetXaxis()->SetBinLabel(3, "Accepted");
  };

  int nLabels = 3;

  // VT histograms
  TH1F* hDx_VT = new TH1F("hDx_VT", "#Delta x (cm);#Delta x;Entries", 200, -0.02, 0.02);
  TH1F* hDy_VT = new TH1F("hDy_VT", "#Delta y (cm);#Delta y;Entries", 200, -0.02, 0.02);
  TH1F* hDz_VT = new TH1F("hDz_VT", "#Delta z (cm);#Delta z;Entries", 200, -0.1, 0.1);
  TH1F* hDpx_VT = new TH1F("hDpx_VT", "#Delta p_{x} (GeV/c);#Delta p_{x};Entries", 200, -1, 1);
  TH1F* hDpy_VT = new TH1F("hDpy_VT", "#Delta p_{y} (GeV/c);#Delta p_{y};Entries", 200, -1, 1);
  TH1F* hDpz_VT = new TH1F("hDpz_VT", "#Delta p_{z} (GeV/c);#Delta p_{z};Entries", 200, -2, 5);
  TH1F* hDp_VT = new TH1F("hDp_VT", "#Delta p (GeV/c);#Delta p;Entries", 200, -2, 2);
  TH2F* hPVsP_VT = new TH2F("hPVsP_VT", ";p_{rec} (GeV/c);#Deltap (GeV/c)", nPBins, minP, maxP, 200, -5, 5);
  TH2F* hDxVsP_VT = new TH2F("hDxVsP_VT", ";p (GeV/c);#Delta x (cm)", nPBins, minP, maxP, 200, -0.02, 0.02);
  TH2F* hDyVsP_VT = new TH2F("hDyVsP_VT", ";p (GeV/c);#Delta y (cm)", nPBins, minP, maxP, 200, -0.02, 0.02);
  TH1F* hP_VT = new TH1F("hP_VT", ";#it{p} (GeV/#it{c});#sigma_{#it{p}} (GeV/#it{c})", nPBins, minP, maxP);
  TH1F* hResX_VT = new TH1F("hResX_VT", ";#it{p} (GeV/#it{c});#sigma_{x} (#mum)", nPBins, minP, maxP);
  TH1F* hResY_VT = new TH1F("hResY_VT", ";#it{p} (GeV/#it{c});#sigma_{y} (#mum)", nPBins, minP, maxP);
  TH1F* hRej_VT = new TH1F("hRej_VT", ";;Counts", nLabels, -0.5, nLabels - 0.5);
  TH1F* hPID_VT = new TH1F("hPID_VT", ";Particle ID;Entries", 10000, -0.5, 9999.5);  
  TH1F* hPIDmuJpsi_VT = new TH1F("hPIDmuJpsi_VT", ";Particle ID;Entries", 10000, -0.5, 9999.5);  

  setRejLabels(hRej_VT);
  // MS histograms
  TH1F* hDx_MS = new TH1F("hDx_MS", "#Delta x (cm);#Delta x;Entries", 200, -5, 5);
  TH1F* hDy_MS = new TH1F("hDy_MS", "#Delta y (cm);#Delta y;Entries", 200, -5, 5);
  TH1F* hDz_MS = new TH1F("hDz_MS", "#Delta z (cm);#Delta z;Entries", 200, -5, 5);
  TH1F* hDpx_MS = new TH1F("hDpx_MS", "#Delta p_{x} (GeV/c);#Delta p_{x};Entries", 200, -1, 1);
  TH1F* hDpy_MS = new TH1F("hDpy_MS", "#Delta p_{y} (GeV/c);#Delta p_{y};Entries", 200, -1, 1);
  TH1F* hDpz_MS = new TH1F("hDpz_MS", "#Delta p_{z} (GeV/c);#Delta p_{z};Entries", 200, -2, 5);
  TH1F* hDp_MS = new TH1F("hDp_MS", "#Delta p (GeV/c);#Delta p;Entries", 200, -5, 5);
  TH2F* hPVsP_MS = new TH2F("hPVsP_MS", ";p_{rec} (GeV/c);#Deltap (GeV/c)", nPBins, minP, maxP, 200, -5, 5);
  TH2F* hDxVsP_MS = new TH2F("hDxVsP_MS", ";p (GeV/c);#Delta x (cm)", nPBins, minP, maxP, 200, -5, 5);
  TH2F* hDyVsP_MS = new TH2F("hDyVsP_MS", ";p (GeV/c);#Delta y (cm)", nPBins, minP, maxP, 200, -5, 5);
  TH1F* hP_MS = new TH1F("hP_MS", ";#it{p} (GeV/#it{c});#sigma_{#it{p}} (GeV/#it{c})", nPBins, minP, maxP);
  TH1F* hResX_MS = new TH1F("hResX_MS", ";#it{p} (GeV/#it{c});#sigma_{x} (#mum)", nPBins, minP, maxP);
  TH1F* hResY_MS = new TH1F("hResY_MS", ";#it{p} (GeV/#it{c});#sigma_{y} (#mum)", nPBins, minP, maxP);
  TH1F* hRej_MS = new TH1F("hRej_MS", ";;Counts", nLabels, -0.5, nLabels - 0.5);
  TH1F* hPID_MS = new TH1F("hPID_MS", ";Particle ID;Entries", 10000, -0.5, 9999.5);
  TH1F* hPIDmuJpsi_MS = new TH1F("hPIDmuJpsi_MS", ";Particle ID;Entries", 10000, -0.5, 9999.5);

  setRejLabels(hRej_MS);
  // MT (Matched) histograms
  TH1F* hDx_MT = new TH1F("hDx_MT", "#Delta x (cm);#Delta x;Entries", 200, -0.02, 0.02);
  TH1F* hDy_MT = new TH1F("hDy_MT", "#Delta y (cm);#Delta y;Entries", 200, -0.02, 0.02);
  TH1F* hDz_MT = new TH1F("hDz_MT", "#Delta z (cm);#Delta z;Entries", 200, -0.1, 0.1);
  TH1F* hDpx_MT = new TH1F("hDpx_MT", "#Delta p_{x} (GeV/c);#Delta p_{x};Entries", 200, -1, 1);
  TH1F* hDpy_MT = new TH1F("hDpy_MT", "#Delta p_{y} (GeV/c);#Delta p_{y};Entries", 200, -1, 1);
  TH1F* hDpz_MT = new TH1F("hDpz_MT", "#Delta p_{z} (GeV/c);#Delta p_{z};Entries", 200, -2, 5);
  TH1F* hDp_MT = new TH1F("hDp_MT", "#Delta p (GeV/c);#Delta p;Entries", 200, -2, 2);
  TH2F* hPVsP_MT = new TH2F("hPVsP_MT", ";p_{rec} (GeV/c);#Deltap (GeV/c)", nPBins, minP, maxP, 200, -5, 5);
  TH2F* hDxVsP_MT = new TH2F("hDxVsP_MT", ";p (GeV/c);#Delta x (cm)", nPBins, minP, maxP, 200, -0.02, 0.02);
  TH2F* hDyVsP_MT = new TH2F("hDyVsP_MT", ";p (GeV/c);#Delta y (cm)", nPBins, minP, maxP, 200, -0.02, 0.02);
  TH1F* hP_MT = new TH1F("hP_MT", ";#it{p} (GeV/#it{c});#sigma_{#it{p}} (GeV/#it{c})", nPBins, minP, maxP);
  TH1F* hResX_MT = new TH1F("hResX_MT", ";#it{p} (GeV/#it{c});#sigma_{x} (#mum)", nPBins, minP, maxP);
  TH1F* hResY_MT = new TH1F("hResY_MT", ";#it{p} (GeV/#it{c});#sigma_{y} (#mum)", nPBins, minP, maxP);
  TH1F* hRej_MT = new TH1F("hRej_MT", ";;Counts", nLabels, -0.5, nLabels - 0.5);
  TH1F* hPID_MT = new TH1F("hPID_MT", ";Particle ID;Entries", 10000, -0.5, 9999.5);
  TH1F* hPIDmuJpsi_MT = new TH1F("hPIDmuJpsi_MT", ";Particle ID;Entries", 10000, -0.5, 9999.5);
  
  setRejLabels(hRej_MT);
  // ---- Fill residuals with cleaner lambda
  auto fillResiduals = [&](TTree* tTrk, std::vector<NA6PTrack>*& vv,
                           TH1F* hdx, TH1F* hdy, TH1F* hdz,
                           TH1F* hdpx, TH1F* hdpy, TH1F* hdpz, TH1F* hdp,
                           TH2F* hPVsP, TH2F* hDxVsP, TH2F* hDyVsP, TH1F* hRej, TH1F* hPID, TH1F* hPIDmuJpsi) {
    Long64_t nEv = tTrk->GetEntries();
    Long64_t nMC = tMC->GetEntries();

    if (nEv != nMC) {
      printf("Warning: nEvents mismatch tracks=%lld mc=%lld\n", nEv, nMC);
    }

    Long64_t n = std::min(nEv, nMC);
    for (Long64_t iev = 0; iev < n; ++iev) {
      if ((iev + 1) % 100 == 0)
        printf("Processing event %lld / %lld\r", iev + 1, n);

      tTrk->GetEntry(iev);
      tMC->GetEntry(iev);

      // Loop over reconstructed tracks
      for (const auto& tr : *vv) {
        int pid = tr.getParticleID();        
        hPID->Fill(pid);

        if (pid < 0) {
          hRej->Fill(0); // fake track
          continue;      // Skip mismatched tracks
        }


        const auto& mc = (*vMC)[pid];
        truth[pid] = &p;
        int motherPdg = -1;
        
        if (p.GetMother(0) >= 0 && p.GetMother(0) < vMC->size()) {
          motherPdg = (*vMC)[p.GetMother(0)].GetPdgCode();
        }
        if (std::abs(p.GetPdgCode()) == 13 && motherPdg == 443){
            const TParticle* mc = it->second;

          double xyz[3], pxyz[3];
          tr.getXYZ(xyz);
          tr.getPXYZ(pxyz);
          double pmag = tr.getP();

          double dx = xyz[0] - mc->Vx();
          double dy = xyz[1] - mc->Vy();
          double dz = xyz[2] - mc->Vz();
          double dpx = pxyz[0] - mc->Px();
          double dpy = pxyz[1] - mc->Py();
          double dpz = pxyz[2] - mc->Pz();
          double dp = pmag - mc->P();

          hdx->Fill(dx);
          hdy->Fill(dy);
          hdz->Fill(dz);
          hdpx->Fill(dpx);
          hdpy->Fill(dpy);
          hdpz->Fill(dpz);
          hdp->Fill(dp);
          hPVsP->Fill(pmag, dp);
          hDxVsP->Fill(pmag, dx);
          hDyVsP->Fill(pmag, dy);
          hRej->Fill(2); // accepted
        }

      }
    }
    printf("\n");
  };

  // ---- Fill all three detectors
  printf("Processing VerTel tracks...\n");
  fillResiduals(tVT, vVT, hDx_VT, hDy_VT, hDz_VT, hDpx_VT, hDpy_VT, hDpz_VT, hDp_VT, hPVsP_VT, hDxVsP_VT, hDyVsP_VT, hRej_VT, hPID_VT, hPIDmuJpsi_VT);

  printf("Processing MuonSpec tracks...\n");
  fillResiduals(tMS, vMS, hDx_MS, hDy_MS, hDz_MS, hDpx_MS, hDpy_MS, hDpz_MS, hDp_MS, hPVsP_MS, hDxVsP_MS, hDyVsP_MS, hRej_MS, hPID_MS, hPIDmuJpsi_MS);

  printf("Processing Matched tracks...\n");
  fillResiduals(tM, vM, hDx_MT, hDy_MT, hDz_MT, hDpx_MT, hDpy_MT, hDpz_MT, hDp_MT, hPVsP_MT, hDxVsP_MT, hDyVsP_MT, hRej_MT, hPID_MT, hPIDmuJpsi_MT);
  // ---- Compute resolutions from 2D histograms
  auto computeResolution = [](TH2F* hPVsP, TH1F* hp, double scale = 1.0) {
    for (int i = 1; i <= hPVsP->GetNbinsX(); ++i) {
      TH1D* py = hPVsP->ProjectionY(Form("py_%s_%d", hPVsP->GetName(), i), i, i);
      if (py->GetEntries() > 10) { // Require minimum statistics
        TF1 f("gaus", "gaus", py->GetMean() - 3 * py->GetRMS(), py->GetMean() + 3 * py->GetRMS());
        py->Fit(&f, "Q");
        hp->SetBinContent(i, f.GetParameter(2) * scale);
        hp->SetBinError(i, f.GetParError(2) * scale);
      }
      delete py;
    }
  };

  printf("Computing resolutions...\n");
  computeResolution(hPVsP_VT, hP_VT);
  computeResolution(hPVsP_MS, hP_MS);
  computeResolution(hPVsP_MT, hP_MT);

  computeResolution(hDxVsP_VT, hResX_VT, 10000); // Convert cm to um
  computeResolution(hDyVsP_VT, hResY_VT, 10000);
  computeResolution(hDxVsP_MS, hResX_MS, 10000);
  computeResolution(hDyVsP_MS, hResY_MS, 10000);
  computeResolution(hDxVsP_MT, hResX_MT, 10000);
  computeResolution(hDyVsP_MT, hResY_MT, 10000);

  // ---- Setup histogram cosmetics
  auto setStyle = [](TH1F* h, Color_t color, Style_t marker) {
    h->SetLineColor(color);
    h->SetMarkerColor(color);
    h->SetMarkerStyle(marker);
    h->SetLineWidth(2);
    h->GetXaxis()->SetTitleSize(0.05);
    h->GetYaxis()->SetTitleSize(0.05);
    h->GetXaxis()->SetLabelSize(0.045);
    h->GetYaxis()->SetLabelSize(0.045);
  };

  setStyle(hP_VT, kRed, 20);
  setStyle(hP_MT, kGreen + 1, 21);
  setStyle(hResX_VT, kRed, 20);
  setStyle(hResX_MT, kGreen + 1, 21);
  setStyle(hResY_VT, kRed, 20);
  setStyle(hResY_MT, kGreen + 1, 21);

  // ---- Create summary plot
  printf("Creating summary plot...\n");
  TCanvas c1("c1", "Resolution Summary", 1200, 1200 * 3);
  c1.Divide(1, 3);

  auto divideByMomentum = [](TH1F* hist) {
    for (int i = 1; i <= hist->GetNbinsX(); ++i) {
      double pres = hist->GetBinCenter(i);
      if (pres != 0) {
        double val = hist->GetBinContent(i);
        double err = hist->GetBinError(i);
        hist->SetBinContent(i, val / pres);
        hist->SetBinError(i, err / pres);
      }
    }
  };
  divideByMomentum(hP_VT);
  divideByMomentum(hP_MT);
  divideByMomentum(hP_MS);
  c1.cd(1);
  gPad->SetLeftMargin(0.14);
  gPad->SetRightMargin(0.05);
  gPad->SetTopMargin(0.08);
  gPad->SetBottomMargin(0.12);
  gPad->SetGrid();
  hP_VT->GetYaxis()->SetRangeUser(0, 0.10);
  hP_VT->SetTitle(";#it{p} (GeV/#it{c});#sigma_{#it{p}}/#it{p}");
  hP_VT->Draw("E1");
  hP_MT->Draw("same E1");
  hP_MS->Draw("same E1");
  TLegend* leg1 = new TLegend(0.5, 0.7, 0.88, 0.88);
  leg1->SetBorderSize(0);
  leg1->SetFillStyle(0);
  leg1->AddEntry(hP_VT, "Vertex telescope only", "lpe");
  leg1->AddEntry(hP_MT, "VT+MS refit", "lpe");
  leg1->AddEntry(hP_MS, "Muon spectrometer only", "lpe");
  leg1->Draw();

  c1.cd(2);
  gPad->SetLeftMargin(0.14);
  gPad->SetRightMargin(0.05);
  gPad->SetTopMargin(0.08);
  gPad->SetBottomMargin(0.12);
  gPad->SetGrid();
  hResX_VT->GetYaxis()->SetRangeUser(0, 50);
  hResX_VT->SetTitle(";#it{p} (GeV/#it{c});#sigma_{d0_{x}} (#mum)");
  hResX_VT->Draw("E1");
  hResX_MT->Draw("same E1");

  c1.cd(3);
  gPad->SetLeftMargin(0.14);
  gPad->SetRightMargin(0.05);
  gPad->SetTopMargin(0.08);
  gPad->SetBottomMargin(0.12);
  gPad->SetGrid();
  hResY_VT->GetYaxis()->SetRangeUser(0, 50);
  hResY_VT->SetTitle(";#it{p} (GeV/#it{c});#sigma_{d0_{y}} (#mum)");
  hResY_VT->Draw("E1");
  hResY_MT->Draw("same E1");

  c1.SaveAs(Form("momentumResolution%s.png", useMC ? "MC" : ""));
  printf("Saved plot: momentumResolution.png\n");

  // ---- Save all histograms to file
  printf("Saving histograms to resolutions%s.root...\n", useMC ? "MC" : "");
  TFile out(Form("resolutions%s.root", useMC ? "MC" : ""), "RECREATE");

  out.mkdir("VerTel");
  out.cd("VerTel");
  hDx_VT->Write();
  hDy_VT->Write();
  hDz_VT->Write();
  hDpx_VT->Write();
  hDpy_VT->Write();
  hDpz_VT->Write();
  hDp_VT->Write();
  hP_VT->Write();
  hPVsP_VT->Write();
  hDxVsP_VT->Write();
  hDyVsP_VT->Write();
  hResX_VT->Write();
  hResY_VT->Write();
  hRej_VT->Write();
  hPID_VT->Write();
  hPIDmuJpsi_VT->Write();

  out.mkdir("MuonSpec");
  out.cd("MuonSpec");
  hDx_MS->Write();
  hDy_MS->Write();
  hDz_MS->Write();
  hDpx_MS->Write();
  hDpy_MS->Write();
  hDpz_MS->Write();
  hDp_MS->Write();
  hP_MS->Write();
  hPVsP_MS->Write();
  hDxVsP_MS->Write();
  hDyVsP_MS->Write();
  hResX_MS->Write();
  hResY_MS->Write();
  hRej_MS->Write();
  hPID_MS->Write();
  hPIDmuJpsi_MS->Write();
  out.mkdir("Matching");
  out.cd("Matching");
  hDx_MT->Write();
  hDy_MT->Write();
  hDz_MT->Write();
  hDpx_MT->Write();
  hDpy_MT->Write();
  hDpz_MT->Write();
  hDp_MT->Write();
  hP_MT->Write();
  hPVsP_MT->Write();
  hDxVsP_MT->Write();
  hDyVsP_MT->Write();
  hResX_MT->Write();
  hResY_MT->Write();
  hRej_MT->Write();
  hPID_MT->Write();
  hPIDmuJpsi_MT->Write();
  out.Close();
  printf("Done! Saved resolutionsMC.root\n");
}