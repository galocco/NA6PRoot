#if !defined(__CINT__) || defined(__MAKECINT__)
#include <TFile.h>
#include <TTree.h>
#include <TH1F.h>
#include <TH3F.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TF1.h>
#include <TLine.h>
#include <TLatex.h>
#include <TStyle.h>
#include <TMath.h>
#include <TDatabasePDG.h>
#include <TGraphAsymmErrors.h>
#include <cmath>
#include <iostream>
#include <vector>
#include <algorithm>
#include <TEfficiency.h>
#include "NA6PLayoutParam.h"
#include "ConfigurableParam.h"
#endif
float effMax = 0.205; // for plotting only
float lambdaBR = 0.639f; // PDG BR for Lambda0 -> p pi- (and c.c.)

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


void PlotLambda0EffPurity(const char* inputFile = "Lambda0-Candidates.root",
                          const char* outPrefix = "Lambda0-EffPurity",
                          const char* layoutFile = "na6pLayout.ini",
                          int nPtBins = 24,
                          float ptMin = 0.f,
                          float ptMax = 3.f,
                          int nYBins = 24,
                          float yMin = 0.f,
                          float yMax = 6.f)
{
  TFile* fin = TFile::Open(inputFile, "READ");
  if (!fin || fin->IsZombie()) {
    std::cerr << "Cannot open input file: " << inputFile << std::endl;
    return;
  }

  TTree* candTree = (TTree*)fin->Get("candidates");
  TTree* genTree = (TTree*)fin->Get("generatedMothers");
  if (!candTree || !genTree) {
    std::cerr << "Missing trees. Need 'candidates' and 'generatedMothers'" << std::endl;
    fin->Close();
    return;
  }

  float mass = 0.f, pt = 0.f, yRec = 0.f, pvzRec = 0.f;
  int isSignal = 0;
  candTree->SetBranchAddress("mass", &mass);
  candTree->SetBranchAddress("pt", &pt);
  candTree->SetBranchAddress("rap", &yRec);
  candTree->SetBranchAddress("pvz", &pvzRec);
  candTree->SetBranchAddress("isSignal", &isSignal);

  float g_pt = 0.f, g_y = 0.f, g_pvz = 0.f;
  genTree->SetBranchAddress("pt", &g_pt);
  genTree->SetBranchAddress("y", &g_y);
  genTree->SetBranchAddress("pvz", &g_pvz);

  na6p::conf::ConfigurableParam::updateFromFile(layoutFile, "", true);
  auto& layoutPar = NA6PLayoutParam::Instance();
  int nTargets = std::max(1, std::min((int)layoutPar.nTargets, (int)NA6PLayoutParam::MaxTargets));
  std::vector<float> zTarget(nTargets, 0.f), halfThick(nTargets, 0.05f);
  for (int it = 0; it < nTargets; ++it) {
    zTarget[it] = layoutPar.posTargetZ[it] + layoutPar.shiftTargets[2];
    halfThick[it] = 0.5f * std::max(1.e-4f, layoutPar.thicknessTarget[it]);
  }

  gStyle->SetOptStat(0);

  double mLambda = 1.115683;
  auto* pLambda = TDatabasePDG::Instance()->GetParticle(3122);
  if (pLambda) {
    mLambda = pLambda->Mass();
  }

  TH1F* hMassSig = new TH1F("hMassSig", ";M_{p#pi} (GeV/#it{c}^{2}); Counts", 240, mLambda - 0.02, mLambda + 0.02);
  TH1F* hMass = new TH1F("hMass", ";M_{p#pi} (GeV/#it{c}^{2}); Counts", 240, mLambda - 0.02, mLambda + 0.02);

  const Long64_t nCand = candTree->GetEntries();
  for (Long64_t i = 0; i < nCand; ++i) {
    candTree->GetEntry(i);
    hMass->Fill(mass);
    if (isSignal == 1) {
      hMassSig->Fill(mass);
    }
  }

  double sigmaMass = hMassSig->GetRMS();
  TF1* fG = new TF1("fG", "gaus", mLambda - 0.012, mLambda + 0.012);
  fG->SetParameters(hMassSig->GetMaximum(), mLambda, 0.002);
  int fitStatus = hMassSig->Fit(fG, "RQ0");
  if (fitStatus == 0 && fG->GetParameter(2) > 1e-5) {
    sigmaMass = fG->GetParameter(2);
  }

  TCanvas* cMass = new TCanvas("cMass", "Lambda0 invariant mass", 700, 500);
  cMass->SetLeftMargin(0.16);
  cMass->SetRightMargin(0.06);
  cMass->SetBottomMargin(0.18);
  cMass->SetTopMargin(0.05);
  auto setAxisStyle1D = [](TH1* hist, float markerSize = 1.25f) {
    hist->SetMarkerSize(markerSize);
    hist->GetXaxis()->SetTitleSize(0.060f);
    hist->GetYaxis()->SetTitleSize(0.060f);
    hist->GetXaxis()->SetLabelSize(0.050f);
    hist->GetYaxis()->SetLabelSize(0.050f);
    hist->GetXaxis()->SetTitleOffset(1.15f);
    hist->GetYaxis()->SetTitleOffset(1.30f);
    hist->GetXaxis()->SetNdivisions(510);
    hist->GetYaxis()->SetNdivisions(510);
  };

  hMass->SetLineColor(kBlue + 1);
  hMass->SetMarkerStyle(20);
  hMass->SetMarkerColor(kBlue + 1);
  setAxisStyle1D(hMass, 1.35f);
  hMassSig->SetLineColor(kRed + 1);
  hMassSig->SetMarkerStyle(20);
  hMassSig->SetMarkerColor(kRed + 1);
  setAxisStyle1D(hMassSig, 1.35f);
  hMass->Draw("hist");
  hMassSig->Draw("hist same");
  cMass->SaveAs(Form("%s-Mass.png", outPrefix));

  const double massLo = mLambda - 3. * sigmaMass;
  const double massHi = mLambda + 3. * sigmaMass;

  TH1F* hPtGen = new TH1F("hPtGen", ";#it{p}_{T} (GeV/#it{c});generated", nPtBins, ptMin, ptMax);
  TH1F* hPtRecSigWin = new TH1F("hPtRecSigWin", ";#it{p}_{T} (GeV/#it{c});reconstructed signal in #pm3#sigma", nPtBins, ptMin, ptMax);
  TH1F* hPtAllWin = new TH1F("hPtAllWin", ";#it{p}_{T} (GeV/#it{c});all candidates in #pm3#sigma", nPtBins, ptMin, ptMax);

  TH1F* hEtaGen = new TH1F("hEtaGen", ";#eta;generated", nYBins, yMin, yMax);
  TH1F* hEtaRecSigWin = new TH1F("hEtaRecSigWin", ";#eta;reconstructed signal in #pm3#sigma", nYBins, yMin, yMax);
  TH1F* hEtaAllWin = new TH1F("hEtaAllWin", ";#eta;all candidates in #pm3#sigma", nYBins, yMin, yMax);

  // 2D histograms: axes are eta (x) and pt (y)
  TH2F* h2Gen       = new TH2F("h2Gen",       ";#eta;#it{p}_{T} (GeV/#it{c})", nYBins, yMin, yMax, nPtBins, ptMin, ptMax);
  TH2F* h2RecSigWin = new TH2F("h2RecSigWin", ";#eta;#it{p}_{T} (GeV/#it{c})", nYBins, yMin, yMax, nPtBins, ptMin, ptMax);
  TH2F* h2AllWin    = new TH2F("h2AllWin",    ";#eta;#it{p}_{T} (GeV/#it{c})", nYBins, yMin, yMax, nPtBins, ptMin, ptMax);

  // 3D histograms: axes are eta (x), pt (y), target index (z)
  TH3F* h3Gen       = new TH3F("h3Gen",       ";#eta;#it{p}_{T} (GeV/#it{c});target index", nYBins, yMin, yMax, nPtBins, ptMin, ptMax, nTargets, -0.5, nTargets - 0.5);
  TH3F* h3RecSigWin = new TH3F("h3RecSigWin", ";#eta;#it{p}_{T} (GeV/#it{c});target index", nYBins, yMin, yMax, nPtBins, ptMin, ptMax, nTargets, -0.5, nTargets - 0.5);

  std::vector<TH1F*> hPtGenT(nTargets), hPtRecSigWinT(nTargets), hPtAllWinT(nTargets);
  std::vector<TH1F*> hEtaGenT(nTargets), hEtaRecSigWinT(nTargets), hEtaAllWinT(nTargets);
  std::vector<TH2F*> h2GenT(nTargets), h2RecSigWinT(nTargets), h2AllWinT(nTargets);
  std::vector<double> nGenT(nTargets, 0.), nSigWinT(nTargets, 0.), nAllWinT(nTargets, 0.);

  for (int it = 0; it < nTargets; ++it) {
    hPtGenT[it]       = new TH1F(Form("hPtGen_t%d", it), ";#it{p}_{T} (GeV/#it{c});generated", nPtBins, ptMin, ptMax);
    hPtRecSigWinT[it] = new TH1F(Form("hPtRecSigWin_t%d", it), ";#it{p}_{T} (GeV/#it{c});reconstructed signal in #pm3#sigma", nPtBins, ptMin, ptMax);
    hPtAllWinT[it]    = new TH1F(Form("hPtAllWin_t%d", it), ";#it{p}_{T} (GeV/#it{c});all candidates in #pm3#sigma", nPtBins, ptMin, ptMax);
    hEtaGenT[it]       = new TH1F(Form("hEtaGen_t%d", it), ";#eta;generated", nYBins, yMin, yMax);
    hEtaRecSigWinT[it] = new TH1F(Form("hEtaRecSigWin_t%d", it), ";#eta;reconstructed signal in #pm3#sigma", nYBins, yMin, yMax);
    hEtaAllWinT[it]    = new TH1F(Form("hEtaAllWin_t%d", it), ";#eta;all candidates in #pm3#sigma", nYBins, yMin, yMax);
    h2GenT[it]       = new TH2F(Form("h2Gen_t%d", it), ";#eta;#it{p}_{T} (GeV/#it{c})", nYBins, yMin, yMax, nPtBins, ptMin, ptMax);
    h2RecSigWinT[it] = new TH2F(Form("h2RecSigWin_t%d", it), ";#eta;#it{p}_{T} (GeV/#it{c})", nYBins, yMin, yMax, nPtBins, ptMin, ptMax);
    h2AllWinT[it]    = new TH2F(Form("h2AllWin_t%d", it), ";#eta;#it{p}_{T} (GeV/#it{c})", nYBins, yMin, yMax, nPtBins, ptMin, ptMax);
  }

  const Long64_t nGen = genTree->GetEntries();
  for (Long64_t i = 0; i < nGen; ++i) {
    genTree->GetEntry(i);
    const double mtGen = std::sqrt(mLambda * mLambda + (double)g_pt * (double)g_pt);
    const double pzGen = mtGen * std::sinh((double)g_y);
    const float etaGen = (std::abs(g_pt) > 1e-9f) ? (float)std::asinh(pzGen / (double)g_pt) : 0.f;
    const int it = GetTargetIndexFromLayout(g_pvz, zTarget, halfThick);
    hPtGen->Fill(g_pt, lambdaBR);
    hEtaGen->Fill(etaGen, lambdaBR);
    h2Gen->Fill(etaGen, g_pt, lambdaBR);
    h3Gen->Fill(etaGen, g_pt, it, lambdaBR);
    hPtGenT[it]->Fill(g_pt, lambdaBR);
    hEtaGenT[it]->Fill(etaGen, lambdaBR);
    h2GenT[it]->Fill(etaGen, g_pt, lambdaBR);
    nGenT[it] += lambdaBR;
  }

  double nSigWin = 0.;
  double nAllWin = 0.;
  for (Long64_t i = 0; i < nCand; ++i) {
    candTree->GetEntry(i);
    if (mass < massLo || mass > massHi) {
      continue;
    }

    const double mtRec = std::sqrt((double)mass * (double)mass + (double)pt * (double)pt);
    const double pzRec = mtRec * std::sinh((double)yRec);
    const float etaRec = (std::abs(pt) > 1e-9f) ? (float)std::asinh(pzRec / (double)pt) : 0.f;
    const int it = GetTargetIndexFromLayout(pvzRec, zTarget, halfThick);

    hPtAllWin->Fill(pt);
    hEtaAllWin->Fill(etaRec);
    h2AllWin->Fill(etaRec, pt);
    hPtAllWinT[it]->Fill(pt);
    hEtaAllWinT[it]->Fill(etaRec);
    h2AllWinT[it]->Fill(etaRec, pt);
    nAllWin += 1.;
    nAllWinT[it] += 1.;

    if (isSignal == 1) {
      hPtRecSigWin->Fill(pt);
      hEtaRecSigWin->Fill(etaRec);
      h2RecSigWin->Fill(etaRec, pt);
      h3RecSigWin->Fill(etaRec, pt, it);
      hPtRecSigWinT[it]->Fill(pt);
      hEtaRecSigWinT[it]->Fill(etaRec);
      h2RecSigWinT[it]->Fill(etaRec, pt);
      nSigWin += 1.;
      nSigWinT[it] += 1.;
    }
  }

  TCanvas* cGen = new TCanvas("cGen", "Lambda0 generated distributions", 1200, 500);
  cGen->Divide(2, 1);
  cGen->cd(1);
  gPad->SetLeftMargin(0.14);
  gPad->SetBottomMargin(0.14);
  hPtGen->SetLineColor(kBlue + 1);
  hPtGen->SetMarkerColor(kBlue + 1);
  hPtGen->SetMarkerStyle(20);
  setAxisStyle1D(hPtGen, 1.2f);
  hPtGen->SetTitle("Generated #Lambda^{0};#it{p}_{T} (GeV/#it{c});Counts");
  hPtGen->Draw("hist");

  cGen->cd(2);
  gPad->SetLeftMargin(0.14);
  gPad->SetBottomMargin(0.14);
  hEtaGen->SetLineColor(kRed + 1);
  hEtaGen->SetMarkerColor(kRed + 1);
  hEtaGen->SetMarkerStyle(21);
  setAxisStyle1D(hEtaGen, 1.2f);
  hEtaGen->SetTitle("Generated #Lambda^{0};#eta;Counts");
  hEtaGen->Draw("hist");

  cGen->SaveAs(Form("%s-Gen.png", outPrefix));

  TEfficiency* hEffVsPt = new TEfficiency(*hPtRecSigWin, *hPtGen);
  TEfficiency* hEffVsEta = new TEfficiency(*hEtaRecSigWin, *hEtaGen);

  TEfficiency* hPurVsPt = new TEfficiency(*hPtRecSigWin, *hPtAllWin);
  TEfficiency* hPurVsEta = new TEfficiency(*hEtaRecSigWin, *hEtaAllWin);

  // 2D efficiency: use TH2F ratio with binomial errors
  // (TEfficiency 2D fails when reco kinematics differ from gen kinematics)
  TH2F* hEff2D = (TH2F*)h2RecSigWin->Clone("hEff2D");
  hEff2D->SetTitle("Efficiency x Acceptance;#eta;#it{p}_{T} (GeV/#it{c})");
  hEff2D->Divide(h2RecSigWin, h2Gen, 1., 1., "B");

  TH3F* hEff3D = (TH3F*)h3RecSigWin->Clone("hEff3D");
  hEff3D->SetTitle("Efficiency x Acceptance;#eta;#it{p}_{T} (GeV/#it{c});target index");
  hEff3D->Divide(h3RecSigWin, h3Gen, 1., 1., "B");



  TH2F* hPur2D = (TH2F*)h2RecSigWin->Clone("hPur2D");
  hPur2D->SetTitle("Purity ;#eta;#it{p}_{T} (GeV/#it{c})");
  hPur2D->Divide(h2RecSigWin, h2AllWin, 1., 1., "B");

  const double effGlobal = (nGen > 0.) ? (nSigWin / nGen / lambdaBR) : 0.;
  const double purGlobal = (nAllWin > 0.) ? (nSigWin / nAllWin) : 0.;

  TCanvas* c = new TCanvas("cEffPur", "Lambda0 efficiency and purity", 1400, 1000);
  c->Divide(2, 2);

  auto drawEff = [](TEfficiency* eff, Color_t color, Style_t markerStyle) {
    eff->SetLineColor(color);
    eff->SetMarkerColor(color);
    eff->SetMarkerStyle(markerStyle);
    eff->SetMarkerSize(1.25);
    eff->Draw("AP");
    gPad->Update();
    if (auto* gr = eff->GetPaintedGraph()) {
      gr->SetMinimum(0.0);
      gr->SetMaximum(effMax);
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

  c->cd(1);
  drawEff(hEffVsPt, kRed + 1, 20);

  c->cd(2);
  drawEff(hEffVsEta, kRed + 1, 20);

  c->cd(3);
  drawEff(hPurVsPt, kBlue + 1, 21);

  c->cd(4);
  drawEff(hPurVsEta, kBlue + 1, 21);

  TLatex lab;
  lab.SetNDC();
  lab.SetTextSize(0.032);
  c->cd(1);
  lab.DrawLatex(0.16, 0.87, Form("m_{#Lambda}=%.6f GeV/#it{c}^{2}, #sigma=%.4f GeV/#it{c}^{2}", mLambda, sigmaMass));
  lab.DrawLatex(0.16, 0.81, Form("Global eff = %.4f", effGlobal));
  c->cd(3);
  lab.DrawLatex(0.16, 0.87, Form("Mass window: [%.6f, %.6f] GeV/#it{c}^{2}", massLo, massHi));
  lab.DrawLatex(0.16, 0.81, Form("Global purity = %.4f", purGlobal));

  c->SaveAs(Form("%s.png", outPrefix));

  // ── 2D efficiency / purity canvas ─────────────────────────────────────────
  TCanvas* c2d = new TCanvas("cEffPur2D", "Lambda0 2D efficiency and purity (#it{p}_{T} vs #eta)", 1200, 500);
  c2d->Divide(2, 1);
  for (int ipad = 1; ipad <= 2; ++ipad) {
    c2d->cd(ipad);
    gPad->SetLeftMargin(0.16);
    gPad->SetRightMargin(0.12);
    gPad->SetBottomMargin(0.18);
    gPad->SetTopMargin(0.07);
  }

  auto draw2DRatio = [](TPad* pad, TH2F* h2) {
    pad->cd();
    gStyle->SetPaintTextFormat(".2f");
    gStyle->SetPalette(kBird);
    h2->SetMinimum(0.);
    h2->SetMaximum(1.);
    h2->GetXaxis()->SetTitleSize(0.060f);
    h2->GetYaxis()->SetTitleSize(0.060f);
    h2->GetXaxis()->SetLabelSize(0.050f);
    h2->GetYaxis()->SetLabelSize(0.050f);
    h2->GetXaxis()->SetTitleOffset(1.20f);
    h2->GetYaxis()->SetTitleOffset(1.20f);
    h2->GetXaxis()->SetNdivisions(510);
    h2->GetYaxis()->SetNdivisions(510);
    h2->Draw("colz text");
    gPad->Update();
  };

  draw2DRatio((TPad*)c2d->cd(1), hEff2D);
  draw2DRatio((TPad*)c2d->cd(2), hPur2D);

  c2d->SaveAs(Form("%s-2D.png", outPrefix));

  std::vector<TH1F*> hEffVsPtT(nTargets), hEffVsEtaT(nTargets), hPurVsPtT(nTargets), hPurVsEtaT(nTargets);
  std::vector<TH2F*> hEff2DT(nTargets), hPur2DT(nTargets);
  TCanvas* cEffOverlay = nullptr;
  TCanvas* cPurOverlay = nullptr;
  for (int it = 0; it < nTargets; ++it) {
    hEffVsPtT[it] = (TH1F*)hPtRecSigWinT[it]->Clone(Form("hEffVsPt_t%d", it));
    hEffVsPtT[it]->SetTitle(Form("Target %d (z=%.3f cm);#it{p}_{T} (GeV/#it{c});Efficiency x Acceptance", it, zTarget[it]));
    hEffVsPtT[it]->Divide(hPtRecSigWinT[it], hPtGenT[it], 1., 1., "B");

    hEffVsEtaT[it] = (TH1F*)hEtaRecSigWinT[it]->Clone(Form("hEffVsEta_t%d", it));
    hEffVsEtaT[it]->SetTitle(Form("Target %d (z=%.3f cm);#eta;Efficiency x Acceptance", it, zTarget[it]));
    hEffVsEtaT[it]->Divide(hEtaRecSigWinT[it], hEtaGenT[it], 1., 1., "B");

    hPurVsPtT[it] = (TH1F*)hPtRecSigWinT[it]->Clone(Form("hPurVsPt_t%d", it));
    hPurVsPtT[it]->SetTitle(Form("Target %d (z=%.3f cm);#it{p}_{T} (GeV/#it{c});Purity", it, zTarget[it]));
    hPurVsPtT[it]->Divide(hPtRecSigWinT[it], hPtAllWinT[it], 1., 1., "B");

    hPurVsEtaT[it] = (TH1F*)hEtaRecSigWinT[it]->Clone(Form("hPurVsEta_t%d", it));
    hPurVsEtaT[it]->SetTitle(Form("Target %d (z=%.3f cm);#eta;Purity", it, zTarget[it]));
    hPurVsEtaT[it]->Divide(hEtaRecSigWinT[it], hEtaAllWinT[it], 1., 1., "B");

    hEff2DT[it] = (TH2F*)h2RecSigWinT[it]->Clone(Form("hEff2D_t%d", it));
    hEff2DT[it]->SetTitle(Form("Target %d (z=%.3f cm) Efficiency x Acceptance;#eta;#it{p}_{T} (GeV/#it{c})", it, zTarget[it]));
    hEff2DT[it]->Divide(h2RecSigWinT[it], h2GenT[it], 1., 1., "B");

    hPur2DT[it] = (TH2F*)h2RecSigWinT[it]->Clone(Form("hPur2D_t%d", it));
    hPur2DT[it]->SetTitle(Form("Target %d (z=%.3f cm) Purity;#eta;#it{p}_{T} (GeV/#it{c})", it, zTarget[it]));
    hPur2DT[it]->Divide(h2RecSigWinT[it], h2AllWinT[it], 1., 1., "B");

    TCanvas* cT = new TCanvas(Form("cEffPur_t%d", it), Form("Target %d efficiency and purity", it), 1400, 1000);
    cT->Divide(2, 2);
    cT->cd(1); hEffVsPtT[it]->SetMinimum(0.); hEffVsPtT[it]->SetMaximum(effMax); hEffVsPtT[it]->SetLineColor(kRed+1); hEffVsPtT[it]->SetMarkerColor(kRed+1); hEffVsPtT[it]->SetMarkerStyle(20); setAxisStyle1D(hEffVsPtT[it], 1.25f); hEffVsPtT[it]->Draw("E1");
    cT->cd(2); hEffVsEtaT[it]->SetMinimum(0.); hEffVsEtaT[it]->SetMaximum(effMax); hEffVsEtaT[it]->SetLineColor(kRed+1); hEffVsEtaT[it]->SetMarkerColor(kRed+1); hEffVsEtaT[it]->SetMarkerStyle(20); setAxisStyle1D(hEffVsEtaT[it], 1.25f); hEffVsEtaT[it]->Draw("E1");
    cT->cd(3); hPurVsPtT[it]->SetMinimum(0.); hPurVsPtT[it]->SetMaximum(1.05); hPurVsPtT[it]->SetLineColor(kBlue+1); hPurVsPtT[it]->SetMarkerColor(kBlue+1); hPurVsPtT[it]->SetMarkerStyle(21); setAxisStyle1D(hPurVsPtT[it], 1.25f); hPurVsPtT[it]->Draw("E1");
    cT->cd(4); hPurVsEtaT[it]->SetMinimum(0.); hPurVsEtaT[it]->SetMaximum(1.05); hPurVsEtaT[it]->SetLineColor(kBlue+1); hPurVsEtaT[it]->SetMarkerColor(kBlue+1); hPurVsEtaT[it]->SetMarkerStyle(21); setAxisStyle1D(hPurVsEtaT[it], 1.25f); hPurVsEtaT[it]->Draw("E1");

    TLatex tPerTarget;
    tPerTarget.SetNDC();
    tPerTarget.SetTextSize(0.030);
    const double effTarget = (nGenT[it] > 0.) ? (nSigWinT[it] / nGenT[it]) : 0.;
    cT->cd(1);
    tPerTarget.DrawLatex(0.14, 0.88, Form("Target %d  z=%.3f cm", it, zTarget[it]));
    tPerTarget.DrawLatex(0.14, 0.82, Form("Eff(target) = %.4f", effTarget));
    tPerTarget.DrawLatex(0.14, 0.76, Form("Eff(global) = %.4f", effGlobal));
    cT->cd(2);
    tPerTarget.DrawLatex(0.14, 0.88, Form("Eff(global) = %.4f", effGlobal));

    TCanvas* c2dT = new TCanvas(Form("cEffPur2D_t%d", it), Form("Target %d 2D efficiency and purity", it), 1200, 500);
    c2dT->Divide(2, 1);
    for (int ipad = 1; ipad <= 2; ++ipad) {
      c2dT->cd(ipad);
      gPad->SetLeftMargin(0.16);
      gPad->SetRightMargin(0.12);
      gPad->SetBottomMargin(0.18);
      gPad->SetTopMargin(0.07);
    }
    draw2DRatio((TPad*)c2dT->cd(1), hEff2DT[it]);
    draw2DRatio((TPad*)c2dT->cd(2), hPur2DT[it]);

    c2dT->cd(1);
    tPerTarget.DrawLatex(0.12, 0.93, Form("Eff(global) = %.4f", effGlobal));

    cT->SaveAs(Form("%s-Target%d.png", outPrefix, it));
    c2dT->SaveAs(Form("%s-Target%d-2D.png", outPrefix, it));
  }

  cEffOverlay = new TCanvas("cEffOverlay", "Lambda0 efficiencies overlap", 1400, 600);
  cEffOverlay->Divide(2, 1);
  for (int ipad = 1; ipad <= 2; ++ipad) {
    cEffOverlay->cd(ipad);
    gPad->SetLeftMargin(0.16);
    gPad->SetRightMargin(0.04);
    gPad->SetBottomMargin(0.16);
    gPad->SetTopMargin(0.05);
  }

  std::vector<Color_t> colors = {kBlack, kRed + 1, kBlue + 1, kGreen + 2, kMagenta + 1, kOrange + 7, kCyan + 1, kViolet + 1, kTeal + 2, kPink + 7};
  auto drawEffOverlay = [](TEfficiency* eff, Color_t color, Style_t markerStyle, float effMax = 0.15f) {
    eff->SetLineColor(color);
    eff->SetMarkerColor(color);
    eff->SetMarkerStyle(markerStyle);
    eff->SetMarkerSize(1.25);
    eff->Draw("AP");
    gPad->Update();
    if (auto* gr = eff->GetPaintedGraph()) {
      gr->SetMinimum(0.0);
      gr->SetMaximum(effMax);
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

  cEffOverlay->cd(1);
  hEffVsEta->SetTitle(";#eta;Efficiency x Acceptance");
  drawEffOverlay(hEffVsEta, colors[0], 20, effMax);

  TLegend* legEta = new TLegend(0.15, 0.60, 0.4, 0.88);
  legEta->SetBorderSize(0);
  legEta->SetFillStyle(0);
  legEta->AddEntry(hEffVsEta, Form("Global eff = %.2f%%", 100. * effGlobal), "lp");

  cEffOverlay->cd(2);
  hEffVsPt->SetTitle(";#it{p}_{T} (GeV/#it{c});Efficiency x Acceptance");
  drawEffOverlay(hEffVsPt, colors[0], 20, effMax);

  TLegend* legPt = new TLegend(0.15, 0.60, 0.4, 0.88);
  legPt->SetBorderSize(0);
  legPt->SetFillStyle(0);
  legPt->AddEntry(hEffVsPt, Form("Global eff = %.2f%%", 100. * effGlobal), "lp");

  for (int it = 0; it < nTargets; ++it) {
    Color_t ccol = colors[(it + 1) % colors.size()];
    Style_t mstyle = 21 + (it % 9);
    const double effTarget = (nGenT[it] > 0.) ? (nSigWinT[it] / nGenT[it]) : 0.;

    hEffVsEtaT[it]->SetLineColor(ccol);
    hEffVsEtaT[it]->SetMarkerColor(ccol);
    hEffVsEtaT[it]->SetMarkerStyle(mstyle);
    hEffVsEtaT[it]->SetMarkerSize(1.25f);
    cEffOverlay->cd(1);
    hEffVsEtaT[it]->SetMinimum(0.);
    hEffVsEtaT[it]->SetMaximum(effMax);
    hEffVsEtaT[it]->Draw("E1 SAME");
    legEta->AddEntry(hEffVsEtaT[it], Form("Target %d eff = %.2f%%", it, 100. * effTarget), "lp");

    hEffVsPtT[it]->SetLineColor(ccol);
    hEffVsPtT[it]->SetMarkerColor(ccol);
    hEffVsPtT[it]->SetMarkerStyle(mstyle);
    hEffVsPtT[it]->SetMarkerSize(1.25f);
    cEffOverlay->cd(2);
    hEffVsPtT[it]->SetMinimum(0.);
    hEffVsPtT[it]->SetMaximum(effMax);
    hEffVsPtT[it]->Draw("E1 SAME");
    legPt->AddEntry(hEffVsPtT[it], Form("Target %d eff = %.2f%%", it, 100. * effTarget), "lp");
  }

  cEffOverlay->cd(1);
  legEta->Draw();
  cEffOverlay->cd(2);
  legPt->Draw();

  cEffOverlay->SaveAs(Form("%s-EffOverlay.png", outPrefix));

  cPurOverlay = new TCanvas("cPurOverlay", "Lambda0 purities overlap", 1400, 600);
  cPurOverlay->Divide(2, 1);
  for (int ipad = 1; ipad <= 2; ++ipad) {
    cPurOverlay->cd(ipad);
    gPad->SetLeftMargin(0.16);
    gPad->SetRightMargin(0.04);
    gPad->SetBottomMargin(0.16);
    gPad->SetTopMargin(0.05);
  }

  cPurOverlay->cd(1);
  hPurVsEta->SetTitle(";#eta;Purity");
  drawEffOverlay(hPurVsEta, colors[0], 20, 1.05);

  TLegend* legPurEta = new TLegend(0.15, 0.2, 0.4, 0.48);
  legPurEta->SetBorderSize(0);
  legPurEta->SetFillStyle(0);
  legPurEta->AddEntry(hPurVsEta, Form("Global pur = %.2f%%", 100. * purGlobal), "lp");

  cPurOverlay->cd(2);
  hPurVsPt->SetTitle(";#it{p}_{T} (GeV/#it{c});Purity");
  drawEffOverlay(hPurVsPt, colors[0], 20, 1.05);

  TLegend* legPurPt = new TLegend(0.15, 0.2, 0.4, 0.48);
  legPurPt->SetBorderSize(0);
  legPurPt->SetFillStyle(0);
  legPurPt->AddEntry(hPurVsPt, Form("Global pur = %.2f%%", 100. * purGlobal), "lp");

  for (int it = 0; it < nTargets; ++it) {
    Color_t ccol = colors[(it + 1) % colors.size()];
    Style_t mstyle = 21 + (it % 9);
    const double purTarget = (nAllWinT[it] > 0.) ? (nSigWinT[it] / nAllWinT[it]) : 0.;

    hPurVsEtaT[it]->SetLineColor(ccol);
    hPurVsEtaT[it]->SetMarkerColor(ccol);
    hPurVsEtaT[it]->SetMarkerStyle(mstyle);
    hPurVsEtaT[it]->SetMarkerSize(1.25f);
    cPurOverlay->cd(1);
    hPurVsEtaT[it]->SetMinimum(0.);
    hPurVsEtaT[it]->SetMaximum(1.05);
    hPurVsEtaT[it]->Draw("E1 SAME");
    legPurEta->AddEntry(hPurVsEtaT[it], Form("Target %d pur = %.2f%%", it, 100. * purTarget), "lp");

    hPurVsPtT[it]->SetLineColor(ccol);
    hPurVsPtT[it]->SetMarkerColor(ccol);
    hPurVsPtT[it]->SetMarkerStyle(mstyle);
    hPurVsPtT[it]->SetMarkerSize(1.25f);
    cPurOverlay->cd(2);
    hPurVsPtT[it]->SetMinimum(0.);
    hPurVsPtT[it]->SetMaximum(1.05);
    hPurVsPtT[it]->Draw("E1 SAME");
    legPurPt->AddEntry(hPurVsPtT[it], Form("Target %d pur = %.2f%%", it, 100. * purTarget), "lp");
  }

  cPurOverlay->cd(1);
  legPurEta->Draw();
  cPurOverlay->cd(2);
  legPurPt->Draw();

  cPurOverlay->SaveAs(Form("%s-PurOverlay.png", outPrefix));

  TFile* fout = new TFile(Form("%s.root", outPrefix), "RECREATE");
  hMassSig->Write();
  hMass->Write();
  fG->Write();
  hPtGen->Write();
  hPtRecSigWin->Write();
  hPtAllWin->Write();
  hEtaGen->Write();
  hEtaRecSigWin->Write();
  hEtaAllWin->Write();
  hEffVsPt->Write();
  hEffVsEta->Write();
  hPurVsPt->Write();
  hPurVsEta->Write();
  h2Gen->Write();
  h2RecSigWin->Write();
  h2AllWin->Write();
  h3Gen->Write();
  h3RecSigWin->Write();
  hEff2D->Write();
  hEff3D->Write();
  hPur2D->Write();
  cMass->Write();
  cGen->Write();
  c->Write();
  c2d->Write();
  if (cEffOverlay) {
    cEffOverlay->Write();
  }
  if (cPurOverlay) {
    cPurOverlay->Write();
  }
  for (int it = 0; it < nTargets; ++it) {
    hPtGenT[it]->Write();
    hPtRecSigWinT[it]->Write();
    hPtAllWinT[it]->Write();
    hEtaGenT[it]->Write();
    hEtaRecSigWinT[it]->Write();
    hEtaAllWinT[it]->Write();
    h2GenT[it]->Write();
    h2RecSigWinT[it]->Write();
    h2AllWinT[it]->Write();
    hEffVsPtT[it]->Write();
    hEffVsEtaT[it]->Write();
    hPurVsPtT[it]->Write();
    hPurVsEtaT[it]->Write();
    hEff2DT[it]->Write();
    hPur2DT[it]->Write();
  }
  fout->Close();

  std::cout << "\n=== Lambda0 efficiency/purity summary ===\n";
  std::cout << "Input file: " << inputFile << "\n";
  std::cout << "mLambda (PDG): " << mLambda << " GeV/#it{c}^2\n";
  std::cout << "sigmaLambda (fit): " << sigmaMass << " GeV/#it{c}^2\n";
  std::cout << "Mass window: [" << massLo << ", " << massHi << "] GeV/#it{c}^2\n";
  std::cout << "N generated: " << nGen << "\n";
  std::cout << "N signal in window: " << nSigWin << "\n";
  std::cout << "N all candidates in window: " << nAllWin << "\n";
  std::cout << "Global efficiency: " << effGlobal << "\n";
  std::cout << "Global purity: " << purGlobal << "\n";
  std::cout << "Targets from layout: " << nTargets << "\n";
  std::cout << "Saved: " << outPrefix << ".png and " << outPrefix << ".root\n\n";

  fin->Close();
}
