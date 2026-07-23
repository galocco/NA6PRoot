#if !defined(__CINT__) || defined(__MAKECINT__)
#include <TTree.h>
#include <TFile.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TF1.h>
#include <TParticle.h>
#include <TCanvas.h>
#include <TPaveStats.h>
#include <TLegend.h>
#include "NA6PTrack.h"
#include "MagneticField.h"
#include "Propagator.h"
#include "NA6PVerTelHit.h"
#include "NA6PMCComposedLabel.h"
#endif

void fillMeanAndRms(TH2F* h2d, TH1F* hMean, TH1F* hRms, TH1F* hSig)
{
  if (!h2d || !hMean || !hRms || !hSig) {
    printf("One of pointers is not set: h2d=%p hMean=%p hRms=%p hSig=%p\n", h2d, hMean, hRms, hSig);
    return;
  }
  int nMomBins = h2d->GetXaxis()->GetNbins();
  hMean->SetStats(0);
  hRms->SetStats(0);
  hSig->SetStats(0);
  for (int jp = 1; jp <= nMomBins; jp++) {
    TH1D* htmp1 = h2d->ProjectionY("htmp1", jp, jp);
    double mean = htmp1->GetMean();
    double emean = htmp1->GetMeanError();
    double rms = htmp1->GetRMS();
    double erms = htmp1->GetRMSError();
    hMean->SetBinContent(jp, mean);
    hMean->SetBinError(jp, emean);
    hRms->SetBinContent(jp, rms);
    hRms->SetBinError(jp, erms);
    double gw = rms;
    double egw = erms;
    if (htmp1->GetEntries() > 20. && htmp1->Fit("gaus", "Q", "", -3. * rms, 3. * rms)) {
      TF1* fg = (TF1*)htmp1->GetListOfFunctions()->FindObject("gaus");
      if (!fg)
        continue;
      gw = fg->GetParameter(2);
      egw = fg->GetParError(2);
    }
    hSig->SetBinContent(jp, gw);
    hSig->SetBinError(jp, egw);
    delete htmp1;
  }
}

void superposHistos(TH1* h1, TH1* h2, TH1* h3, int col1 = 1, int col2 = kGreen + 1, int col3 = 2)
{
  h1->SetLineColor(col1);
  h1->SetLineWidth(2);
  if (h2->GetMaximum() > h1->GetMaximum())
    h1->SetMaximum(1.02 * h2->GetMaximum());
  if (h3->GetMaximum() > h1->GetMaximum())
    h1->SetMaximum(1.02 * h3->GetMaximum());
  h1->SetMinimum(0.5);
  h1->Draw();
  h2->SetLineColor(col2);
  h2->SetLineWidth(2);
  h2->Draw("sames");
  h3->SetLineColor(col3);
  h3->SetLineWidth(2);
  h3->Draw("sames");
  gPad->Update();
  TPaveStats* st1 = (TPaveStats*)h1->GetListOfFunctions()->FindObject("stats");
  if (st1) {
    st1->SetY1NDC(0.72);
    st1->SetY2NDC(0.92);
    st1->SetTextColor(h1->GetLineColor());
  }
  TPaveStats* st2 = (TPaveStats*)h2->GetListOfFunctions()->FindObject("stats");
  if (st2) {
    st2->SetY1NDC(0.51);
    st2->SetY2NDC(0.71);
    st2->SetTextColor(h2->GetLineColor());
  }
  TPaveStats* st3 = (TPaveStats*)h3->GetListOfFunctions()->FindObject("stats");
  if (st3) {
    st3->SetY1NDC(0.30);
    st3->SetY2NDC(0.50);
    st3->SetTextColor(h3->GetLineColor());
  }
  gPad->Modified();
}

void plotTracks(const char* dirSimu = ".")
{
  auto magField = new MagneticField();
  magField->loadField();
  magField->setAsGlobalField();

  TFile* ft = new TFile(Form("%s/TracksVerTel.root", dirSimu));
  if (!ft)
    return;
  TTree* trTree = (TTree*)ft->Get("tracksVerTel");
  std::vector<NA6PTrack>* trArr = nullptr;
  std::vector<NA6PMCComposedLabel>* labArr = nullptr;
  trTree->SetBranchAddress("VerTel", &trArr);
  trTree->SetBranchAddress("VerTelMCTruth", &labArr);

  TFile* fk = new TFile(Form("%s/MCKine.root", dirSimu));
  TTree* mcTree = (TTree*)fk->Get("mckine");
  std::vector<TParticle>* mcArr = nullptr;
  mcTree->SetBranchAddress("tracks", &mcArr);

  TFile* fh = new TFile(Form("%s/HitsVerTel.root", dirSimu));
  TTree* th = (TTree*)fh->Get("hitsVerTel");
  std::vector<NA6PVerTelHit> vtHits, *vtHitsPtr = &vtHits;
  th->SetBranchAddress("VerTel", &vtHitsPtr);

  int nMomBins = 20;
  double maxP = 10.;
  TH1F* hMomGen = new TH1F("hMomGen", ";p_{gen} (GeV/c);counts", nMomBins, 0., maxP);
  TH1F* hEtaGen = new TH1F("hEtaGen", ";#eta_{gen};counts", 20, 1., 5.);
  TH1F* hMomTrackable = new TH1F("hMomTrackable", ";p_{gen} (GeV/c);counts", nMomBins, 0., maxP);
  TH1F* hEtaTrackable = new TH1F("hEtaTrackable", ";#eta_{gen};counts", 20, 1., 5.);
  TH1F* hMomTracked = new TH1F("hMomTracked", ";p_{gen} (GeV/c);counts", nMomBins, 0., maxP);
  TH1F* hEtaTracked = new TH1F("hEtaTracked", ";#eta_{gen};counts", 20, 1., 5.);
  TH1F* hMomAllReco = new TH1F("hMomAllReco", ";p (GeV/c); counts", nMomBins, 0., maxP);
  TH1F* hMomGoodReco = new TH1F("hMomGoodReco", ";p (GeV/c); counts", nMomBins, 0., maxP);
  TH1F* hMomFakeReco = new TH1F("hMomFakeReco", ";p (GeV/c); counts", nMomBins, 0., maxP);
  TH1F* hEtaAllReco = new TH1F("hEtaAllReco", ";#eta; counts", 20, 1., 5.);
  TH1F* hEtaGoodReco = new TH1F("hEtaGoodReco", ";#eta; counts", 20, 1., 5.);
  TH1F* hEtaFakeReco = new TH1F("hEtaFakeReco", ";#eta; counts", 20, 1., 5.);
  TH1F* hNclu = new TH1F("hNclu", ";n_{ITSclus};counts", 7, -0.5, 6.5);
  double impMax = 500.;  // in microns
  double deltaMax = 0.5; // in GeV/c
  TH2F* hImpParXVsP = new TH2F("hImpParXVsP", ";p (GeV/c);Track Imp. Par. X (#mum)};counts", nMomBins, 0., maxP, 100, -impMax, impMax);
  TH2F* hImpParYVsP = new TH2F("hImpParYVsP", ";p (GeV/c);Track Imp. Par. Y (#mum)};counts", nMomBins, 0., maxP, 100, -impMax, impMax);
  TH2F* hDeltaPxVsP = new TH2F("hDeltaPxVsP", ";p (GeV/c);p_{x}^{rec}-p_{x}^{gen} (GeV/c);counts", nMomBins, 0., maxP, 500, -deltaMax, deltaMax);
  TH2F* hDeltaPyVsP = new TH2F("hDeltaPyVsP", ";p (GeV/c);p_{y}^{rec}-p_{y}^{gen} (GeV/c);counts", nMomBins, 0., maxP, 500, -deltaMax, deltaMax);
  TH2F* hDeltaPzVsP = new TH2F("hDeltaPzVsP", ";p (GeV/c);p_{z}^{rec}-p_{z}^{gen} (GeV/c);counts", nMomBins, 0., maxP, 500, -deltaMax, deltaMax);
  TH2F* hRelDeltaPxVsP = new TH2F("hRelDeltaPxVsP", ";p (GeV/c);(p_{x}^{rec}-p_{x}^{gen})/p_{x}^{gen};counts", nMomBins, 0., maxP, 100, -deltaMax, deltaMax);
  TH2F* hRelDeltaPyVsP = new TH2F("hRelDeltaPyVsP", ";p (GeV/c);(p_{y}^{rec}-p_{y}^{gen})/p_{y}^{gen};counts", nMomBins, 0., maxP, 100, -deltaMax, deltaMax);
  TH2F* hRelDeltaPzVsP = new TH2F("hRelDeltaPzVsP", ";p (GeV/c);(p_{z}^{rec}-p_{z}^{gen})/p_{z}^{gen};counts", nMomBins, 0., maxP, 100, -deltaMax, deltaMax);

  int nEv = trTree->GetEntries();
  printf("Number of events = %d\n", nEv);
  for (int jEv = 0; jEv < nEv; jEv++) {
    mcTree->GetEvent(jEv);
    trTree->GetEvent(jEv);
    th->GetEvent(jEv);
    int nPart = mcArr->size();
    int nTracks = trArr->size();
    int nHits = vtHits.size();
    printf("Event %d particles = %d tracks = %d\n", jEv, nPart, nTracks);
    double xvert = 0;
    double yvert = 0;
    double zvert = 0;

    for (int jp = 0; jp < nPart; jp++) {
      auto curPart = mcArr->at(jp);
      if (curPart.IsPrimary()) {
        xvert = curPart.Vx();
        yvert = curPart.Vy();
        zvert = curPart.Vz();
      }
      int jDau = curPart.GetFirstDaughter();
      double zDecay = 999.;
      double zOrig = curPart.Vz();
      if (jDau >= 0) {
        TParticle dauPart = mcArr->at(jDau);
        zDecay = dauPart.Vz();
      }
      double pxPart = curPart.Px();
      double pyPart = curPart.Py();
      double pzPart = curPart.Pz();
      double momPart = curPart.P();
      double phiPart = curPart.Phi();
      double thetaPart = std::acos(pzPart / momPart);
      double etaPart = -std::log(std::tan(thetaPart / 2.));
      if (zOrig < 7 && zDecay > 40 && etaPart > 1 && etaPart < 4) {
        hMomGen->Fill(momPart);
        hEtaGen->Fill(etaPart);
      }
      int maskHits = 0;
      for (int jHit = 0; jHit < nHits; ++jHit) {
        const auto& hit = vtHits.at(jHit);
        int idPart = hit.getTrackID();
        if (idPart == jp) {
          int nLay = hit.getDetectorID() / 4;
          maskHits += (1 << nLay);
        }
      }
      if (maskHits == 31) {
        hMomTrackable->Fill(momPart);
        hEtaTrackable->Fill(etaPart);
      }
    }
    for (int jTr = 0; jTr < nTracks; ++jTr) {
      NA6PTrack tr = trArr->at(jTr);
      Propagator::PropOpt propOpt;
      propOpt.matCorr = Propagator::MatCorrType::USEMatCorrNONE;
      Propagator::Instance()->propagateToZ(tr, zvert, propOpt);
      NA6PMCComposedLabel mcCompLabel = labArr->at(jTr);
      int nClusters = tr.getNHits();
      hNclu->Fill(nClusters);
      if (nClusters < 5)
        continue;
      auto pxyzReco = tr.getPXYZ<double>();
      double pxReco = pxyzReco[0];
      double pyReco = pxyzReco[1];
      double pzReco = pxyzReco[2];
      double ptReco = std::sqrt(pxyzReco[0] * pxyzReco[0] + pxyzReco[1] * pxyzReco[1]);
      double momReco = tr.getP();
      double thetaReco = std::acos(pxyzReco[2] / momReco);
      double etaReco = -std::log(std::tan(thetaReco / 2.));
      double impparX = tr.getX() - xvert;
      double impparY = tr.getY() - yvert;
      hImpParXVsP->Fill(momReco, impparX * 1e4);
      hImpParYVsP->Fill(momReco, impparY * 1e4);
      hMomAllReco->Fill(momReco);
      hEtaAllReco->Fill(etaReco);
      int mcLabel = mcCompLabel.getTrackID();
      if (!mcCompLabel.isFake()) {
        hMomGoodReco->Fill(momReco);
        hEtaGoodReco->Fill(etaReco);
      } else {
        hMomFakeReco->Fill(momReco);
        hEtaFakeReco->Fill(etaReco);
      }

      TParticle part = mcArr->at(mcLabel);
      double pxPart = part.Px();
      double pyPart = part.Py();
      double pzPart = part.Pz();
      double momPart = part.P();
      double phiPart = part.Phi();
      double thetaPart = std::acos(pzPart / momPart);
      double etaPart = -std::log(std::tan(thetaPart / 2.));
      hMomTracked->Fill(momPart);
      hEtaTracked->Fill(etaPart);
      hDeltaPxVsP->Fill(momReco, pxReco - pxPart);
      hDeltaPyVsP->Fill(momReco, pyReco - pyPart);
      hDeltaPzVsP->Fill(momReco, pzReco - pzPart);
      hRelDeltaPxVsP->Fill(momReco, (pxReco - pxPart) / pxPart);
      hRelDeltaPyVsP->Fill(momReco, (pyReco - pyPart) / pyPart);
      hRelDeltaPzVsP->Fill(momReco, (pzReco - pzPart) / pzPart);
    }
  }

  TH1F* hPurityPt = (TH1F*)hMomGoodReco->Clone("hPurityPt");
  hPurityPt->GetYaxis()->SetTitle("purity");
  for (int iBin = 1; iBin <= hMomGoodReco->GetNbinsX(); iBin++) {
    double cg = hMomGoodReco->GetBinContent(iBin);
    double ct = hMomAllReco->GetBinContent(iBin);
    if (ct == 0) {
      hPurityPt->SetBinContent(iBin, 0.);
      hPurityPt->SetBinError(iBin, 0.);
    } else {
      double p = cg / ct;
      double ep = std::sqrt(p * (1 - p) / ct);
      hPurityPt->SetBinContent(iBin, p);
      hPurityPt->SetBinError(iBin, ep);
    }
  }
  TH1F* hPurityEta = (TH1F*)hEtaGoodReco->Clone("hPurityEta");
  hPurityEta->GetYaxis()->SetTitle("purity");
  for (int iBin = 1; iBin <= hEtaGoodReco->GetNbinsX(); iBin++) {
    double cg = hEtaGoodReco->GetBinContent(iBin);
    double ct = hEtaAllReco->GetBinContent(iBin);
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

  TH1F* hImpParXMean = new TH1F("hImpParXMean", ";p (GeV/c);<Imp Par X> (#mum)", nMomBins, 0., maxP);
  TH1F* hImpParYMean = new TH1F("hImpParYMean", ";p (GeV/c);<Imp Par Y> (#mum)", nMomBins, 0., maxP);
  TH1F* hImpParXRms = new TH1F("hImpParXRms", ";p (GeV/c);rms (Imp Par X) (#mum)", nMomBins, 0., maxP);
  TH1F* hImpParYRms = new TH1F("hImpParYRms", ";p (GeV/c);rms (Imp Par Y) (#mum)", nMomBins, 0., maxP);
  TH1F* hImpParXSig = new TH1F("hImpParXSig", ";p (GeV/c);#sigma(Imp Par X) (#mum)", nMomBins, 0., maxP);
  TH1F* hImpParYSig = new TH1F("hImpParYSig", ";p (GeV/c);#sigma(Imp Par Y) (#mum)", nMomBins, 0., maxP);
  fillMeanAndRms(hImpParXVsP, hImpParXMean, hImpParXRms, hImpParXSig);
  fillMeanAndRms(hImpParYVsP, hImpParYMean, hImpParYRms, hImpParYSig);

  TH1F* hDeltaPxMean = new TH1F("hDeltaPxMean", ";p (GeV/c);<p_{x}^{rec}-p_{x}^{gen}> (GeV/c)", nMomBins, 0., maxP);
  TH1F* hDeltaPxRms = new TH1F("hDeltaPxRms", ";p (GeV/c);rms (p_{x}^{rec}-p_{x}^{gen}) (GeV/c)", nMomBins, 0., maxP);
  TH1F* hDeltaPxSig = new TH1F("hDeltaPxSig", ";p (GeV/c);#sigma(p_{x}^{rec}-p_{x}^{gen}) (GeV/c)", nMomBins, 0., maxP);
  TH1F* hRelDeltaPxMean = new TH1F("hRelDeltaPxMean", ";p (GeV/c);<(p_{x}^{rec}-p_{x}^{gen})/p_{x}^{gen}>", nMomBins, 0., maxP);
  TH1F* hRelDeltaPxRms = new TH1F("hRelDeltaPxRms", ";p (GeV/c);rms (p_{x}^{rec}-p_{x}^{gen})/p_{x}^{gen}", nMomBins, 0., maxP);
  TH1F* hRelDeltaPxSig = new TH1F("hRelDeltaPxSig", ";p (GeV/c);#sigma(p_{x}^{rec}-p_{x}^{gen})/p_{x}^{gen} (GeV/c)", nMomBins, 0., maxP);
  fillMeanAndRms(hDeltaPxVsP, hDeltaPxMean, hDeltaPxRms, hDeltaPxSig);
  fillMeanAndRms(hRelDeltaPxVsP, hRelDeltaPxMean, hRelDeltaPxRms, hRelDeltaPxSig);

  TH1F* hDeltaPyMean = new TH1F("hDeltaPyMean", ";p (GeV/c);<p_{y}^{rec}-p_{y}^{gen}> (GeV/c)", nMomBins, 0., maxP);
  TH1F* hDeltaPyRms = new TH1F("hDeltaPyRms", ";p (GeV/c);rms (p_{y}^{rec}-p_{y}^{gen}) (GeV/c)", nMomBins, 0., maxP);
  TH1F* hDeltaPySig = new TH1F("hDeltaPySig", ";p (GeV/c);#sigma(p_{y}^{rec}-p_{y}^{gen}) (GeV/c)", nMomBins, 0., maxP);
  TH1F* hRelDeltaPyMean = new TH1F("hRelDeltaPyMean", ";p (GeV/c);<(p_{y}^{rec}-p_{y}^{gen})/p_{y}^{gen}>", nMomBins, 0., maxP);
  TH1F* hRelDeltaPyRms = new TH1F("hRelDeltaPyRms", ";p (GeV/c);rms (p_{y}^{rec}-p_{y}^{gen})/p_{y}^{gen}", nMomBins, 0., maxP);
  TH1F* hRelDeltaPySig = new TH1F("hRelDeltaPySig", ";p (GeV/c);#sigma(p_{y}^{rec}-p_{y}^{gen})/p_{y}^{gen} (GeV/c)", nMomBins, 0., maxP);
  fillMeanAndRms(hDeltaPyVsP, hDeltaPyMean, hDeltaPyRms, hDeltaPySig);
  fillMeanAndRms(hRelDeltaPyVsP, hRelDeltaPyMean, hRelDeltaPyRms, hRelDeltaPySig);

  TH1F* hDeltaPzMean = new TH1F("hDeltaPzMean", ";p (GeV/c);<p_{z}^{rec}-p_{z}^{gen}> (GeV/c)", nMomBins, 0., maxP);
  TH1F* hDeltaPzRms = new TH1F("hDeltaPzRms", ";p (GeV/c);rms (p_{z}^{rec}-p_{z}^{gen}) (GeV/c)", nMomBins, 0., maxP);
  TH1F* hDeltaPzSig = new TH1F("hDeltaPzSig", ";p (GeV/c);#sigma(p_{z}^{rec}-p_{z}^{gen}) (GeV/c)", nMomBins, 0., maxP);
  TH1F* hRelDeltaPzMean = new TH1F("hRelDeltaPzMean", ";p (GeV/c);<(p_{z}^{rec}-p_{z}^{gen})/p_{z}^{gen}>", nMomBins, 0., maxP);
  TH1F* hRelDeltaPzRms = new TH1F("hRelDeltaPzRms", ";p (GeV/c);rms (p_{z}^{rec}-p_{z}^{gen})/p_{z}^{gen}", nMomBins, 0., maxP);
  TH1F* hRelDeltaPzSig = new TH1F("hRelDeltaPzSig", ";p (GeV/c);#sigma(p_{z}^{rec}-p_{z}^{gen})/p_{z}^{gen} (GeV/c)", nMomBins, 0., maxP);
  fillMeanAndRms(hDeltaPzVsP, hDeltaPzMean, hDeltaPzRms, hDeltaPzSig);
  fillMeanAndRms(hRelDeltaPzVsP, hRelDeltaPzMean, hRelDeltaPzRms, hRelDeltaPzSig);

  TCanvas* cef = new TCanvas("cef", "Efficiency", 1400, 800);
  cef->Divide(2, 2);
  cef->cd(1);
  superposHistos(hMomGen, hMomTrackable, hMomTracked, kMagenta + 1, kBlue + 1, 1);
  cef->cd(2);
  TH1F* hEffMom = (TH1F*)hMomTracked->Clone("hEffMom");
  hEffMom->Divide(hMomTracked, hMomTrackable, 1., 1., "B");
  hEffMom->GetYaxis()->SetTitle("Efficiency");
  hEffMom->SetStats(0);
  hEffMom->Draw();
  cef->cd(3);
  superposHistos(hEtaGen, hEtaTrackable, hEtaTracked, kMagenta + 1, kBlue + 1, 1);
  cef->cd(4);
  TH1F* hEffEta = (TH1F*)hEtaTracked->Clone("hEffEta");
  hEffEta->Divide(hEtaTracked, hEtaTrackable, 1., 1., "B");
  hEffEta->GetYaxis()->SetTitle("Efficiency");
  hEffEta->SetStats(0);
  hEffEta->Draw();

  TCanvas* cpu = new TCanvas("cpu", "Purity", 1400, 800);
  cpu->Divide(2, 2);
  cpu->cd(1);
  gPad->SetLogy();
  superposHistos(hMomAllReco, hMomGoodReco, hMomFakeReco);
  cpu->cd(2);
  hPurityPt->GetYaxis()->SetTitle("Purity");
  hPurityPt->SetMinimum(0.8);
  hPurityPt->SetStats(0);
  hPurityPt->SetLineWidth(2);
  hPurityPt->Draw();
  cpu->cd(3);
  superposHistos(hEtaAllReco, hEtaGoodReco, hEtaFakeReco);
  cpu->cd(4);
  hPurityEta->GetYaxis()->SetTitle("Purity");
  hPurityEta->SetMinimum(0.8);
  hPurityEta->SetStats(0);
  hPurityEta->SetLineWidth(2);
  hPurityEta->Draw();

  TCanvas* cip = new TCanvas("cip", "Impact Parameter", 1200, 800);
  cip->Divide(2, 2);
  cip->cd(1);
  gPad->SetTickx();
  gPad->SetTicky();
  hImpParXMean->SetMarkerStyle(25);
  hImpParXMean->SetMarkerColor(1);
  hImpParXMean->SetLineColor(1);
  hImpParXMean->SetMinimum(-20);
  hImpParXMean->SetMaximum(20);
  hImpParXMean->Draw("P");
  cip->cd(2);
  gPad->SetTickx();
  gPad->SetTicky();
  hImpParXSig->SetMarkerStyle(25);
  hImpParXSig->SetMarkerColor(1);
  hImpParXSig->SetLineColor(1);
  hImpParXSig->SetMinimum(0);
  hImpParXSig->SetMaximum(120);
  hImpParXSig->Draw("P");
  cip->cd(3);
  gPad->SetTickx();
  gPad->SetTicky();
  hImpParYMean->SetMarkerStyle(25);
  hImpParYMean->SetMarkerColor(1);
  hImpParYMean->SetLineColor(1);
  hImpParYMean->SetMinimum(-20);
  hImpParYMean->SetMaximum(20);
  hImpParYMean->Draw("P");
  cip->cd(4);
  gPad->SetTickx();
  gPad->SetTicky();
  hImpParYSig->SetMarkerStyle(25);
  hImpParYSig->SetMarkerColor(1);
  hImpParYSig->SetLineColor(1);
  hImpParYSig->SetMinimum(0);
  hImpParYSig->SetMaximum(120);
  hImpParYSig->Draw("P");

  TCanvas* cmom = new TCanvas("cmom", "Momentum Resolution", 1400, 900);
  cmom->Divide(3, 3);
  cmom->cd(1);
  gPad->SetTickx();
  gPad->SetTicky();
  hDeltaPxMean->SetMarkerStyle(25);
  hDeltaPxMean->SetMarkerColor(1);
  hDeltaPxMean->SetLineColor(1);
  hDeltaPxMean->SetMinimum(-0.05);
  hDeltaPxMean->SetMaximum(0.05);
  hDeltaPxMean->Draw("P");
  cmom->cd(2);
  gPad->SetTickx();
  gPad->SetTicky();
  hDeltaPxSig->SetMarkerStyle(25);
  hDeltaPxSig->SetMarkerColor(1);
  hDeltaPxSig->SetLineColor(1);
  hDeltaPxSig->SetMinimum(0);
  hDeltaPxSig->SetMaximum(0.05);
  hDeltaPxSig->Draw("P");
  cmom->cd(3);
  gPad->SetTickx();
  gPad->SetTicky();
  hRelDeltaPxSig->SetMarkerStyle(25);
  hRelDeltaPxSig->SetMarkerColor(1);
  hRelDeltaPxSig->SetLineColor(1);
  hRelDeltaPxSig->SetMinimum(0);
  hRelDeltaPxSig->SetMaximum(0.1);
  hRelDeltaPxSig->Draw("P");
  cmom->cd(4);
  gPad->SetTickx();
  gPad->SetTicky();
  hDeltaPyMean->SetMarkerStyle(25);
  hDeltaPyMean->SetMarkerColor(1);
  hDeltaPyMean->SetLineColor(1);
  hDeltaPyMean->SetMinimum(-0.05);
  hDeltaPyMean->SetMaximum(0.05);
  hDeltaPyMean->Draw("P");
  cmom->cd(5);
  gPad->SetTickx();
  gPad->SetTicky();
  hDeltaPySig->SetMarkerStyle(25);
  hDeltaPySig->SetMarkerColor(1);
  hDeltaPySig->SetLineColor(1);
  hDeltaPySig->SetMinimum(0);
  hDeltaPySig->SetMaximum(0.05);
  hDeltaPySig->Draw("P");
  cmom->cd(6);
  gPad->SetTickx();
  gPad->SetTicky();
  hRelDeltaPySig->SetMarkerStyle(25);
  hRelDeltaPySig->SetMarkerColor(1);
  hRelDeltaPySig->SetLineColor(1);
  hRelDeltaPySig->SetMinimum(0);
  hRelDeltaPySig->SetMaximum(0.1);
  hRelDeltaPySig->Draw("P");
  cmom->cd(7);
  gPad->SetTickx();
  gPad->SetTicky();
  hDeltaPzMean->SetMarkerStyle(25);
  hDeltaPzMean->SetMarkerColor(1);
  hDeltaPzMean->SetLineColor(1);
  hDeltaPzMean->SetMinimum(-0.05);
  hDeltaPzMean->SetMaximum(0.05);
  hDeltaPzMean->Draw("P");
  cmom->cd(8);
  gPad->SetTickx();
  gPad->SetTicky();
  hDeltaPzSig->SetMarkerStyle(25);
  hDeltaPzSig->SetMarkerColor(1);
  hDeltaPzSig->SetLineColor(1);
  hDeltaPzSig->SetMinimum(0);
  hDeltaPzSig->SetMaximum(0.2);
  hDeltaPzSig->Draw("P");
  cmom->cd(9);
  hRelDeltaPzSig->SetMarkerStyle(25);
  hRelDeltaPzSig->SetMarkerColor(1);
  hRelDeltaPzSig->SetLineColor(1);
  hRelDeltaPzSig->SetMinimum(0);
  hRelDeltaPzSig->SetMaximum(0.05);
  hRelDeltaPzSig->Draw("P");

  TCanvas* c2d = new TCanvas("c2d", "", 1500, 500);
  c2d->Divide(3, 1);
  c2d->cd(1);
  hDeltaPxVsP->Draw("colz");
  c2d->cd(2);
  hDeltaPyVsP->Draw("colz");
  c2d->cd(3);
  hDeltaPzVsP->Draw("colz");

  TFile* outRoot = new TFile("TrackingPerformance.root", "recreate");
  hImpParXVsP->Write();
  hImpParYVsP->Write();
  hDeltaPxVsP->Write();
  hDeltaPyVsP->Write();
  hDeltaPzVsP->Write();
  hRelDeltaPxVsP->Write();
  hRelDeltaPyVsP->Write();
  hRelDeltaPzVsP->Write();
  hMomGen->Write();
  hEtaGen->Write();
  hMomTrackable->Write();
  hEtaTrackable->Write();
  hMomTracked->Write();
  hEtaTracked->Write();
  hMomAllReco->Write();
  hMomGoodReco->Write();
  hMomFakeReco->Write();
  hEtaAllReco->Write();
  hEtaGoodReco->Write();
  hEtaFakeReco->Write();
  outRoot->Close();
}
