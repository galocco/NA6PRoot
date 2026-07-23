#include <TFile.h>
#include <TTree.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TStyle.h>
#include <TProfile.h>
#include <TCanvas.h>
#include <TParticle.h>
#include "NA6PVerTelHit.h"
#include "NA6PMCTruthContainer.h"
#include "NA6PVerTelDigit.h"
#include "NA6PVerTelSegmentation.h"
#include "NA6PGeometryManager.h"

void getGlobalCoord(const NA6PGeometryManager& geoMan, int modID, float xloc, float yloc, double xyzGlo[3])
{
  auto& matrix = geoMan.getMatrix(modID);

  double xyzLoc[3];
  xyzLoc[0] = xloc - geoMan.getModuleHalfX(modID);
  xyzLoc[1] = yloc - geoMan.getModuleHalfY(modID);
  xyzLoc[2] = 0.;
  matrix.LocalToMaster(xyzLoc, xyzGlo);
}

void plotVerTelDigits(const char* dirSimu = ".")
{

  NA6PGeometryManager geoMan;
  geoMan.loadGeometry(Form("%s/geometry.root", dirSimu));
  NA6PVerTelSegmentation vt;

  TFile* fk = new TFile(Form("%s/MCKine.root", dirSimu));
  TTree* mcTree = (TTree*)fk->Get("mckine");
  std::vector<TParticle>* mcArr = nullptr;
  mcTree->SetBranchAddress("tracks", &mcArr);

  TFile* fh = new TFile(Form("%s/HitsVerTel.root", dirSimu));
  TTree* th = (TTree*)fh->Get("hitsVerTel");
  std::vector<NA6PVerTelHit> vtHits, *vtHitsPtr = &vtHits;
  th->SetBranchAddress("VerTel", &vtHitsPtr);

  TFile* fd = TFile::Open(Form("%s/DigitsVerTel.root", dirSimu));
  TTree* td = (TTree*)fd->Get("digitsVerTel");
  std::vector<NA6PVerTelDigit> vtDigits, *vtDigitsPtr = &vtDigits;
  NA6PMCTruthContainer vtMCLabels, *vtMCLabelsPtr = &vtMCLabels;
  td->SetBranchAddress("VerTel", &vtDigitsPtr);
  td->SetBranchAddress("VerTelMCTruth", &vtMCLabelsPtr);

  TH1F* hMCLabelsPerDigit = new TH1F("hMCLabelsPerDigit", ";n MC labels per digit; counts", 11, -0.5, 10.5);
  TH1F* hMCFirstLabels = new TH1F("hMCFirstLabels", "First label only;MC label; counts", 200, 0., 4000.);
  TH1F* hMCAllLabels = new TH1F("hMCAllLabels", "All labels;MC label; counts", 200, 0., 4000.);
  TH1F* hDigitsPerHitAll = new TH1F("hDigitsPerHitAll", "All particles;n digits per hit; counts", 11, -0.5, 10.5);
  TH1F* hDigitsPerHitPrim = new TH1F("hDigitsPerHitPrim", "Primary particles;n digits per hit; counts", 11, -0.5, 10.5);
  TH2F** hDigitsPerHitPrimVsX = new TH2F*[5];
  TH2F** hDigitsPerHitPrimVsY = new TH2F*[5];
  TProfile* pDigitsPerHitPrimVsX[5];
  TProfile* pDigitsPerHitPrimVsY[5];
  for (int jLay = 0; jLay < 5; ++jLay) {
    hDigitsPerHitPrimVsX[jLay] = new TH2F(Form("hDigitsPerHitPrimVsX%d", jLay), "Primary particles;x (cm);n digits per hit; counts", 30, -15., 15., 8, -0.5, 7.5);
    hDigitsPerHitPrimVsY[jLay] = new TH2F(Form("hDigitsPerHitPrimVsY%d", jLay), "Primary particles;y (cm);n digits per hit; counts", 30, -15., 15., 8, -0.5, 7.5);
  }

  TH1F* hRSU = new TH1F("hRSU", ";rsu;counts", 42, -0.5, 41.5);
  TH1F* hTile = new TH1F("hTile", ";tile;counts", 12, -0.5, 11.5);
  TH1F* hRow = new TH1F("hRow", ";row;counts", 444, -0.5, 443.5);
  TH1F* hCol = new TH1F("hCol", ";col;counts", 156, -0.5, 155.5);
  TH1F* hDeltaXloc1 = new TH1F("hDeltaXloc1", "1-digit hits;x_{loc}^{digit}-x_{loc}^{hit} (#mum);counts", 100, -30., 30.);
  TH1F* hDeltaYloc1 = new TH1F("hDeltaYloc1", "1-digit hits;y_{loc}^{digit}-y_{loc}^{hit} (#mum);counts", 100, -30., 30.);
  TH1F* hDeltaXloc2 = new TH1F("hDeltaXloc2", "2-digit hits;x_{loc}^{digit}-x_{loc}^{hit} (#mum);counts", 100, -30., 30.);
  TH1F* hDeltaYloc2 = new TH1F("hDeltaYloc2", "2-digit hits;y_{loc}^{digit}-y_{loc}^{hit} (#mum);counts", 100, -30., 30.);
  TH2F** hXYHit = new TH2F*[5];
  TH1F** hZHit = new TH1F*[5];
  TH2F** hXYGloDig = new TH2F*[5];
  TH1F** hZGloDig = new TH1F*[5];
  for (int jLay = 0; jLay < 5; jLay++) {
    hXYHit[jLay] = new TH2F(Form("hXYHitLay%d", jLay), Form("hits layer %d;x_{glo} (cm);y_{glo} (cm)", jLay), 100, -15., 15., 100, -15., 15.);
    hZHit[jLay] = new TH1F(Form("hZHitLay%d", jLay), Form("hits layer %d;z_{glo} (cm);counts", jLay), 100, 0., 40.);
    hXYGloDig[jLay] = new TH2F(Form("hXYGloDigLay%d", jLay), Form("digits layer %d;x_{glo} (cm);y_{glo} (cm)", jLay), 100, -15., 15., 100, -15., 15.);
    hZGloDig[jLay] = new TH1F(Form("hZGloDigLay%d", jLay), Form("digits layer %d;z_{glo} (cm);counts", jLay), 100, 0., 40.);
    hXYHit[jLay]->SetStats(0);
    hZHit[jLay]->SetStats(0);
    hXYGloDig[jLay]->SetStats(0);
    hZGloDig[jLay]->SetStats(0);
  }
  TH2F** hXYLocDig = new TH2F*[20];
  for (int jMod = 0; jMod < 20; jMod++) {
    hXYLocDig[jMod] = new TH2F(Form("hXYLocDigMod%d", jMod), Form("digits mod %d;x_{loc} (cm);y_{loc} (cm)", jMod), 100, 0., 15., 100, 0., 15.);
    hXYLocDig[jMod]->SetStats(0);
  }

  int nEv = mcTree->GetEntries();
  printf("Number of events = %d\n", nEv);
  int countEx5 = 0;
  for (int jEv = 0; jEv < nEv; jEv++) {
    mcTree->GetEvent(jEv);
    th->GetEvent(jEv);
    td->GetEvent(jEv);
    int nPart = mcArr->size();
    int nHits = vtHits.size();
    int nDigits = vtDigits.size();
    int nMClabels = vtMCLabels.getNElements();
    printf("Event %d particles = %d hits = %d digits = %d MClabels = %d\n", jEv, nPart, nHits, nDigits, nMClabels);
    for (const auto& hit : vtHits) {
      int nDet = hit.getDetectorID();
      int lay = nDet / 4;
      if (lay < 0 || lay >= 5) {
        printf("ERROR wrong layer %d\n", lay);
        continue;
      }
      hXYHit[lay]->Fill(hit.getX(), hit.getY());
      hZHit[lay]->Fill(hit.getZ());
    }
    for (int jDig = 0; jDig < nDigits; ++jDig) {
      const auto& dig = vtDigits.at(jDig);
      int rsu = dig.getRSU();
      int tile = dig.getTile();
      int row = dig.getRow();
      int col = dig.getCol();
      int modId = dig.getDetectorID();
      int lay = modId / 4;
      hRSU->Fill(rsu);
      hTile->Fill(tile);
      hRow->Fill(row);
      hCol->Fill(col);
      float xpix, ypix;
      bool locOk = vt.indicesToLocal(rsu, tile, row, col, xpix, ypix);
      if (locOk) {
        hXYLocDig[modId]->Fill(xpix, ypix);
        double xyzGlo[3];
        getGlobalCoord(geoMan, dig.getDetectorID(), xpix, ypix, xyzGlo);
        hXYGloDig[lay]->Fill(xyzGlo[0], xyzGlo[1]);
        hZGloDig[lay]->Fill(xyzGlo[2]);
      }

      std::span labels = vtMCLabels.getLabels(jDig);
      int nLabels = labels.size();
      hMCLabelsPerDigit->Fill(nLabels);
      for (int jLab = 0; jLab < nLabels; jLab++) {
        NA6PMCComposedLabel lbl = labels[jLab];
        if (jLab == 0)
          hMCFirstLabels->Fill(lbl.getTrackID());
        hMCAllLabels->Fill(lbl.getTrackID());
      }
    }
    for (int jHit = 0; jHit < nHits; ++jHit) {
      const auto& hit = vtHits.at(jHit);
      int idPartH = hit.getTrackID();
      auto curPart = mcArr->at(idPartH);
      int sensH = hit.getDetectorID();
      int digperhit = 0;
      float sumx = 0.;
      float sumy = 0.;
      float sumgooddigs = 0.;
      for (int jDig = 0; jDig < nDigits; ++jDig) {
        const auto& dig = vtDigits.at(jDig);
        int sensD = dig.getDetectorID();
        if (sensH == sensD) {
          std::span labels = vtMCLabels.getLabels(jDig);
          int nLabels = labels.size();
          for (int jLab = 0; jLab < nLabels; jLab++) {
            NA6PMCComposedLabel lbl = labels[jLab];
            int idPartD = lbl.getTrackID();
            if (idPartD == idPartH) {
              digperhit++;
              int rsu = dig.getRSU();
              int tile = dig.getTile();
              int row = dig.getRow();
              int col = dig.getCol();
              float xpix, ypix;
              bool locOk = vt.indicesToLocal(rsu, tile, row, col, xpix, ypix);
              if (locOk) {
                sumgooddigs += 1.;
                sumx += xpix;
                sumy += ypix;
              }
            }
          }
        }
      }
      hDigitsPerHitAll->Fill(digperhit);
      float xdig = sumx / sumgooddigs;
      float ydig = sumy / sumgooddigs;
      auto modID = hit.getDetectorID();
      auto& matrix = geoMan.getMatrix(modID);
      double xyzGloS[3] = {hit.getXIn(), hit.getYIn(), hit.getZIn()};
      double xyzGloE[3] = {hit.getXOut(), hit.getYOut(), hit.getZOut()};
      double xyzLocS[3], xyzLocE[3];
      matrix.MasterToLocal(xyzGloS, xyzLocS);
      xyzLocS[0] += geoMan.getModuleHalfX(modID);
      xyzLocS[1] += geoMan.getModuleHalfY(modID);
      matrix.MasterToLocal(xyzGloE, xyzLocE);
      xyzLocE[0] += geoMan.getModuleHalfX(modID);
      xyzLocE[1] += geoMan.getModuleHalfY(modID);
      double xhit = 0.5f * (xyzLocS[0] + xyzLocE[0]);
      double yhit = 0.5f * (xyzLocS[1] + xyzLocE[1]);
      if (digperhit == 1) {
        hDeltaXloc1->Fill((xdig - xhit) * 1e4);
        hDeltaYloc1->Fill((ydig - yhit) * 1e4);
      } else if (digperhit == 2) {
        hDeltaXloc2->Fill((xdig - xhit) * 1e4);
        hDeltaYloc2->Fill((ydig - yhit) * 1e4);
      }

      int layH = sensH / 4;
      if (curPart.IsPrimary()) {
        hDigitsPerHitPrim->Fill(digperhit);
        hDigitsPerHitPrimVsX[layH]->Fill(hit.getX(), digperhit);
        hDigitsPerHitPrimVsY[layH]->Fill(hit.getY(), digperhit);
      }
    }
  }

  gStyle->SetOptStat(11111111);

  TCanvas* cdig = new TCanvas("cdig", "", 1400, 800);
  cdig->Divide(2, 2);
  cdig->cd(1);
  hRSU->Draw();
  cdig->cd(2);
  hTile->Draw();
  cdig->cd(3);
  hRow->Draw();
  cdig->cd(4);
  hCol->Draw();

  TCanvas* cxyh = new TCanvas("cxyh", "Hits XY", 1400, 800);
  cxyh->Divide(3, 2);
  for (int jLay = 0; jLay < 5; jLay++) {
    cxyh->cd(jLay + 1);
    hXYHit[jLay]->Draw("colz");
  }
  TCanvas* czh = new TCanvas("czh", "Hits Z", 1400, 800);
  czh->Divide(3, 2);
  for (int jLay = 0; jLay < 5; jLay++) {
    czh->cd(jLay + 1);
    hZHit[jLay]->Draw();
  }
  TCanvas* cxyl = new TCanvas("cxyl", "Digits XY", 1400, 800);
  cxyl->Divide(5, 4);
  for (int jMod = 0; jMod < 20; jMod++) {
    int jLay = jMod / 4;
    int j4 = jMod % 4;
    int jPad = jLay + 1 + j4 * 5;
    cxyl->cd(jPad);
    hXYLocDig[jMod]->Draw("colz");
  }
  TCanvas* cxyg = new TCanvas("cxyg", "Digits XY", 1400, 800);
  cxyg->Divide(3, 2);
  for (int jLay = 0; jLay < 5; jLay++) {
    cxyg->cd(jLay + 1);
    hXYGloDig[jLay]->Draw("colz");
  }
  TCanvas* czg = new TCanvas("czg", "Digits Z", 1400, 800);
  czg->Divide(3, 2);
  for (int jLay = 0; jLay < 5; jLay++) {
    czg->cd(jLay + 1);
    hZGloDig[jLay]->Draw("colz");
  }

  TCanvas* clab = new TCanvas("clab", "labels", 1400, 500);
  clab->Divide(3, 1);
  clab->cd(1);
  gPad->SetLogy();
  hMCLabelsPerDigit->Draw();
  clab->cd(2);
  hMCFirstLabels->Draw();
  clab->cd(3);
  hMCAllLabels->Draw();

  TCanvas* chd = new TCanvas("chd", "", 1400, 800);
  chd->Divide(2, 1);
  chd->cd(1);
  gPad->SetLogy();
  hDigitsPerHitAll->Draw();
  chd->cd(2);
  gPad->SetLogy();
  hDigitsPerHitPrim->Draw();

  TCanvas* cres = new TCanvas("cres", "", 1400, 800);
  cres->Divide(2, 2);
  cres->cd(1);
  hDeltaXloc1->Draw();
  cres->cd(2);
  hDeltaYloc1->Draw();
  cres->cd(3);
  hDeltaXloc2->Draw();
  cres->cd(4);
  hDeltaYloc2->Draw();

  TCanvas* clay = new TCanvas("clay", "", 1400, 800);
  clay->Divide(5, 2);
  for (int jLay = 0; jLay < 5; ++jLay) {
    pDigitsPerHitPrimVsX[jLay] = hDigitsPerHitPrimVsX[jLay]->ProfileX(Form("profileXlay%d\n", jLay));
    pDigitsPerHitPrimVsY[jLay] = hDigitsPerHitPrimVsY[jLay]->ProfileX(Form("profileYlay%d\n", jLay));
    clay->cd(jLay + 1);
    gPad->SetLogz();
    hDigitsPerHitPrimVsX[jLay]->SetStats(0);
    hDigitsPerHitPrimVsX[jLay]->Draw("colz");
    pDigitsPerHitPrimVsX[jLay]->SetMarkerStyle(20);
    pDigitsPerHitPrimVsX[jLay]->SetMarkerSize(0.5);
    pDigitsPerHitPrimVsX[jLay]->SetLineWidth(2);
    pDigitsPerHitPrimVsX[jLay]->Draw("psame");
    clay->cd(jLay + 6);
    gPad->SetLogz();
    hDigitsPerHitPrimVsY[jLay]->SetStats(0);
    hDigitsPerHitPrimVsY[jLay]->Draw("colz");
    pDigitsPerHitPrimVsY[jLay]->SetMarkerStyle(20);
    pDigitsPerHitPrimVsY[jLay]->SetMarkerSize(0.5);
    pDigitsPerHitPrimVsY[jLay]->SetLineWidth(2);
    pDigitsPerHitPrimVsY[jLay]->Draw("psame");
  }
}
