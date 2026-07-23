#if !defined(__CINT__) || defined(__MAKECINT__)
#include <Riostream.h>
#include <TString.h>
#include <TMath.h>
#include <TTree.h>
#include <TLegend.h>
#include <TLegendEntry.h>
#include <TH1F.h>
#include <TF1.h>
#include <TH2F.h>
#include <TFile.h>
#include <TCanvas.h>
#include <TPad.h>
#include <TPaveStats.h>
#include <TParticle.h>
#include <TParticlePDG.h>
#include <TVectorD.h>
#include <TRandom3.h>
#include <TStopwatch.h>
#include "NA6PBaseCluster.h"
#include "NA6PVerTelCluster.h"
#include "NA6PVertex.h"
#include "NA6PMCTruthContainer.h"
#include "NA6PVertexerTracklets.h"
#endif

void inspectVertexerTracklets(int firstEv = 0,
                              int lastEv = 999999,
                              const char* dirSimu = ".")
// ../testN6Proot/pions/latesttag")
//			const char *dirSimu = "Angantyr")
{
  TFile* fk = new TFile(Form("%s/MCKine.root", dirSimu));
  TTree* mcTree = (TTree*)fk->Get("mckine");
  int nEv = mcTree->GetEntries();
  std::vector<TParticle>* mcArr = nullptr;
  mcTree->SetBranchAddress("tracks", &mcArr);

  TFile* fc = new TFile(Form("%s/ClustersVerTel.root", dirSimu));
  printf("Open cluster file: %s\n", fc->GetName());
  TTree* tc = (TTree*)fc->Get("clustersVerTel");
  std::vector<NA6PVerTelCluster> vtClus, *vtClusPtr = &vtClus;
  NA6PMCTruthContainer vtCluMCLabels, *vtCluMCLabelsPtr = &vtCluMCLabels;
  tc->SetBranchAddress("VerTel", &vtClusPtr);
  tc->SetBranchAddress("VerTelMCTruth", &vtCluMCLabelsPtr);

  if (lastEv > nEv || lastEv < 0)
    lastEv = nEv;
  if (firstEv < 0)
    firstEv = 0;

  TH1F* hdz = new TH1F("hdz", "", 100, -0.5, 0.5);
  TH1F* hncontr = new TH1F("hncontr", "", 100, -0.5, 999.5);
  TH1F* hnvert = new TH1F("hnvert", "", 11, -0.5, 10.5);
  TH1D* histocheck = new TH1D("histocheck", "; z_{intersection} (cm); counts", 250, -20., 5.);
  const int kMaxPileupVertices = NA6PVertexerTracklets::kMaxPileupVertices;
  TH1D* histocheck2[kMaxPileupVertices];
  for (int jPil = 0; jPil < kMaxPileupVertices; jPil++) {
    histocheck2[jPil] = new TH1D(Form("histocheckPil%d", jPil + 1), "", 250, -20., 5.);
  }
  TH1F* hSigmaZ = new TH1F("hSigmaZ", "", 100, 0., 10.);

  NA6PVertexerTracklets* vertxr = new NA6PVertexerTracklets();
  // configurations for the vertexer
  vertxr->setVerbosity(true);
  vertxr->configureFromRecoParam(/*"testpp.ini"*/);
  // vertxr->setUseKDEForPeakFinding();
  // vertxr->setMultiVertexInOneGo();
  // vertxr->setMaxDeltaPhiTracklet(0.05);
  // vertxr->configurePairClusterization(0.15, 0.2, 0.15, 0.2);
  for (int jEv = firstEv; jEv < lastEv; jEv++) {
    mcTree->GetEvent(jEv);
    tc->GetEvent(jEv);
    // set cluster indices
    for (size_t i = 0; i < vtClus.size(); ++i) {
      vtClus[i].setClusterIndex(static_cast<int>(i));
    }

    int nPart = mcArr->size();
    double xVertGen = 0;
    double yVertGen = 0;
    double zVertGen = 0;
    // get primary vertex position from the Kine Tree
    for (int jp = 0; jp < nPart; jp++) {
      auto curPart = mcArr->at(jp);
      if (curPart.IsPrimary()) {
        xVertGen = curPart.Vx();
        yVertGen = curPart.Vy();
        zVertGen = curPart.Vz();
        break;
      }
    }
    std::vector<TracklIntersection> zIntersec;
    int nClus = vtClus.size();
    printf("--- Process event %d, nClus = %d, primary vertex in %.2f %.2f %.2f cm ---\n", jEv, nClus, xVertGen, yVertGen, zVertGen);
    vertxr->setBeamX(xVertGen);
    vertxr->setBeamY(yVertGen);
    vertxr->resetClusters(nClus);
    vertxr->setClusterMCTruth(&vtCluMCLabels);
    std::vector<int> firstCluPerLay;
    std::vector<int> lastCluPerLay;
    vertxr->sortClustersByLayerAndEta(vtClus, firstCluPerLay, lastCluPerLay);
    std::vector<TrackletForVertex> allTracklets;
    std::vector<TrackletForVertex> selTracklets;
    vertxr->computeLayerTracklets(vtClus, firstCluPerLay, lastCluPerLay, allTracklets);
    std::vector<int> firstTrklPerLay;
    std::vector<int> lastTrklPerLay;
    vertxr->sortTrackletsByLayerAndIndex(allTracklets, firstTrklPerLay, lastTrklPerLay);
    vertxr->printStats(allTracklets, vtClus, "tracklets");
    vertxr->selectTracklets(allTracklets, firstTrklPerLay, lastTrklPerLay, selTracklets);
    vertxr->printStats(selTracklets, vtClus, "selected tracklets");
    //
    std::vector<NA6PVertex> zVertices;
    bool retCode;
    if (vertxr->getMethodForPeakFinding() == NA6PVertexerTracklets::kPairs) {
      retCode = vertxr->compute3DVertices(selTracklets, vtClus, zVertices);
    } else {
      vertxr->computeIntersections(selTracklets, vtClus, zIntersec);
      if (vertxr->getMethodForPeakFinding() == NA6PVertexerTracklets::kKDE)
        retCode = vertxr->findVertexKDE(zIntersec, zVertices);
      else
        retCode = vertxr->findVertexHistoPeak(zIntersec, zVertices);
      histocheck->Reset("MICE");
      for (auto zi : zIntersec) {
        histocheck->Fill(zi.zeta);
        hSigmaZ->Fill(zi.sigmazeta);
      }
    }
    int nVertices = zVertices.size();
    int jv = 0;
    for (auto vert : zVertices) {
      double zRec = vert.getZ();
      printf("Vertex %d, z = %f contrib = %d\n", jv++, zRec, vert.getNContributors());
      hdz->Fill(zRec - zVertGen);
      hncontr->Fill(vert.getNContributors());
    }
    // pileup detection
    if (vertxr->getMultiVertexMode() == NA6PVertexerTracklets::kMultiVertIterative) {
      std::vector<TrackletForVertex> remainingTracklets;
      for (int jPil = 0; jPil < kMaxPileupVertices; jPil++) {
        vertxr->filterOutUsedTracklets(selTracklets, remainingTracklets);
        vertxr->printStats(remainingTracklets, vtClus, "remaining tracklets");
        if (vertxr->getMethodForPeakFinding() == NA6PVertexerTracklets::kPairs) {
          retCode = vertxr->compute3DVertices(remainingTracklets, vtClus, zVertices);
        } else {
          vertxr->computeIntersections(remainingTracklets, vtClus, zIntersec);
          if (vertxr->getMethodForPeakFinding() == NA6PVertexerTracklets::kKDE)
            retCode = vertxr->findVertexKDE(zIntersec, zVertices);
          else
            retCode = vertxr->findVertexHistoPeak(zIntersec, zVertices);
          histocheck2[jPil]->Reset("MICES");
          for (auto zi : zIntersec)
            histocheck2[jPil]->Fill(zi.zeta);
        }
        if (!retCode)
          break;
        selTracklets.swap(remainingTracklets);
      }
    }
    nVertices = zVertices.size();
    hnvert->Fill(nVertices);
    printf("After pileup search: nVertices = %d\n", nVertices);
    jv = 0;
    for (auto vert : zVertices) {
      double zRec = vert.getZ();
      printf("Vertex %d, z = %f nContrib = %d\n", jv++, zRec, vert.getNContributors());
    }
  }

  TCanvas* cev = new TCanvas("cev", "", 1200, 500);
  histocheck->SetLineWidth(2);
  histocheck->Draw();
  for (int jPil = 0; jPil < kMaxPileupVertices; jPil++) {
    if (histocheck2[jPil]->GetEntries() > 0) {
      histocheck2[jPil]->SetLineColor(jPil + 2);
      histocheck2[jPil]->SetLineWidth(2);
      histocheck2[jPil]->Draw("sames");
    }
  }

  TCanvas* coutp = new TCanvas("coutp", "", 1200, 800);
  coutp->Divide(2, 2);
  coutp->cd(1);
  hnvert->Draw();
  coutp->cd(2);
  hdz->Draw();
  coutp->cd(3);
  hncontr->Draw();
  coutp->cd(4);
  hSigmaZ->Draw();
}
