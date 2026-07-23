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
#include <TTree.h>
#include <TFile.h>
#include <TCanvas.h>
#include <TPad.h>
#include <TPaveStats.h>
#include <TParticle.h>
#include <TParticlePDG.h>
#include <TVectorD.h>
#include <TRandom3.h>
#include <TStopwatch.h>
#include "NA6PLayoutParam.h"
#include "NA6PVerTelCluster.h"
#include "NA6PVertex.h"
#include "NA6PVerTelReconstruction.h"
#include "NA6PVertexerTracklets.h"
#include "NA6PMCEventHeader.h"
#include "NA6PMCTruthContainer.h"
#include "NA6PMCGenHeader.h"
#endif

void runVertexerTracklets(int firstEv = 0,
                          int lastEv = 10,
                          const char* dirSimu = ".",
                          const char* na6pLayoutFile = "./na6pLayout.ini",
                          const char* fOutName = "vtx_fixed_new_test.root")
{
  TFile* fk = new TFile(Form("%s/MCKine.root", dirSimu));
  TTree* mcTree = (TTree*)fk->Get("mckine");
  int nEv = mcTree->GetEntries();
  std::vector<TParticle>* mcArr = nullptr;
  NA6PMCEventHeader* mcHead = nullptr;
  mcTree->SetBranchAddress("tracks", &mcArr);
  mcTree->SetBranchAddress("header", &mcHead);

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

  int nVtx, targId;
  double genX, genY, genZ;
  int nContrib[20];
  double recsX[20], recsY[20], recsZ[20];

  TFile* fOut = new TFile(fOutName, "RECREATE");
  TTree* tree = new TTree("vertex", "vertex");
  tree->Branch("nVtx", &nVtx, "nVtx/I");
  tree->Branch("targId", &targId, "targId/I");
  tree->Branch("genX", &genX, "genX/D");
  tree->Branch("genY", &genY, "genY/D");
  tree->Branch("genZ", &genZ, "genZ/D");
  tree->Branch("nContrib", nContrib, "nContrib[20]/I");
  tree->Branch("recsX", recsX, "recsX[20]/D");
  tree->Branch("recsY", recsY, "recsY[20]/D");
  tree->Branch("recsZ", recsZ, "recsZ[20]/D");

  TH2F* hxgenrec = new TH2F("hxgenrec", ";x_{gen} (cm);x_{rec} (cm)", 1000, -5, 5, 1000, -5, 5);
  TH2F* hygenrec = new TH2F("hygenrec", ";y_{gen} (cm);y_{rec} (cm)", 1000, -5, 5, 1000, -5, 5);
  TH2F* hzgenrec = new TH2F("hzgenrec", ";z_{gen} (cm);z_{rec} (cm)", 1000, -9, 1, 1000, -9, 1);
  TH1F* hxgen = new TH1F("hxgen", "; x_{gen} (cm); counts", 1000, -5, 5);
  TH1F* hygen = new TH1F("hygen", "; y_{gen} (cm); counts", 1000, -5, 5);
  TH1F* hzgen = new TH1F("hzgen", "; x_{gen} (cm); counts", 1000, -9, 1);
  TH1F* hxrec = new TH1F("hxrec", "; x_{rec} (cm); counts", 1000, -5, 5);
  TH1F* hyrec = new TH1F("hyrec", "; y_{rec} (cm); counts", 1000, -5, 5);
  TH1F* hzrec = new TH1F("hzrec", "; x_{rec} (cm); counts", 1000, -9, 1);
  TH1F* hdx = new TH1F("hdx", "; x_{rec} - x_{gen} (cm); counts", 100, -0.5, 0.5);
  TH1F* hdy = new TH1F("hdy", "; y_{rec} - y_{gen} (cm); counts", 100, -0.5, 0.5);
  TH1F* hdz = new TH1F("hdz", "; z_{rec} - z_{gen} (cm); counts", 100, -0.5, 0.5);
  TH1F* hncontr = new TH1F("hncontr", "; N_{contributors}; counts", 100, -0.5, 999.5);
  TH1F* hnvert = new TH1F("hnvert", "; N_{vertices}; counts", 11, -0.5, 10.5);
  TH2F* hnvertcontrib = new TH2F("hnvertcontrib", "; N_{vertices}; N_{contributors}", 11, -0.5, 10.5, 100, -0.5, 999.5);

  TH1F* hzgens[5]; // Generated vertex not primary
  TH1F *hxrecs[5], *hyrecs[5], *hzrecs[5];
  TH1F *hdxs[5], *hdys[5], *hdzs[5];
  for (int i = 0; i < 5; i++) {
    hzgens[i] = new TH1F(Form("hzgen_target_%i", i), "; z_{rec} cm; counts", 1000, -9, 1);

    hxrecs[i] = new TH1F(Form("hxrec_target_%i", i), "; x_{rec} cm; counts", 1000, -5, 5);
    hyrecs[i] = new TH1F(Form("hyrec_target_%i", i), "; y_{rec} cm; counts", 1000, -5, 5);
    hzrecs[i] = new TH1F(Form("hzrec_target_%i", i), "; z_{rec} cm; counts", 1000, -9, 1);

    hdxs[i] = new TH1F(Form("hdx_target_%i", i), "; x_{rec} - x_{gen} cm; counts", 100, -0.5, 0.5);
    hdys[i] = new TH1F(Form("hdy_target_%i", i), "; y_{rec} - y_{gen} cm; counts", 100, -0.5, 0.5);
    hdzs[i] = new TH1F(Form("hdz_target_%i", i), "; z_{rec} - z_{gen} cm; counts", 100, -0.5, 0.5);
  }

  int indexTarg = -999;

  NA6PVerTelReconstruction* vtrec = new NA6PVerTelReconstruction();
  vtrec->initVertexer();
  vtrec->getVertexerTracklets()->setVerbosity(true);
  for (int jEv = firstEv; jEv < lastEv; jEv++) {
    mcTree->GetEvent(jEv);
    std::cout << "NColl: " << mcHead->getNColl() << " NPart: " << mcHead->getNPart() << " Impact Parameter: " << mcHead->getImpPar() << std::endl;
    tc->GetEvent(jEv);
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

        hxgen->Fill(xVertGen);
        hygen->Fill(yVertGen);
        hzgen->Fill(zVertGen);
        genX = xVertGen;
        genY = yVertGen;
        genZ = zVertGen;
        break;
      }
    }

    na6p::conf::ConfigurableParam::updateFromFile(na6pLayoutFile, "", true);
    const auto& layoutPar = NA6PLayoutParam::Instance();
    int nTargs = int(layoutPar.nTargets);

    // Find the position of the primary target
    for (int iTarg = 0; iTarg < nTargs; iTarg++) {
      double minTargZ = layoutPar.posTargetZ[iTarg] - ((layoutPar.thicknessTarget[iTarg]) / 2.);
      double maxTargZ = layoutPar.posTargetZ[iTarg] + ((layoutPar.thicknessTarget[iTarg]) / 2.);
      if (zVertGen > minTargZ && zVertGen < maxTargZ) {
        indexTarg = iTarg;
      }
    }
    targId = indexTarg;

    // Check if there are other interactions
    std::vector<int> tmpTargs;
    for (int i = 0; i < 5; i++) {
      if (i != indexTarg) {
        tmpTargs.push_back(i);
      }
    }

    for (int jp = 0; jp < nPart; jp++) {
      auto curPart = mcArr->at(jp);
      double zGenTmp = curPart.Vz();
      if (!curPart.IsPrimary() && zGenTmp < 1 && zGenTmp > -9) {
        for (int i = tmpTargs.size() - 1; i >= 0; --i) {
          int tmpTarg = tmpTargs[i];
          double minTargZ = layoutPar.posTargetZ[tmpTarg] - ((layoutPar.thicknessTarget[tmpTarg]) / 2.);
          double maxTargZ = layoutPar.posTargetZ[tmpTarg] + ((layoutPar.thicknessTarget[tmpTarg]) / 2.);
          if (zGenTmp > minTargZ && zGenTmp < maxTargZ) {
            // std::cout << "vertex in target " << tmpTarg << " found" << std::endl;
            hzgens[indexTarg]->Fill(zGenTmp);
            tmpTargs.erase(tmpTargs.begin() + i);
          }
        }
      }
    }

    vtrec->setClusters(vtClus);
    vtrec->setClustersMCLabels(vtCluMCLabels);
    vtrec->runVertexerTracklets();
    std::vector<NA6PVertex> zVertices = vtrec->getVertices();
    int nVertices = zVertices.size();
    nVtx = nVertices;
    hnvert->Fill(nVertices);
    int jv = 0;
    std::cout << "*********** N. vertices: " << nVertices << std::endl;
    for (auto vert : zVertices) {
      double xRec = vert.getX();
      double yRec = vert.getY();
      double zRec = vert.getZ();
      if (jv == 0) {
        hxrec->Fill(xRec);
        hyrec->Fill(yRec);
        hzrec->Fill(zRec);

        hxgenrec->Fill(xVertGen, xRec);
        hygenrec->Fill(yVertGen, yRec);
        hzgenrec->Fill(zVertGen, zRec);

        hdx->Fill(xRec - xVertGen);
        hdy->Fill(yRec - yVertGen);
        hdz->Fill(zRec - zVertGen);

        hxrecs[indexTarg]->Fill(xRec);
        hyrecs[indexTarg]->Fill(yRec);
        hzrecs[indexTarg]->Fill(zRec);

        hdxs[indexTarg]->Fill(xRec - xVertGen);
        hdys[indexTarg]->Fill(yRec - yVertGen);
        hdzs[indexTarg]->Fill(zRec - zVertGen);

        hncontr->Fill(vert.getNContributors());
      }
      hnvertcontrib->Fill(jv, vert.getNContributors());

      nContrib[jv] = vert.getNContributors();
      recsX[jv] = xRec;
      recsY[jv] = yRec;
      recsZ[jv] = zRec;

      printf("Vertex %d, z = %f contrib = %d\n", jv++, zRec, vert.getNContributors());
    }
    tree->Fill();
  }
  vtrec->closeVerticesOutput();

  TCanvas* coutp = new TCanvas("coutp", "", 1200, 600);
  coutp->Divide(3, 1);
  coutp->cd(1);
  hnvert->SetLineWidth(2);
  hnvert->Draw();
  coutp->cd(3);
  hdz->SetLineWidth(2);
  hdz->Draw();
  coutp->cd(2);
  hncontr->SetLineWidth(2);
  hncontr->Draw();

  // TFile *fOut = new TFile(Form("%s/%s", dirSimu, fOutName), "RECREATE");
  fOut->cd();
  coutp->Write();
  hnvert->Write();
  hxgen->Write();
  hygen->Write();
  hzgen->Write();
  hxrec->Write();
  hyrec->Write();
  hzrec->Write();
  hxgenrec->Write();
  hygenrec->Write();
  hzgenrec->Write();
  hdx->Write();
  hdy->Write();
  hdz->Write();

  for (int i = 0; i < 5; i++) {
    hzgens[i]->Write();
    hxrecs[i]->Write();
    hyrecs[i]->Write();
    hzrecs[i]->Write();
    hdxs[i]->Write();
    hdys[i]->Write();
    hdzs[i]->Write();
  }

  hncontr->Write();
  hnvertcontrib->Write();
  tree->Write();
  fOut->Close();
}
