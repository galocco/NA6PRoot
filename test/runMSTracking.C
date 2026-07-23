#if !defined(__CINT__) || defined(__MAKECINT__)
#include <TTree.h>
#include <TFile.h>
#include <TParticle.h>
#include <TStopwatch.h>
#include "NA6PMuonSpecCluster.h"
#include "NA6PMuonSpecReconstruction.h"
#include "Propagator.h"
#endif

void runMSTracking(int firstEv = 0,
                   int lastEv = 99999,
                   const char* dirSimu = ".")
{

  // kine file needed to get the primary vertex position (temporary)
  TFile* fk = new TFile(Form("%s/MCKine.root", dirSimu));
  TTree* mcTree = (TTree*)fk->Get("mckine");
  std::vector<TParticle>* mcArr = nullptr;
  mcTree->SetBranchAddress("tracks", &mcArr);

  TFile* fc = new TFile(Form("%s/ClustersMuonSpec.root", dirSimu));
  printf("Open cluster file: %s\n", fc->GetName());
  TTree* tc = (TTree*)fc->Get("clustersMuonSpec");
  std::vector<NA6PMuonSpecCluster> msClus, *msClusPtr = &msClus;
  NA6PMCTruthContainer msCluMCLabels, *msCluMCLabelsPtr = &msCluMCLabels;
  tc->SetBranchAddress("MuonSpec", &msClusPtr);
  tc->SetBranchAddress("MuonSpecMCTruth", &msCluMCLabelsPtr);

  if (!Propagator::loadField() || !Propagator::loadGeometry(Form("%s/geometry.root", dirSimu))) {
    return;
  }

  // if needed, configure the recoparams like
  // NA6PRecoParam::setValue("reco.msMaxChi2ndfTracksCA[0]", std::to_string(5.));
  NA6PMuonSpecReconstruction* msrec = new NA6PMuonSpecReconstruction();
  msrec->initTracker();
  auto msTracker = msrec->getTracker();
  int nEv = tc->GetEntries();
  if (lastEv > nEv || lastEv < 0)
    lastEv = nEv;
  if (firstEv < 0)
    firstEv = 0;

  TStopwatch timer;
  timer.Start();

  for (int jEv = firstEv; jEv < lastEv; jEv++) {
    msrec->clearEvent();
    mcTree->GetEvent(jEv);
    int nPart = mcArr->size();
    double zvert = 0;
    // get primary vertex position from the Kine Tree
    for (int jp = 0; jp < nPart; jp++) {
      auto curPart = mcArr->at(jp);
      if (curPart.IsPrimary()) {
        zvert = curPart.Vz();
        break;
      }
    }
    tc->GetEvent(jEv);
    NA6PVertex pvert;
    pvert.setXYZ(0., 0., zvert);
    msrec->setPrimaryVertex(&pvert);
    msrec->setClusters(msClus);
    msrec->setClustersMCLabels(msCluMCLabels);
    msrec->runTracking();
  }
  msrec->closeTracksOutput();

  timer.Stop();
  timer.Print();
}
